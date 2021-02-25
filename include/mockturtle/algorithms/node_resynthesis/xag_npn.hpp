/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file xag_npn.hpp
  \brief Replace with size-optimum XAGs and AIGs from NPN (from ABC rewrite)

  \author Heinz Riener
  \author Mathias Soeken
  \author Siang-Yun (Sonia) Lee
*/

#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include <fmt/format.h>
#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/npn.hpp>
#include <kitty/print.hpp>
#include <kitty/static_truth_table.hpp>

#include "../../algorithms/simulation.hpp"
#include "../../networks/xag.hpp"
#include "../../utils/index_list.hpp"
#include "../../utils/node_map.hpp"
#include "../../utils/stopwatch.hpp"

namespace mockturtle
{

struct xag_npn_resynthesis_params
{
  /*! \brief Be verbose. */
  bool verbose{false};
};

struct xag_npn_resynthesis_stats
{
  stopwatch<>::duration time_classes{0};
  stopwatch<>::duration time_db{0};

  uint32_t db_size;
  uint32_t covered_classes;

  void report() const
  {
    std::cout << fmt::format( "[i] build classes time = {:>5.2f} secs\n", to_seconds( time_classes ) );
    std::cout << fmt::format( "[i] build db time      = {:>5.2f} secs\n", to_seconds( time_db ) );
  }
};

/*! \brief Resynthesis function based on pre-computed AIGs.
 *
 * This resynthesis function can be passed to ``cut_rewriting``.  It will
 * produce a network based on pre-computed XAGs with up to at most 4 variables.
 * Consequently, the nodes' fan-in sizes in the input network must not exceed
 * 4.
 *
   \verbatim embed:rst

   Example

   .. code-block:: c++

      const aig_network aig = ...;
      xag_npn_resynthesis<aig_network> resyn;
      aig = cut_rewriting( aig, resyn );

   .. note::

      The implementation of this algorithm was heavily inspired by the rewrite
      command in AIG.  It uses the same underlying database of subcircuits.
   \endverbatim
 */
template<class Ntk, class DatabaseNtk = xag_network>
class xag_npn_resynthesis
{
public:
  xag_npn_resynthesis( xag_npn_resynthesis_params const& ps = {}, xag_npn_resynthesis_stats* pst = nullptr )
      : ps( ps ),
        pst( pst ),
        _repr( 1u << 16u )
  {
    static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
    static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
    static_assert( has_create_and_v<Ntk>, "Ntk does not implement the create_and method" );
    static_assert( has_create_xor_v<Ntk>, "Ntk does not implement the create_xor method" );
    static_assert( has_create_not_v<Ntk>, "Ntk does not implement the create_not method" );

    static_assert( is_network_type_v<DatabaseNtk>, "DatabaseNtk is not a network type" );
    static_assert( has_get_node_v<DatabaseNtk>, "DatabaseNtk does not implement the get_node method" );
    static_assert( has_is_complemented_v<DatabaseNtk>, "DatabaseNtk does not implement the is_complemented method" );
    static_assert( has_is_xor_v<DatabaseNtk>, "DatabaseNtk does not implement the is_xor method" );
    static_assert( has_size_v<DatabaseNtk>, "DatabaseNtk does not implement the size method" );
    static_assert( has_create_pi_v<DatabaseNtk>, "DatabaseNtk does not implement the create_pi method" );
    static_assert( has_create_and_v<DatabaseNtk>, "DatabaseNtk does not implement the create_and method" );
    static_assert( has_create_xor_v<DatabaseNtk>, "DatabaseNtk does not implement the create_xor method" );
    static_assert( has_foreach_fanin_v<DatabaseNtk>, "DatabaseNtk does not implement the foreach_fanin method" );
    static_assert( has_foreach_node_v<DatabaseNtk>, "DatabaseNtk does not implement the foreach_node method" );
    static_assert( has_make_signal_v<DatabaseNtk>, "DatabaseNtk does not implement the make_signal method" );

    build_classes();
    build_db();
  }

  virtual ~xag_npn_resynthesis()
  {
    if ( ps.verbose )
    {
      st.report();
    }

    if ( pst )
    {
      *pst = st;
    }
  }

  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& ntk, kitty::dynamic_truth_table const& function, LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
    kitty::static_truth_table<4u> tt = kitty::extend_to<4u>( function );

    /* get representative of function */
    const auto [repr, phase, perm] = _repr[*tt.cbegin()];

    /* check if representative has circuits */
    const auto it = _repr_to_signal.find( repr );
    if ( it == _repr_to_signal.end() )
    {
      return;
    }

    std::vector<signal<Ntk>> pis( 4, ntk.get_constant( false ) );
    std::copy( begin, end, pis.begin() );

    std::unordered_map<node<DatabaseNtk>, signal<Ntk>> db_to_ntk;
    db_to_ntk.insert( {0, ntk.get_constant( false )} );
    for ( auto i = 0; i < 4; ++i )
    {
      db_to_ntk.insert( {i + 1, ( phase >> perm[i] & 1 ) ? ntk.create_not( pis[perm[i]] ) : pis[perm[i]]} );
    }

    for ( auto const& cand : it->second )
    {
      const auto f = copy_db_entry( ntk, _db.get_node( cand ), db_to_ntk );
      if ( !fn( _db.is_complemented( cand ) != ( phase >> 4 & 1 ) ? ntk.create_not( f ) : f ) )
      {
        return;
      }
    }
  }

private:
  signal<Ntk>
  copy_db_entry( Ntk& ntk, node<DatabaseNtk> const& n, std::unordered_map<node<DatabaseNtk>, signal<Ntk>>& db_to_ntk ) const
  {
    if ( const auto it = db_to_ntk.find( n ); it != db_to_ntk.end() )
    {
      return it->second;
    }

    std::array<signal<Ntk>, 2> fanin{};
    _db.foreach_fanin( n, [&]( auto const& f, auto i ) {
      const auto ntk_f = copy_db_entry( ntk, _db.get_node( f ), db_to_ntk );
      fanin[i] = _db.is_complemented( f ) ? ntk.create_not( ntk_f ) : ntk_f;
    } );

    const auto f = _db.is_xor( n ) ? ntk.create_xor( fanin[0], fanin[1] ) : ntk.create_and( fanin[0], fanin[1] );
    db_to_ntk.insert( {n, f} );
    return f;
  }

  void build_classes()
  {
    stopwatch t( st.time_classes );

    kitty::static_truth_table<4u> tt;
    do {
      _repr[*tt.cbegin()] = kitty::exact_npn_canonization( tt );
      kitty::next_inplace( tt );
    } while ( !kitty::is_const0( tt ) );
  }

  void build_db()
  {
    stopwatch t( st.time_db );

    decode( _db, xag_index_list{std::vector<uint32_t>{subgraphs, subgraphs + sizeof subgraphs / sizeof subgraphs[0]}} );
    const auto sim_res = simulate_nodes<kitty::static_truth_table<4u>>( _db );

    _db.foreach_node( [&]( auto n ) {
      if ( std::get<0>( _repr[*sim_res[n].cbegin()] ) == sim_res[n] )
      {
        if ( _repr_to_signal.count( sim_res[n] ) == 0 )
        {
          _repr_to_signal.insert( {sim_res[n], {_db.make_signal( n )}} );
        }
        else
        {
          _repr_to_signal[sim_res[n]].push_back( _db.make_signal( n ) );
        }
      }
      else
      {
        const auto f = ~sim_res[n];
        if ( std::get<0>( _repr[*f.cbegin()] ) == f )
        {
          if ( _repr_to_signal.count( f ) == 0 )
          {
            _repr_to_signal.insert( {f, {!_db.make_signal( n )}} );
          }
          else
          {
            _repr_to_signal[f].push_back( !_db.make_signal( n ) );
          }
        }
      }
    } );

    st.db_size = _db.size();
    st.covered_classes = static_cast<uint32_t>( _repr_to_signal.size() );
  }

  xag_npn_resynthesis_params ps;
  xag_npn_resynthesis_stats st;
  xag_npn_resynthesis_stats* pst{nullptr};

  std::vector<std::tuple<kitty::static_truth_table<4u>, uint32_t, std::vector<uint8_t>>> _repr;
  std::unordered_map<kitty::static_truth_table<4u>, std::vector<signal<DatabaseNtk>>, kitty::hash<kitty::static_truth_table<4u>>> _repr_to_signal;

  DatabaseNtk _db;

  // clang-format off
  inline static const uint32_t subgraphs[] = {1780 << 16 | 0 << 8 | 4,
      2, 4, 2, 5, 3, 4, 3, 5, 4, 2, 2, 6, 2, 7, 3, 6, 3, 7, 6, 2, 4, 6,
      4, 7, 5, 6, 5, 7, 6, 4, 2, 8, 2, 9, 3, 8, 3, 9, 8, 2, 4, 8,
      4, 9, 5, 8, 5, 9, 8, 4, 6, 8, 6, 9, 7, 8, 7, 9, 8, 6, 5, 11,
      6, 10, 6, 11, 7, 10, 7, 11, 10, 6, 8, 10, 8, 11, 9, 10, 9, 11, 10, 8,
      6, 12, 6, 13, 7, 12, 7, 13, 12, 6, 9, 12, 9, 13, 12, 8, 5, 15, 6, 14,
      6, 15, 7, 14, 7, 15, 14, 6, 8, 14, 8, 15, 9, 14, 9, 15, 14, 8, 6, 16,
      6, 17, 7, 16, 7, 17, 16, 6, 8, 16, 8, 17, 9, 16, 9, 17, 16, 8, 6, 18,
      6, 19, 7, 18, 7, 19, 18, 6, 8, 19, 9, 18, 9, 19, 18, 8, 4, 20, 4, 21,
      5, 20, 7, 21, 8, 20, 9, 21, 20, 8, 11, 21, 20, 10, 15, 21, 20, 14, 17, 21,
      19, 21, 4, 22, 4, 23, 5, 22, 9, 22, 9, 23, 22, 8, 22, 12, 15, 23, 17, 23,
      18, 23, 4, 24, 7, 25, 9, 25, 24, 8, 11, 25, 13, 25, 15, 25, 24, 14, 19, 25,
      4, 26, 4, 27, 5, 26, 5, 27, 26, 4, 8, 27, 9, 26, 9, 27, 26, 8, 11, 27,
      13, 27, 17, 27, 26, 16, 19, 27, 4, 28, 28, 4, 9, 28, 9, 29, 28, 8, 11, 28,
      11, 29, 13, 29, 17, 29, 18, 29, 19, 28, 19, 29, 2, 30, 2, 31, 3, 30, 5, 31,
      7, 31, 8, 30, 8, 31, 9, 30, 9, 31, 30, 8, 13, 31, 17, 31, 19, 31, 23, 31,
      27, 31, 29, 31, 2, 32, 2, 33, 5, 33, 32, 6, 8, 33, 32, 8, 13, 33, 17, 33,
      21, 33, 25, 33, 27, 33, 28, 33, 32, 28, 2, 34, 3, 35, 34, 4, 7, 35, 34, 8,
      11, 35, 15, 35, 19, 35, 34, 18, 23, 35, 27, 35, 33, 35, 2, 36, 2, 37, 3, 36,
      3, 37, 36, 2, 8, 36, 8, 37, 9, 36, 9, 37, 36, 8, 11, 37, 15, 37, 17, 37,
      18, 37, 19, 37, 21, 37, 25, 37, 27, 37, 29, 37, 2, 38, 3, 38, 38, 2, 8, 38,
      8, 39, 9, 38, 9, 39, 38, 8, 11, 38, 11, 39, 15, 38, 15, 39, 17, 39, 18, 38,
      19, 38, 19, 39, 21, 39, 23, 38, 25, 39, 27, 38, 27, 39, 28, 38, 29, 38, 29, 39,
      4, 40, 4, 41, 6, 40, 9, 41, 13, 41, 15, 41, 19, 41, 23, 41, 25, 41, 29, 41,
      31, 41, 33, 41, 35, 41, 36, 41, 37, 41, 40, 36, 39, 41, 4, 42, 4, 43, 5, 42,
      6, 43, 7, 42, 17, 43, 27, 43, 30, 43, 31, 42, 31, 43, 32, 43, 33, 42, 42, 32,
      36, 43, 37, 42, 37, 43, 42, 36, 39, 42, 39, 43, 42, 38, 7, 45, 9, 45, 11, 45,
      21, 45, 31, 45, 44, 32, 36, 44, 36, 45, 39, 45, 44, 38, 4, 46, 4, 47, 5, 46,
      6, 47, 7, 46, 46, 6, 13, 47, 19, 47, 23, 47, 31, 46, 46, 30, 32, 47, 33, 47,
      34, 47, 35, 47, 36, 46, 36, 47, 37, 46, 37, 47, 46, 36, 38, 47, 39, 47, 4, 49,
      48, 4, 6, 49, 48, 6, 15, 48, 19, 48, 19, 49, 25, 48, 28, 49, 29, 48, 29, 49,
      31, 49, 33, 48, 35, 48, 36, 49, 39, 48, 48, 38, 2, 50, 2, 51, 6, 50, 7, 51,
      9, 51, 13, 51, 19, 51, 21, 51, 23, 51, 25, 51, 26, 51, 50, 26, 31, 51, 35, 51,
      39, 51, 47, 51, 48, 51, 2, 53, 3, 52, 6, 52, 6, 53, 17, 53, 22, 52, 23, 52,
      23, 53, 26, 53, 27, 53, 37, 53, 45, 53, 3, 55, 7, 55, 9, 55, 11, 55, 21, 55,
      22, 55, 23, 55, 26, 54, 26, 55, 31, 55, 43, 55, 53, 55, 2, 56, 3, 56, 6, 57,
      7, 56, 56, 6, 11, 57, 15, 57, 19, 57, 21, 56, 56, 20, 23, 57, 24, 57, 25, 56,
      25, 57, 26, 56, 27, 56, 27, 57, 56, 26, 33, 57, 41, 57, 2, 59, 3, 59, 58, 2,
      6, 59, 7, 58, 7, 59, 58, 6, 13, 59, 17, 59, 19, 59, 58, 20, 25, 59, 26, 59,
      27, 59, 58, 28, 35, 58, 58, 34, 38, 58, 38, 59, 39, 58, 43, 59, 47, 59, 2, 60,
      4, 60, 4, 61, 5, 61, 60, 4, 9, 61, 10, 61, 11, 61, 13, 61, 15, 61, 16, 61,
      17, 61, 18, 61, 19, 61, 23, 61, 27, 61, 33, 61, 39, 61, 43, 61, 47, 61, 48, 61,
      60, 52, 57, 61, 58, 61, 2, 63, 4, 62, 4, 63, 12, 63, 13, 62, 17, 63, 19, 63,
      27, 63, 37, 63, 45, 63, 55, 63, 3, 65, 5, 65, 9, 65, 11, 65, 16, 64, 16, 65,
      18, 65, 21, 65, 31, 65, 43, 65, 53, 65, 57, 65, 63, 65, 2, 66, 2, 67, 3, 66,
      3, 67, 66, 2, 4, 67, 5, 66, 66, 4, 10, 67, 11, 66, 66, 10, 13, 67, 14, 67,
      15, 67, 16, 66, 17, 66, 17, 67, 66, 16, 18, 66, 19, 66, 19, 67, 66, 18, 25, 67,
      35, 67, 41, 67, 51, 67, 57, 67, 3, 69, 68, 2, 4, 68, 4, 69, 5, 68, 68, 4,
      11, 69, 16, 69, 17, 68, 17, 69, 68, 16, 18, 68, 18, 69, 68, 18, 23, 69, 27, 69,
      68, 32, 37, 69, 39, 68, 43, 69, 47, 69, 57, 69, 58, 68, 70, 68, 9, 73, 37, 73,
      41, 73, 45, 73, 51, 73, 55, 73, 61, 73, 65, 73, 74, 2, 74, 4, 74, 16, 74, 18,
      33, 75, 41, 75, 74, 46, 48, 75, 51, 75, 58, 75, 67, 75, 8, 77, 9, 77, 76, 8,
      17, 77, 35, 77, 51, 77, 61, 77, 68, 77, 69, 77, 76, 68, 78, 2, 8, 79, 9, 78,
      9, 79, 78, 8, 17, 79, 78, 16, 31, 79, 63, 79, 9, 80, 17, 80, 48, 81, 51, 81,
      6, 83, 7, 83, 17, 83, 68, 83, 75, 83, 7, 85, 21, 85, 6, 87, 17, 87, 27, 87,
      37, 87, 61, 87, 65, 87, 67, 87, 69, 87, 6, 89, 7, 88, 88, 6, 23, 89, 88, 22,
      25, 88, 33, 89, 61, 89, 90, 74, 9, 93, 15, 93, 9, 94, 9, 95, 15, 95, 94, 14,
      23, 95, 9, 97, 25, 97, 31, 97, 51, 97, 9, 99, 98, 14, 21, 99, 35, 99, 45, 98,
      47, 99, 98, 46, 69, 99, 7, 103, 33, 103, 35, 103, 37, 103, 39, 103, 51, 103, 61, 103,
      6, 105, 21, 104, 27, 104, 27, 105, 104, 26, 31, 104, 37, 105, 104, 36, 55, 105, 65, 105,
      69, 105, 108, 68, 112, 4, 13, 113, 33, 113, 58, 112, 61, 115, 97, 115, 103, 115, 116, 4,
      8, 117, 9, 117, 57, 117, 116, 56, 63, 117, 13, 118, 13, 119, 55, 119, 61, 121, 98, 121,
      122, 2, 23, 125, 61, 125, 6, 127, 126, 6, 21, 126, 27, 127, 37, 126, 37, 127, 126, 36,
      65, 127, 77, 127, 126, 78, 128, 20, 77, 131, 9, 133, 11, 133, 132, 10, 79, 133, 88, 133,
      8, 135, 9, 134, 9, 135, 134, 8, 11, 135, 43, 135, 134, 42, 53, 135, 61, 135, 63, 135,
      73, 135, 87, 135, 134, 86, 134, 88, 136, 4, 136, 6, 9, 136, 9, 137, 136, 8, 47, 137,
      53, 136, 55, 136, 57, 137, 61, 137, 63, 137, 69, 137, 136, 68, 75, 137, 89, 137, 105, 137,
      127, 137, 9, 139, 138, 8, 11, 138, 11, 139, 88, 139, 21, 141, 31, 141, 43, 141, 63, 141,
      73, 141, 133, 141, 7, 143, 25, 143, 27, 143, 35, 143, 39, 143, 47, 143, 67, 143, 75, 143,
      78, 143, 79, 143, 95, 143, 97, 143, 101, 143, 131, 143, 6, 145, 7, 144, 7, 145, 144, 6,
      73, 145, 144, 78, 143, 145, 146, 4, 7, 146, 146, 6, 27, 146, 27, 147, 69, 147, 135, 147,
      7, 148, 148, 6, 25, 149, 61, 149, 81, 149, 101, 149, 131, 149, 150, 90, 9, 153, 143, 153,
      9, 154, 154, 62, 61, 157, 68, 157, 156, 68, 9, 159, 158, 8, 61, 158, 68, 158, 83, 158,
      143, 159, 149, 159, 68, 161, 151, 161, 158, 161, 7, 162, 21, 164, 28, 164, 31, 164, 38, 164,
      164, 134, 166, 6, 68, 167, 166, 72, 158, 167, 9, 169, 168, 166, 170, 48, 174, 58, 158, 177,
      15, 178, 19, 178, 180, 14, 9, 186, 188, 8, 9, 192, 192, 68, 68, 195, 194, 68, 104, 197,
      202, 78, 9, 211, 210, 8, 216, 58, 41, 226, 48, 226, 143, 229, 149, 229, 41, 230, 51, 230,
      61, 231, 178, 235, 9, 236, 236, 42, 236, 86, 236, 164, 238, 58, 178, 241, 19, 243, 159, 243,
      186, 243, 192, 243, 211, 243, 213, 243, 229, 243, 5, 244, 19, 249, 149, 249, 159, 249, 186, 249,
      192, 249, 211, 249, 143, 257, 149, 257, 243, 257, 249, 257, 68, 261, 260, 180, 9, 263, 262, 8,
      61, 262, 68, 262, 83, 262, 143, 263, 149, 263, 177, 262, 243, 263, 249, 263, 19, 264, 264, 36,
      38, 264, 268, 4, 143, 271, 149, 271, 272, 58, 9, 280, 9, 283, 41, 282, 48, 282, 51, 282,
      58, 282, 61, 282, 68, 282, 282, 166, 9, 285, 286, 58, 290, 68, 292, 58, 158, 295, 262, 295,
      296, 38, 19, 300, 300, 36, 97, 300, 300, 134, 199, 300, 223, 300, 300, 236, 252, 300, 302, 36,
      302, 222, 246, 304, 9, 308, 308, 68, 243, 308, 249, 308, 310, 58, 104, 312, 314, 58, 68, 317,
      316, 68, 126, 319, 320, 68, 322, 8, 326, 6, 326, 210, 330, 204, 246, 331, 330, 248, 332, 48,
      126, 335, 334, 128, 143, 339, 149, 339, 41, 341, 48, 341, 346, 8, 348, 58, 350, 4, 350, 252,
      178, 353, 352, 180, 41, 354, 41, 356, 51, 356, 61, 357, 41, 359, 48, 359, 104, 361, 360, 106,
      362, 106, 364, 8, 300, 367, 300, 369, 9, 370, 370, 42, 370, 86, 370, 300, 41, 373, 48, 373,
      372, 48, 41, 374, 48, 374, 300, 375, 376, 104, 376, 178, 376, 264, 376, 300, 19, 379, 21, 379,
      31, 379, 48, 379, 101, 379, 103, 379, 159, 379, 182, 379, 185, 379, 192, 379, 207, 379, 263, 379,
      267, 379, 271, 379, 277, 379, 301, 379, 308, 379, 339, 379, 369, 379, 374, 379, 3, 380, 379, 381,
      186, 382, 211, 382, 19, 385, 31, 385, 48, 385, 93, 385, 101, 385, 103, 385, 153, 385, 159, 385,
      173, 385, 185, 385, 192, 385, 207, 385, 384, 210, 384, 254, 263, 385, 267, 385, 271, 385, 279, 385,
      308, 385, 339, 385, 343, 385, 374, 385, 41, 387, 178, 388, 388, 220, 390, 248, 243, 393, 126, 396,
      398, 128, 400, 8, 400, 148, 379, 403, 385, 403, 404, 302, 379, 405, 385, 405, 143, 407, 149, 407,
      243, 407, 249, 407, 9, 409, 408, 8, 61, 408, 68, 408, 83, 408, 143, 409, 149, 409, 177, 408,
      243, 409, 249, 409, 295, 408, 379, 409, 385, 409, 19, 414, 28, 414, 414, 134, 414, 236, 414, 370,
      413, 415, 418, 2, 143, 418, 243, 418, 422, 48, 51, 425, 41, 426, 48, 426, 9, 429, 379, 429,
      411, 429, 418, 429, 385, 431, 9, 432, 9, 435, 41, 434, 48, 434, 51, 434, 58, 434, 61, 434,
      68, 434, 436, 48, 41, 440, 48, 440, 51, 443, 9, 445, 411, 445, 418, 445, 9, 446, 385, 449,
      9, 451, 41, 450, 48, 450, 51, 450, 58, 450, 61, 450, 68, 450, 68, 453, 158, 453, 262, 453,
      408, 453, 454, 270, 454, 338, 454, 344, 158, 457, 262, 457, 408, 457, 458, 134, 458, 154, 458, 236,
      458, 370, 334, 460, 373, 460, 397, 460, 35, 462, 39, 462, 239, 462, 313, 462, 360, 462, 25, 464,
      29, 464, 35, 464, 39, 464, 464, 78, 354, 466, 387, 466, 19, 468, 39, 468, 307, 468, 328, 468,
      385, 468, 19, 470, 39, 470, 379, 472, 385, 472, 75, 474, 224, 474, 251, 474, 15, 476, 19, 476,
      171, 476, 191, 476, 208, 476, 478, 218, 478, 300, 15, 484, 19, 484, 25, 484, 29, 484, 486, 134,
      486, 236, 486, 370, 224, 489, 239, 489, 274, 489, 313, 489, 360, 489, 379, 491, 385, 491, 307, 493,
      328, 493, 35, 497, 39, 497, 39, 499, 385, 499, 379, 503, 385, 503, 61, 511, 512, 298, 512, 414,
      514, 416, 13, 517, 287, 517, 367, 517, 518, 298, 379, 521, 51, 525, 61, 525, 385, 525, 528, 134,
      528, 154, 528, 236, 528, 370, 517, 535, 517, 541, 379, 543, 501, 545, 7, 550, 51, 553, 462, 553,
      464, 553, 497, 553, 507, 553, 549, 553, 5, 554, 51, 557, 385, 557, 462, 557, 464, 557, 497, 557,
      564, 478, 385, 567, 462, 573, 497, 573, 578, 300, 580, 300, 472, 584, 51, 587, 61, 587, 143, 588,
      243, 588, 68, 591, 545, 591, 592, 6, 68, 593, 545, 595, 596, 4, 39, 598, 600, 154, 151, 603,
      39, 604, 537, 611, 379, 613, 75, 614, 171, 616, 379, 619, 15, 620, 25, 620, 51, 620, 61, 620,
      61, 623, 259, 623, 395, 623, 561, 623, 68, 625, 158, 625, 262, 625, 408, 625, 626, 420, 626, 438,
      158, 629, 262, 629, 408, 629, 630, 46, 630, 56, 630, 162, 632, 134, 632, 154, 632, 236, 632, 370,
      636, 78, 553, 636, 557, 636, 574, 640, 535, 642, 243, 648, 67, 650, 77, 650, 195, 650, 317, 650,
      533, 650, 562, 650, 650, 608, 652, 42, 48, 652, 103, 652, 201, 652, 491, 652, 495, 652, 313, 654,
      360, 654, 39, 656, 334, 659, 373, 659, 397, 659, 662, 512, 33, 665, 191, 665, 208, 665, 569, 665,
      614, 665, 670, 36, 670, 218, 15, 677, 171, 677, 570, 677, 67, 679, 77, 679, 195, 679, 317, 679,
      533, 679, 562, 679, 669, 679, 334, 681, 397, 681, 243, 683, 484, 683, 665, 685, 686, 134, 686, 154,
      686, 236, 686, 370, 476, 691, 677, 691, 694, 522, 61, 701, 468, 701, 499, 701, 700, 606, 313, 703,
      360, 703, 704, 6, 379, 707, 385, 707, 652, 707, 7, 708, 51, 711, 143, 711, 464, 711, 636, 711,
      3, 712, 51, 715, 143, 715, 149, 715, 464, 715, 636, 715, 623, 721, 722, 644, 379, 725, 517, 726,
      143, 729, 379, 731, 474, 732, 736, 630, 584, 739, 583, 740, 638, 740, 61, 742, 468, 742, 499, 742,
      334, 744, 397, 744, 249, 749, 750, 6, 68, 751, 754, 46, 753, 755, 468, 757, 499, 757, 557, 757,
      758, 2, 39, 763, 151, 765, 61, 766, 61, 769, 243, 771, 39, 773, 61, 775, 259, 775, 395, 775,
      561, 775, 721, 775, 776, 508, 776, 564, 647, 781, 782, 532, 784, 468, 39, 789, 557, 789, 158, 791,
      262, 791, 408, 791, 158, 793, 262, 793, 408, 793, 796, 66, 796, 178, 796, 264, 798, 504, 800, 134,
      800, 154, 800, 236, 800, 370, 802, 68, 804, 66, 135, 804, 237, 804, 371, 804, 804, 416, 577, 804,
      623, 804, 735, 804, 775, 804, 570, 806, 691, 806, 808, 98, 535, 808, 808, 684, 726, 808, 810, 66,
      814, 66, 689, 818, 718, 818, 820, 568, 822, 42, 57, 822, 822, 502, 531, 822, 558, 822, 822, 674,
      737, 822, 824, 42, 48, 824, 103, 824, 201, 824, 491, 824, 495, 824, 707, 824, 826, 32, 307, 828,
      328, 828, 39, 830, 832, 570, 33, 834, 191, 834, 208, 834, 836, 508, 354, 839, 387, 839, 840, 512,
      35, 843, 224, 843, 251, 843, 616, 843, 642, 845, 846, 126, 75, 849, 574, 849, 732, 849, 151, 851,
      689, 853, 718, 853, 57, 855, 531, 855, 558, 855, 737, 855, 354, 857, 387, 857, 35, 859, 224, 859,
      251, 859, 143, 861, 484, 861, 843, 863, 864, 134, 864, 154, 864, 236, 864, 370, 474, 867, 640, 867,
      849, 867, 143, 871, 483, 873, 647, 873, 851, 873, 634, 875, 667, 875, 51, 877, 462, 877, 497, 877,
      307, 879, 328, 879, 33, 881, 191, 881, 208, 881, 882, 570, 884, 4, 652, 887, 824, 887, 808, 889,
      816, 889, 5, 890, 892, 480, 808, 895, 816, 895, 61, 897, 243, 897, 816, 897, 3, 898, 61, 901,
      243, 901, 468, 901, 499, 901, 816, 901, 904, 810, 904, 814, 61, 907, 517, 908, 642, 908, 243, 911,
      476, 912, 677, 912, 804, 915, 650, 917, 679, 917, 916, 796, 916, 810, 61, 919, 249, 919, 557, 919,
      584, 919, 740, 919, 61, 920, 583, 920, 739, 920, 804, 920, 804, 921, 924, 814, 61, 927, 61, 928,
      634, 930, 667, 930, 51, 932, 462, 932, 497, 932, 354, 934, 387, 934, 35, 936, 224, 936, 251, 936,
      938, 582, 149, 941, 919, 941, 942, 4, 947, 949, 950, 2, 763, 953, 769, 953, 789, 953, 919, 953,
      143, 955, 39, 959, 161, 965, 157, 969, 195, 969, 763, 969, 769, 969, 789, 969, 957, 969, 39, 971,
      591, 971, 157, 975, 157, 979, 763, 979, 769, 979, 789, 979, 957, 979, 39, 981, 591, 981, 157, 983,
      988, 78, 135, 988, 237, 988, 371, 988, 407, 988, 988, 538, 988, 696, 988, 868, 41, 991, 48, 991,
      379, 992, 379, 995, 243, 996, 243, 999, 143, 1000, 143, 1003, 1004, 488, 1006, 658, 83, 1009, 143, 1011,
      243, 1011, 379, 1011, 41, 1012, 48, 1012, 33, 1014, 685, 1014, 385, 1017, 789, 1017, 33, 1018, 527, 1020,
      1022, 46, 143, 1024, 1026, 6, 1028, 962, 1030, 132, 9, 1033, 1032, 8, 243, 1033, 249, 1033, 379, 1033,
      385, 1033, 51, 1034, 35, 1036, 645, 1038, 957, 1043, 143, 1046, 379, 1046, 67, 1049, 283, 1049, 435, 1049,
      451, 1049, 1050, 810, 133, 1052, 1052, 810, 73, 1055, 1054, 132, 763, 1055, 769, 1055, 789, 1055, 919, 1055,
      9, 1057, 1056, 8, 9, 1058, 143, 1059, 51, 1061, 645, 1063, 1064, 154, 161, 1067, 1068, 796, 1070, 178,
      1072, 18, 1072, 716, 1074, 88, 1074, 126, 1076, 74, 1080, 16, 61, 1083, 1084, 812, 137, 1087, 553, 1089,
      557, 1089, 711, 1089, 715, 1089, 33, 1091, 685, 1091, 1090, 854, 23, 1093, 527, 1093, 1092, 852, 1094, 156,
      137, 1097, 1098, 582, 957, 1101, 143, 1103, 761, 1103, 787, 1103, 1104, 810, 83, 1107, 149, 1107, 919, 1107,
      379, 1108, 1110, 794, 379, 1113, 243, 1114, 1114, 672, 137, 1116, 1120, 36, 243, 1122, 1124, 126, 1126, 124,
      1126, 546, 1126, 660, 1128, 58, 9, 1131, 1132, 58, 31, 1134, 143, 1136, 379, 1136, 9, 1138, 493, 1140,
      1142, 796, 61, 1145, 143, 1147, 379, 1147, 9, 1149, 61, 1151, 489, 1153, 489, 1154, 897, 1154, 945, 1157,
      1158, 794, 665, 1160, 843, 1162, 1164, 8, 51, 1166, 61, 1166, 39, 1168, 39, 1170, 507, 1170, 462, 1173,
      497, 1173, 1174, 796, 31, 1176, 197, 1179, 197, 1180, 897, 1180, 27, 1182, 381, 1185, 39, 1189, 39, 1191,
      507, 1191, 945, 1193, 41, 1197, 48, 1197, 1198, 48, 1200, 48, 761, 1203, 787, 1203, 747, 1204, 1206, 8,
      379, 1208, 51, 1210, 1212, 336, 1214, 154, 659, 1217, 659, 1218, 889, 1218, 747, 1221, 9, 1223, 61, 1225,
      61, 1227, 1228, 98, 61, 1231, 39, 1232, 1234, 324, 693, 1236, 143, 1239, 157, 1239, 77, 1241, 83, 1241,
      195, 1241, 317, 1241, 453, 1241, 625, 1241, 37, 1242, 51, 1245, 121, 1245, 21, 1246, 319, 1249, 319, 1250,
      889, 1250, 693, 1253, 1254, 6, 9, 1261, 1260, 8, 11, 1262, 79, 1262, 1262, 134, 9, 1264, 161, 1267,
      9, 1268, 1270, 134, 31, 1273, 47, 1273, 73, 1273, 111, 1273, 133, 1273, 145, 1273, 169, 1273, 211, 1273,
      215, 1273, 285, 1273, 289, 1273, 565, 1273, 915, 1273, 804, 1275, 73, 1276, 169, 1276, 285, 1276, 31, 1279,
      41, 1279, 51, 1279, 61, 1279, 133, 1279, 143, 1279, 211, 1279, 61, 1280, 301, 1283, 179, 1287, 11, 1288,
      903, 1288, 89, 1291, 9, 1292, 67, 1295, 133, 1298, 9, 1301, 1302, 8, 650, 1305, 679, 1305, 1304, 796,
      1306, 814, 1308, 6, 261, 1309, 591, 1309, 665, 1310, 493, 1313, 665, 1315, 493, 1316, 89, 1318, 11, 1321,
      903, 1321, 11, 1323, 9, 1327, 61, 1328, 261, 1331, 405, 1333, 11, 1334, 83, 1339, 9, 1340, 143, 1341,
      665, 1345, 493, 1347, 39, 1349, 75, 1351, 61, 1353, 1354, 802, 1356, 56, 1356, 144, 97, 1358, 199, 1358,
      301, 1358, 367, 1358, 381, 1358, 713, 1358, 899, 1358, 1360, 16, 179, 1362, 47, 1364, 39, 1366, 75, 1368,
      67, 1370, 1372, 144, 1374, 16, 25, 1378, 245, 1382, 555, 1382, 891, 1382, 143, 1385, 497, 1385, 804, 1387,
      1388, 816, 143, 1391, 497, 1391, 379, 1392, 143, 1394, 1396, 78, 39, 1399, 650, 1401, 679, 1401, 1400, 796,
      39, 1403, 1404, 796, 233, 1407, 699, 1407, 11, 1409, 245, 1411, 555, 1411, 891, 1411, 61, 1413};
  // clang-format off
}; // namespace mockturtle

} /* namespace mockturtle */