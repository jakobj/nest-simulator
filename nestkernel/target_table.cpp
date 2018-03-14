/*
 *  target_table.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "target_table.h"

// Includes from nestkernel:
#include "kernel_manager.h"

void
nest::TargetTable::initialize()
{
  const thread num_threads = kernel().vp_manager.get_num_threads();
  targets_.resize( num_threads, NULL );
  secondary_send_buffer_pos_.resize( num_threads, NULL );

#pragma omp parallel
  {
    const thread tid = kernel().vp_manager.get_thread_id();
    targets_[ tid ] =
      new std::vector< std::vector< Target > >( 0, std::vector< Target >( 0, Target() ) );
    secondary_send_buffer_pos_[ tid ] =
      new std::vector< std::vector< std::vector< size_t > > > ( 0 );
  } // of omp parallel
}

void
nest::TargetTable::finalize()
{
  for ( std::vector< std::vector< std::vector< Target > >* >::iterator it =
          targets_.begin();
        it != targets_.end();
        ++it )
  {
    ( *it )->clear();
    delete *it;
  }
  targets_.clear();

  for ( std::vector< std::vector< std::vector< std::vector< size_t > > >* >::iterator it =
          secondary_send_buffer_pos_.begin();
        it != secondary_send_buffer_pos_.end();
        ++it )
  {
    ( *it )->clear();
    delete *it;
  }
  secondary_send_buffer_pos_.clear();
}

void
nest::TargetTable::prepare( const thread tid )
{
  // add one to max_num_local_nodes to avoid possible overflow in case
  // of rounding errors
  const size_t num_local_nodes = kernel().node_manager.get_max_num_local_nodes() + 1;

  targets_[ tid ]->resize( num_local_nodes,
    std::vector< Target >( 0, Target() ) );

  ( *secondary_send_buffer_pos_[ tid ] ).resize( num_local_nodes );

  for ( size_t lid = 0; lid < num_local_nodes; ++lid )
  {
    // resize to maximal possible synapse-type index
    ( *secondary_send_buffer_pos_[ tid ] )[ lid ].resize( kernel().model_manager.get_num_synapse_prototypes() );
  }
}

void
nest::TargetTable::compress_secondary_send_buffer_pos( const thread tid )
{
  for ( std::vector< std::vector< std::vector< size_t > > >::iterator it = ( *secondary_send_buffer_pos_[ tid ] ).begin();
        it != ( *secondary_send_buffer_pos_[ tid ] ).end(); ++it )
  {
    for ( std::vector< std::vector< size_t > >::iterator iit = ( *it ).begin(); iit != ( *it ).end(); ++iit )
    {
      std::sort( iit->begin(), iit->end() );
      const std::vector< size_t >::iterator new_it = std::unique( iit->begin(), iit->end() );
      iit->resize( std::distance( iit->begin(), new_it) );
    }
  }
}

void
nest::TargetTable::add_target( const thread tid, const thread target_rank, const TargetData& target_data )
{
  const index lid = target_data.get_source_lid();

  // use 1.5 growth strategy (see, e.g.,
  // https://github.com/facebook/folly/blob/master/folly/docs/FBVector.md)
  if ( ( *targets_[ tid ] )[ lid ].size() == ( *targets_[ tid ] )[ lid ].capacity() )
  {
    ( *targets_[ tid ] )[ lid ].reserve( ( ( *targets_[ tid ] )[ lid ].size() * 3 + 1 ) / 2 );
  }

  if ( target_data.is_primary() )
  {
    const TargetDataFields& target_fields = target_data.target_data;

    ( *targets_[ tid ] )[ lid ].push_back(
      Target( target_fields.get_tid(),
              target_rank,
              target_fields.get_syn_id(),
              target_fields.get_lcid() ) );
  }
  else
  {
    const SecondaryTargetDataFields& secondary_fields =
        target_data.secondary_data;
    const size_t send_buffer_pos = secondary_fields.get_recv_buffer_pos()
      + kernel().mpi_manager.get_send_displacement_secondary_events_in_int( target_rank );
    const synindex syn_id =
      secondary_fields.get_syn_id();

    assert( syn_id < ( *secondary_send_buffer_pos_[ tid ] )[ lid ].size() );
    ( *secondary_send_buffer_pos_[ tid ] )[ lid ][ syn_id ].push_back(
      send_buffer_pos );
  }
}

void
nest::TargetTable::print_targets( const thread tid ) const
{
  std::cout << "-------------TARGETS-------------------\n";
  for ( std::vector< std::vector< Target > >::const_iterator cit =
          ( *targets_[ tid ] ).begin();
        cit != ( *targets_[ tid ] ).end();
        ++cit )
  {
    for ( std::vector< Target >::const_iterator ciit = ( *cit ).begin();
          ciit != ( *cit ).end();
          ++ciit )
    {
      std::cout << ( *ciit ).get_lcid() << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "---------------------------------------\n";
}

void
nest::TargetTable::print_secondary_send_buffer_pos( const thread tid ) const
{
  std::cout << "-------------SENDBUFFERPOS-------------------\n";
  for ( std::vector< std::vector< std::vector< size_t > > >::const_iterator cit =
          ( *secondary_send_buffer_pos_[ tid ] ).begin();
        cit != ( *secondary_send_buffer_pos_[ tid ] ).end();
        ++cit )
  {
    for ( std::vector< std::vector< size_t > >::const_iterator ciit = ( *cit ).begin();
          ciit != ( *cit ).end();
          ++ciit )
    {
      for ( std::vector< size_t >::const_iterator ciiit =  ( *ciit ).begin();
            ciiit != ( *ciit ).end();
          ++ciiit )
      {
        std::cout << *ciiit << ", ";
      }
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "---------------------------------------\n";
}
