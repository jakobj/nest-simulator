/*
 *  static_connection.h
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


/* BeginDocumentation
  Name: static_synapse - Synapse type for static connections.

  Description:
   static_synapse does not support any kind of plasticity. It simply stores
   the parameters target, weight, delay and receiver port for each connection.

  FirstVersion: October 2005
  Author: Jochen Martin Eppler, Moritz Helias

  Transmits: SpikeEvent, RateEvent, CurrentEvent, ConductanceEvent,
  DoubleDataEvent, DataLoggingRequest

  Remarks: Refactored for new connection system design, March 2007

  SeeAlso: synapsedict, tsodyks_synapse, stdp_synapse
*/

#ifndef DIRECTSTATICCONNECTION_H
#define DIRECTSTATICCONNECTION_H

// Includes from nestkernel:
#include "connection.h"

namespace nest
{

/**
 * Class representing a static connection. A static connection has the
 * properties weight, delay and receiver port. A suitable Connector containing
 * these connections can be obtained from the template GenericConnector.
 */


template < typename targetidentifierT >
class DirectStaticConnection : public Connection< targetidentifierT >
{
  double weight_;
  RingBuffer* target_buffer_;

public:
  // this line determines which common properties to use
  typedef CommonSynapseProperties CommonPropertiesType;

  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  DirectStaticConnection()
    : ConnectionBase()
    , weight_( 1.0 )
    , target_buffer_( NULL )
  {
  }

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  DirectStaticConnection( const DirectStaticConnection& rhs )
    : ConnectionBase( rhs )
    , weight_( rhs.weight_ )
    , target_buffer_( rhs.target_buffer_ )
  {
  }

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;


  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( RateEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( DataLoggingRequest&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( CurrentEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( ConductanceEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( DoubleDataEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( DSSpikeEvent&, rport )
    {
      return invalid_port_;
    }
    port
    handles_test_event( DSCurrentEvent&, rport )
    {
      return invalid_port_;
    }
  };

  void
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
    target_buffer_ = t.get_target_buffer( weight_ );
  }

  void
  send( Event& e, const thread tid, const CommonSynapseProperties& )
  {
    e.set_weight( weight_ );
    e.set_delay( get_delay_steps() );
    e.set_receiver( *get_target( tid ) );
    e.set_rport( get_rport() );
    // e();
    target_buffer_->add_value(
      e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() );
  }

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  void
  set_weight( double w )
  {
    weight_ = w;
  }
};

template < typename targetidentifierT >
void
DirectStaticConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
DirectStaticConnection< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
}

} // namespace

#endif /* #ifndef STATICCONNECTION_H */
