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

#ifndef EPROPCONNECTION_H
#define EPROPCONNECTION_H

// Includes from nestkernel:
#include "connection.h"

namespace nest
{

/** @BeginDocumentation
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
template < typename targetidentifierT >
class EPropConnection : public Connection< targetidentifierT >
{
  double weight_;
  Time time_lastspike_;
  Time time_second_to_lastspike_;
  double z_hat_;
  std::vector< double > elg_;
  std::vector< double > times_;
  double tau_m_;

public:
  // this line determines which common properties to use
  typedef CommonSynapseProperties CommonPropertiesType;

  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  EPropConnection()
    : ConnectionBase()
    , weight_( 1.0 )
    , time_lastspike_( Time() )
    , time_second_to_lastspike_( Time() )
    , z_hat_( 0.0 )
    , elg_( std::vector< double >( 0 ) )
    , times_( std::vector< double >( 0 ) )
    , tau_m_( 20. )
  {
  }

  /**
   * Copy constructor from a property object.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  EPropConnection( const EPropConnection& rhs )
    : ConnectionBase( rhs )
    , weight_( rhs.weight_ )
    , time_lastspike_( rhs.time_lastspike_ )
    , time_second_to_lastspike_( rhs.time_second_to_lastspike_ )
    , z_hat_( rhs.z_hat_ )
    , elg_( rhs.elg_ )
    , times_( rhs.times_ )
    , tau_m_( rhs.tau_m_ )
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
  }

  void
  send( Event& e, const thread tid, const CommonSynapseProperties& )
  {
    update_eligibility_trace_( tid );

    e.set_weight( weight_ );
    e.set_delay_steps( get_delay_steps() );
    e.set_receiver( *get_target( tid ) );
    e.set_rport( get_rport() );
    e();

    time_second_to_lastspike_ = time_lastspike_;
    time_lastspike_ = e.get_stamp();
  }

  void
  update_eligibility_trace_( const thread tid )
  {
    if ( time_lastspike_.get_steps() == 0 )
    {
      return;
    }

    const double z_hat = z_hat_ * std::exp( -( time_lastspike_.get_ms() - time_second_to_lastspike_.get_ms() ) / tau_m_ ) + 1.;

    const Archiving_Node* post_node = static_cast< Archiving_Node* >( get_target( tid ) );
    const double pseudo_derivative = post_node->get_pseudo_derivative( time_lastspike_, get_delay_steps() );
    double elg = pseudo_derivative * z_hat;
    if ( elg_.size() > 0 )
    {
      elg += std::exp( -( time_lastspike_.get_ms() - time_second_to_lastspike_.get_ms() ) / tau_m_ ) * elg_.back();
    }
    elg_.push_back( elg );
    times_.push_back( time_lastspike_.get_ms() );

    z_hat_ = z_hat;
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
EPropConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{

  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< long >( d, names::size_of, sizeof( *this ) );

  initialize_property_doublevector( d, "elg" );
  append_property( d, "elg", elg_ );

  initialize_property_doublevector( d, "times" );
  append_property( d, "times", times_ );
}

template < typename targetidentifierT >
void
EPropConnection< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
}

} // namespace

#endif /* #ifndef EPROPCONNECTION_H */
