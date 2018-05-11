/*
 *  time_driven_static_connection.h
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


#ifndef TIME_DRIVEN_STATIC_CONNECTION_H
#define TIME_DRIVEN_STATIC_CONNECTION_H

/* @BeginDocumentation
   Name: time_driven_static_synapse - Synapse type for time-driven static connections.

   Description:

   time_driven_static_synapse does not support any kind of plasticity. It simply stores
   the parameters target, weight, delay and receiver port for each connection.

   FirstVersion: April 2019

   Author: Jakob Jordan

   Transmits: TimeDrivenSpikeEvent

   SeeAlso: static_synapse
*/

// Includes from nestkernel:
#include "connection.h"

namespace nest
{
template < typename targetidentifierT >
class TimeDrivenStaticConnection : public Connection< targetidentifierT >
{

public:
  // this line determines which common properties to use
  typedef CommonSynapseProperties CommonPropertiesType;

  typedef Connection< targetidentifierT > ConnectionBase;
  typedef TimeDrivenSpikeEvent EventType;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  TimeDrivenStaticConnection()
    : ConnectionBase()
    , weight_( 1.0 )
  {
  }

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not
  // automatically
  // found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  void
  check_connection( Node& s,
    Node& t,
    rport receptor_type,
    const CommonPropertiesType& )
  {
    EventType ge;

    s.sends_secondary_event( ge );
    ge.set_sender( s );
    Connection< targetidentifierT >::target_.set_rport(
      t.handles_test_event( ge, receptor_type ) );
    Connection< targetidentifierT >::target_.set_target( &t );
  }

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param p The port under which this connection is stored in the Connector.
   */
  void
  send( Event& e, thread t, const CommonSynapseProperties& )
  {
    e.set_weight( weight_ );
    e.set_delay( get_delay_steps() );
    e.set_receiver( *get_target( t ) );
    e.set_rport( get_rport() );
    // e();

    // instead of calling handle() in the TimeDrivenSpikeEvent, we
    // loop over all lags and send out a SpikeEvent for every
    // timestep, in which the TimeDrivenSpikeEvent contains an entry;
    // for this we need to copy all information from the original to
    // the new event

    EventType re = *static_cast< EventType* >( &e );

    size_t lag = 0;
    std::vector< unsigned int >::iterator it = re.begin();
    // The call to get_coeffvalue( it ) in this loop also advances the iterator it
    while ( it != re.end() )
    {
      const unsigned int v = re.get_coeffvalue( it );
      if ( v > 0 )
      {
        SpikeEvent se( re, lag );
        se();
      }
      ++lag;
    }
  }

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:
  double weight_; //!< connection weight
};

template < typename targetidentifierT >
void
TimeDrivenStaticConnection< targetidentifierT >::get_status(
  DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
TimeDrivenStaticConnection< targetidentifierT >::set_status(
  const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
}

} // namespace

#endif /* #ifndef TIME_DRIVEN_STATIC_CONNECTION_H */
