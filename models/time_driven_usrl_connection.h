/*
 *  time_driven_usrl_connection.h
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


#ifndef TIME_DRIVEN_USRL_CONNECTION_H
#define TIME_DRIVEN_USRL_CONNECTION_H

/* @BeginDocumentation TODO fix
   Name: time_driven_usrl_synapse - Synapse type for time-driven usrl connections.

   Description:

   time_driven_usrl_synapse does not support any kind of plasticity. It simply stores
   the parameters target, weight, delay and receiver port for each connection.

   FirstVersion: April 2019

   Author: Jakob Jordan

   Transmits: TimeDrivenSpikeEvent

   SeeAlso: time_driven_static_synapse
*/

// Includes from nestkernel:
#include "connection.h"

// Includes from libnestutil:
#include "numerics.h"
#include "propagator_stability.h"

// Includes from nestkernel:
#include "ring_buffer.h"

namespace nest
{
template < typename targetidentifierT >
class TimeDrivenUSRLConnection : public Connection< targetidentifierT >
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
  TimeDrivenUSRLConnection()
    : ConnectionBase()
    , weight_( 1.0 )
    , E_( 0. )
    , tau_M_( 500. )
    , tau_( 2.0 )
    , Tau_( 10.0 )
    , C_( 250.0 )
    , Theta_( -55. )
    , rho_( 0.01 )
    , delta_( 5. )
    , PSP_( 0. )
    , i_syn_( 0. )
    , initialized_( false )
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
  void send( Event& e, thread t, const CommonSynapseProperties& )
  {
    if ( not initialized_ )
    {
      init_buffers();
      init_propagators();
      initialized_ = true;
    }

    e.set_weight( weight_ );
    e.set_delay_steps( get_delay_steps() );
    e.set_receiver( *get_target( t ) );
    e.set_rport( get_rport() );
    // e();

    EventType re = *static_cast< EventType* >( &e );

    size_t lag = 0;
    std::vector< unsigned int >::iterator it = re.begin();
    // The call to get_coeffvalue( it ) in this loop also advances the iterator it
    while ( it != re.end() )
      {

        // propagate PSP
        PSP_ = PSP_ * P22_ + i_syn_ * P21_;
        i_syn_ *= P11_;
        // weighted_spikes_ = spikes_.get_value( lag );
        // i_syn_ += weighted_spikes_;
        i_syn_ += spikes_.get_value( lag );

        const Archiving_Node* post_node = static_cast< Archiving_Node* >( &re.get_receiver() );
        const double s = post_node->get_activity( lag );
        const double u = post_node->get_u( lag );

        const double h = Time::get_resolution().get_ms();

        // update elgibility trace
        E_ = E_ * PE_ + (1. - PE_ ) * 1. / delta_ * PSP_ * ( s - phi_( u ) * h ) * 1e9;  // TODO where does the factor 1e9 come from?

        const unsigned int v = re.get_coeffvalue( it );
        if ( v > 0 )
          {
            spikes_.add_value( re.get_delay_steps() + lag, 1. );
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

  void init_buffers();

  void init_propagators();

private:
  double weight_; //!< connection weight
  double E_; //!< eligibility trace
  double tau_M_; //!< time constant of eligibility trace

  // propagator for eligibility trace
  double PE_;

  // parameters of postsynaptic neuron
  double tau_;
  double Tau_;
  double C_;
  double Theta_;
  double rho_;
  double delta_;

  // propagators, required for computing PSP
  double P11_;
  double P21_;
  double P22_;

  // dynamic variables, required for computing PSP
  double PSP_;
  double i_syn_;
  // double weighted_spikes_;

  bool initialized_;

  // buffers to store incoming spikes and process them after a delay
  RingBuffer spikes_;

  double phi_( const double u ) const;
};

template < typename targetidentifierT >
void
TimeDrivenUSRLConnection< targetidentifierT >::get_status(
  DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< long >( d, names::size_of, sizeof( *this ) );
  def< double >( d, "E", E_ );
  def< double >( d, "PSP", PSP_ );
}

template < typename targetidentifierT >
void
TimeDrivenUSRLConnection< targetidentifierT >::set_status(
  const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  // updateValue< double >( d, "E", E_ );
  // updateValue< double >( d, "PSP", PSP_ );
  // updateValue< double >( d, "i_syn", i_syn_ );
}

template < typename targetidentifierT >
void
TimeDrivenUSRLConnection< targetidentifierT >::init_buffers()
{
  spikes_.clear(); // includes resize
}

template < typename targetidentifierT >
void
TimeDrivenUSRLConnection< targetidentifierT >::init_propagators()
{
  const double h = Time::get_resolution().get_ms();
  P11_ = std::exp( -h / tau_ );
  P22_ = std::exp( -h / Tau_ );
  P21_ = propagator_32( tau_, Tau_, C_, h );

  PE_ = std::exp( -h / tau_M_ );
}

template < typename targetidentifierT >
inline double
TimeDrivenUSRLConnection< targetidentifierT >::phi_( const double u ) const
{
  return rho_ * std::exp( 1 / delta_ * ( u - Theta_ ) );
}

} // namespace

#endif /* #ifndef TIME_DRIVEN_USRL_CONNECTION_H */
