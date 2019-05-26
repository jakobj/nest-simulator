/*
 *  hebbian_rate_connection.h
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
Name: hebbian_rate_connection.h - Synapse type for rate connections with delay.

Description:
 hebbian_rate_connection.h is a connector to create connections with delay
 between rate model neurons.

 To create instantaneous rate connections please use
 the synapse type rate_connection.

Transmits: DelayedRateConnectionEvent

References:

 Hahne, J., Dahmen, D., Schuecker, J., Frommer, A.,
 Bolten, M., Helias, M. and Diesmann, M. (2016).
 Integration of continuous-time dynamics in a spiking
 neural network simulator arXiv:1610.09990 [q-bio.NC]

Author: David Dahmen, Jan Hahne, Jannis Schuecker
SeeAlso: rate_connection, rate_neuron_ipn, rate_neuron_opn
*/


#ifndef HEBBIAN_RATE_CONNECTION_H
#define HEBBIAN_RATE_CONNECTION_H

#include "connection.h"

namespace nest
{

/**
 * Class representing a delay-rate-connection. A delay-rate-connection
 * has the properties weight, delay and receiver port.
 */

template < typename targetidentifierT >
class HebbianRateConnection : public Connection< targetidentifierT >
{

  double A_;
  double weight_;
  double Wmax_;
  double Wmin_;
  rate_volume_transmitter_ipn* vt_; 
  double n_threshold_;
  double post_threshold_;
  double weight_decay_constant_;
  std::vector<double> pre_rate_hist_;
  unsigned int trace_delay_steps_;

public:
  // this line determines which common properties to use
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;
  typedef DelayedRateConnectionEvent EventType;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  HebbianRateConnection()
    : ConnectionBase()
    , A_( 0.01 )
    , weight_( 0.5 )
    , Wmax_( 1.0 )
    , Wmin_( 0.0 )
    , vt_( nullptr )
    , n_threshold_( 1.0 )
    , post_threshold_( 0.0 )
    , weight_decay_constant_( 0. )
    , pre_rate_hist_( 0 )
    , trace_delay_steps_( 1 )
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
    assert( vt_ != nullptr );

//    double t_now = e.get_stamp().get_ms();
    DelayedRateConnectionEvent* re = static_cast< DelayedRateConnectionEvent* >( &e );
    
    const Archiving_Node* const target = static_cast<Archiving_Node*>(get_target( t ));

    const std::vector< double >& post_hist = target->get_rate_hist();

    std::vector< unsigned int >::iterator it = re->begin();
    while ( it != re->end() )
    {
        const double pre_now = re->get_coeffvalue(it);
        pre_rate_hist_.push_back(pre_now);
        if (pre_rate_hist_.size() > trace_delay_steps_){
            pre_rate_hist_.erase(pre_rate_hist_.begin());
        }

        const double post = post_hist.front();
        const double pre = pre_rate_hist_.front();
        //weight_ += A_ * pre * std::pow(target->get_rate(), 2) * ( vt_->get_rate() - n_threshold_) - 0.0001 * weight_;
        // weight_ += A_ * pre * ( target->get_rate() - post_threshold_ ) * ( vt_->get_rate() - n_threshold_) - weight_decay_constant_ * weight_;
        // weight_ += A_ * pre * ( vt_->get_rate() - n_threshold_) - weight_decay_constant_ * weight_;
        //if ( target->get_rate() - post_threshold_  > 0.)
        if ( post - post_threshold_  > 0.)
        {
          weight_ += A_ * pre * ( 1. ) * ( vt_->get_rate() - n_threshold_) - weight_decay_constant_ * weight_;
        }
        // else
        // {
        //   // "homeostasis"
        //   weight_ += A_ * pre * ( -0.1 ) * ( vt_->get_rate() - n_threshold_) - weight_decay_constant_ * weight_;
        // }
    }

    if (weight_ > Wmax_){
        weight_ = Wmax_;
    }
    if (weight_ < Wmin_){
        weight_ = Wmin_;
    }

    e.set_weight( weight_ );
    e.set_delay_steps( get_delay_steps() );
    e.set_receiver( *get_target( t ) );
    e.set_rport( get_rport() );
    e();
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
HebbianRateConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );

  if ( vt_ != nullptr )
  {
    def< long >( d, "vt", vt_->get_gid() );
  }
  else
  {
    def< long >( d, "vt", -1 );
  }

  def< double >( d, "A", A_);
  def< double >( d, names::weight, weight_ );
  def< double >( d, "Wmax", Wmax_);
  def< double >( d, "Wmin", Wmin_);
  def< double >( d, "n_threshold", n_threshold_);
  def< double >( d, "post_threshold", post_threshold_);
  def< double >( d, "weight_decay_constant", weight_decay_constant_);
  def< double >( d, "trace_delay_steps", trace_delay_steps_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
HebbianRateConnection< targetidentifierT >::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );

  long vtgid;
  if ( updateValue< long >( d, "vt", vtgid ) )
  {

    vt_ = dynamic_cast< rate_volume_transmitter_ipn* >( kernel().node_manager.get_node(
      vtgid, kernel().vp_manager.get_thread_id() ) );

    if ( vt_ == nullptr )
      throw BadProperty( "Dopamine source must be of type rate_volume_transmitter_ipn" );
  }

  updateValue< double >( d, "A", A_);
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, "Wmax", Wmax_);
  updateValue< double >( d, "Wmin", Wmin_);
  updateValue< double >( d, "n_threshold", n_threshold_);
  updateValue< double >( d, "post_threshold", post_threshold_);
  updateValue< double >( d, "weight_decay_constant", weight_decay_constant_);
  updateValue< double >( d, "trace_delay_steps", trace_delay_steps_ );
}

} // namespace

#endif /* #ifndef HEBBIAN_RATE_CONNECTION_H */
