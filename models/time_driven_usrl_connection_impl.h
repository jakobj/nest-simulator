/*
 *  time_driven_usrl_connection_impl.h
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

#include "time_driven_usrl_connection.h"

#include "archiving_node.h"

namespace nest
{

template < typename targetidentifierT >
TimeDrivenUSRLConnection< targetidentifierT >::TimeDrivenUSRLConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , E_( 0. )
  , tau_M_( 500. )
  , tau_ex_( 2.0 )
  , Tau_( 10.0 )
  , C_( 250.0 )
  , Theta_( -55. )
  , rho_( 0.01 )
  , delta_( 5. )
  , V_m_( 0. )
  , i_syn_ex_( 0. )
{
}

template < typename targetidentifierT >
void
TimeDrivenUSRLConnection< targetidentifierT >::send( Event& e, thread t, const CommonSynapseProperties& )
{
  e.set_weight( weight_ );
  e.set_delay( get_delay_steps() );
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
    V_m_ = V_m_ * P22_ + i_syn_ex_ * P21ex_;
    i_syn_ex_ *= P11ex_;
    // weighted_spikes_ex_ = spikes_ex_.get_value( lag );
    // i_syn_ex_ += weighted_spikes_ex_;
    i_syn_ex_ += spikes_ex_.get_value( lag );

    const Archiving_Node* post_node = static_cast< Archiving_Node* >( &re.get_receiver() );
    const double s = post_node->get_activity( lag );
    const double u = post_node->get_u( lag );

    const double h = Time::get_resolution().get_ms();

    // update elgibility trace
    E_ = E_ * PE_ + (1. - PE_ ) * 1. / delta_ * V_m_ * ( s - phi_( u ) * h ) * 1e9;

    const unsigned int v = re.get_coeffvalue( it );
    if ( v > 0 )
    {
      spikes_ex_.add_value( re.get_delay() + lag, 1. );
      SpikeEvent se( re, lag );
      se();
    }
    ++lag;
  }
}

} // namespace nest
