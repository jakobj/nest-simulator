/*
 *  rate_volume_transmitter.h
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

#ifndef RATE_VOLUME_TRANSMITTER_H
#define RATE_VOLUME_TRANSMITTER_H

// Includes from models:
#include "rate_neuron_ipn.h"
#include "rate_neuron_ipn_impl.h"
#include "rate_neuron_opn.h"
#include "rate_neuron_opn_impl.h"


namespace nest
{
/* BeginDocumentation
Name: rate_volume_transmitter - Linear rate model (Ornstein-Uhlenbeck process).

Description:

 rate_volume_transmitter is an implementation of a linear rate model with time evolution of an
 Ornstein-Uhlenbeck process.

 Supports connections to other rate_volume_transmitter models with either zero or non-zero
delay.
 Uses the secondary_event concept introduced with the gap-junction framework.

Parameters:

 The following parameters can be set in the status dictionary.

 rate       double - Rate (unitless)
 tau        double - Time constant in ms.
 mean       double - Mean of Gaussian white noise.
 std        double - Standard deviation of Gaussian white noise.

References:

 Hahne, J., Helias, M., Kunkel, S., Igarashi, J.,
 Bolten, M., Frommer, A. and Diesmann, M.,
 A unified framework for spiking and gap-junction interactions
 in distributed neuronal network simulations,
 Front. Neuroinform. 9:22. (2015),
 doi: 10.3389/fninf.2015.00022

Sends: RateNeuronEvent, DelayRateNeuronEvent

Receives: RateNeuronEvent, DelayRateNeuronEvent, DataLoggingRequest

Author: David Dahmen, Jan Hahne, Jannis Schuecker
SeeAlso: rate_connection, delay_rate_connection
*/

class nonlinearities_rate_volume_transmitter
{
private:
  /** gain factor of gain function */
  double g_;

  /** threshold of gain function */
  double theta_;

public:
  /** sets default parameters */
  nonlinearities_rate_volume_transmitter()
    : g_( 1.0 )
    , theta_( 0.0 )
  {
  }

  void get( DictionaryDatum& ) const; //!< Store current values in dictionary
  void set( const DictionaryDatum& ); //!< Set values from dicitonary

  double input( double h );               // non-linearity on input
  double mult_coupling_ex( double rate ); // factor of multiplicative coupling
  double mult_coupling_in( double rate ); // factor of multiplicative coupling
};

inline double nonlinearities_rate_volume_transmitter::input( double h )
{
  return g_ * ( h - theta_ );
}

inline double nonlinearities_rate_volume_transmitter::mult_coupling_ex( double rate )
{
  return rate;
}

inline double nonlinearities_rate_volume_transmitter::mult_coupling_in( double rate )
{
  return rate;
}

typedef rate_neuron_ipn< nest::nonlinearities_rate_volume_transmitter > rate_volume_transmitter_ipn;
typedef rate_neuron_opn< nest::nonlinearities_rate_volume_transmitter > rate_volume_transmitter_opn;

template <>
void RecordablesMap< rate_volume_transmitter_ipn >::create();
template <>
void RecordablesMap< rate_volume_transmitter_opn >::create();

} // namespace nest


#endif /* #ifndef RATE_VOLUME_TRANSMITTER_H */
