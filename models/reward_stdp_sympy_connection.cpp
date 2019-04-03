#include "reward_stdp_sympy_connection.h"

// Includes from nestkernel:
#include "kernel_manager.h"

// Includes from sli:
#include "dictdatum.h"

namespace nest
{

RewardSTDPSympyConnectionCommonProperties::RewardSTDPSympyConnectionCommonProperties()
  : vt_( nullptr )
{
}

void
RewardSTDPSympyConnectionCommonProperties::get_status( DictionaryDatum& d ) const
{
  CommonSynapseProperties::get_status( d );
  if ( vt_ != nullptr )
  {
    def< long >( d, names::vt, vt_->get_gid() );
  }
  else
  {
    def< long >( d, names::vt, -1 );
  }
}

void
RewardSTDPSympyConnectionCommonProperties::set_status( const DictionaryDatum& d,
  ConnectorModel& cm )
{
  long vtgid;
  if ( updateValue< long >( d, names::vt, vtgid ) )
  {
    vt_ = dynamic_cast< iaf_psc_delta* >( kernel().node_manager.get_node(
      vtgid, kernel().vp_manager.get_thread_id() ) );
    if ( vt_ == 0 )
    {
      throw BadProperty( "Dopamine source must be iaf psc delta" );
    }
  }
}

} // namespace nest
