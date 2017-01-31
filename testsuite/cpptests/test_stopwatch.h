/*
 *  test_stopwatch.h
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

#ifndef TEST_STOPWATCH_H
#define TEST_STOPWATCH_H

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// C++ includes:
#include <time.h>
#include <unistd.h>

// Includes from libnestutil:
#include "stopwatch.h"

namespace nest
{

BOOST_AUTO_TEST_SUITE( test_stopwatch )

BOOST_AUTO_TEST_CASE( test_continue )
{
  Stopwatch sw;
  sw.start();
  usleep( 30000 );
  sw.stop();
  sw.start();
  usleep( 30000 );
  sw.stop();
  sw.start();
  usleep( 30000 );
  sw.stop();
  BOOST_CHECK( sw.elapsed( sw.MILLISEC ) > 60. );
}

BOOST_AUTO_TEST_CASE( test_isrunning )
{
  Stopwatch sw;
  sw.start();
  BOOST_CHECK( sw.isRunning() == true );
  sw.stop();
}

BOOST_AUTO_TEST_CASE( test_reset )
{
  Stopwatch sw;
  sw.start();
  usleep( 30000 );
  sw.stop();
  sw.reset();
  BOOST_CHECK( sw.elapsed() == 0. );
}

BOOST_AUTO_TEST_SUITE_END()

} // of namespace nest

#endif /* TEST_STOPWATCH_H */
