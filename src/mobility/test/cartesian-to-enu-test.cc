/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2014 University of Washington
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Michael Rademacher <michael.rademacher@h-brs.de>
 */

#include <ns3/test.h>
#include <ns3/log.h>
#include <ns3/geographic-positions.h>

NS_LOG_COMPONENT_DEFINE ("CartesianToEnuTest");

using namespace ns3;

// 10 meter tolerance for testing

const double TOLERANCE = 10;

const double PY_LATITUDE[20] = {40.4515, -58.33198, 46.21359, -68.30818, -29.38809, -39.48692, -21.62076, -13.11518, 8.16865, 2.0474, -12.93731, -16.61893, -87.24257, 2.98975, 39.06312, -68.89821, 74.71039, 39.7921, -83.9569, -87.7697};
const double PY_LONGITUDE[20] = {40.76942, -98.4833, -137.98293, -106.34823, -46.38852, 21.31267, -173.65946, 99.34394, 25.02673, -52.78749, -130.37588, 179.47591, -84.37122, -102.00597, -60.80065, 75.56805, -73.16124, -45.62193, 25.17815, 11.06719};
const double PY_ALTITUDE[20] = {24.53377, 87.47744, 2.19083, 27.75906, 9.35936, 56.21888, 55.96524, 49.02098, 44.20026, 69.82811, 51.55482, 72.24857, 70.39217, 11.91279, 98.79186, 49.67058, 93.58012, 21.66397, 27.11736, 76.56395};

const double PY_LATITUDE0[1] = {56.89573};
const double PY_LONGITUDE0[1] = {68.40773};
const double PY_ALTITUDE0[1] = {3.58296};

const double PY_X[20] = {-2254660.2317729373, -761311.1881176669, 1965181.203847663, -216092.56587409764, -5049088.257814381, -3610541.269871612, 5241030.169872816, 3193942.525112406, -4336677.322626001, -5452506.072052102, 2001946.9495030446, 5704789.518119741, -140827.5309734901, -1060736.0639813072, -3842484.441901084, 287063.9969980133, -1048710.4341403996, -4482217.878141689, -461436.17500234826, -209673.31085852173};
const double PY_Y[20] = {-1339088.6954630343, -194055.39827353787, 5839531.359554926, -1232690.2324435487, 274185.7734686434, -4994763.814163696, 1071804.4083751882, -5229732.912232256, -3332800.000589893, 2908824.1749714636, 4175499.8211402334, 870594.3830910032, -3218908.7573753223, 5461161.274102118, 4828993.458071813, -5132203.447627224, 4474928.411350561, 3911175.9092143173, -3844036.156907589, -3562227.829910987};
const double PY_Z[20] = {-563377.0811454075, -12676588.68400391, -4688203.491668948, -12594610.537037417, -10243527.52274681, -7909753.591259146, -9837191.083282571, -4657046.728393177, -3102696.5221898155, -7976718.052676314, -10766353.677236367, -9081701.872832824, -11831431.265356543, -9516567.019614203, -4726364.843332922, -10080823.473582584, -1949682.4499756778, -4053378.1962837875, -11390191.705734283, -11610653.56688637};



/**
 * \ingroup mobility
 * \defgroup mobility-test mobility module tests
 */


/**
 * \ingroup mobility-test
 * \ingroup tests
 *
 * \brief Geo To Cartesian Test Case
 *
 * This test verifies the accuracy of the GeographicToEnu
 * method in the GeographicPositions class, which converts earth
 * geographic/geodetic coordinates to ENU coordinates. To do so, it
 * compares the values generated from the method to values generated from the
 * pymap3d function geodetic2enu(), which can be found here https://github.com/scivision/pymap3d
 */
class CartesianToEnuTestCase : public TestCase
{
public:
  /**
   * Constructor
   *
   * \param latitude latitude
   * \param longitude longitude
   * \param altitude altitude
   * \param latitude0 latitude0
   * \param longitude0 longitude0
   * \param altitude0 altitude0
   * \param sphType sphere type
   * \param i index
   */
  CartesianToEnuTestCase (double latitude,
                          double longitude,
                          double altitude,
                          double latitude0,
                          double longitude0,
                          double altitude0,
                          GeographicPositions::EarthSpheroidType sphType,
                          int i);
  virtual ~CartesianToEnuTestCase ();

private:
  virtual void DoRun (void);
  /**
   * Name function
   * \param latitude the latitude
   * \param longitude the longitude
   * \param altitude the altitude
   * \param latitude0 the latitude center
   * \param longitude0 the longitude center
   * \param altitude0 the altitude at the center
   * \param sphType the sphere type
   * \returns the name string
   */
  static std::string Name (double latitude,
                           double longitude,
                           double altitude,
                           double latitude0,
                           double longitude0,
                           double altitude0,
                           GeographicPositions::EarthSpheroidType sphType);
  double m_latitude; ///< latitude
  double m_longitude; ///< longitude
  double m_altitude; ///< altitude
  double m_latitude0; ///< latitude
  double m_longitude0; ///< longitude
  double m_altitude0; ///< altitude
  GeographicPositions::EarthSpheroidType m_sphType; ///< spheroid type
  int m_i; ///< index
};

std::string 
CartesianToEnuTestCase::Name (double latitude,
                              double longitude,
                              double altitude,
                              double latitude0,
                              double longitude0,
                              double altitude0,
                              GeographicPositions::EarthSpheroidType sphType)
{
  std::ostringstream oss;
  oss << "latitude = " << latitude << " degrees, "
      << "longitude = " << longitude << " degrees, "
      << "altitude = " << altitude << " meters, "
      << "latitude = " << latitude0 << " degrees, "
      << "longitude = " << longitude0 << " degrees, "
      << "altitude = " << altitude0 << " degrees, "
      << "earth spheroid type = " << sphType;
  return oss.str();
}

CartesianToEnuTestCase::CartesianToEnuTestCase (double latitude,
                                                double longitude,
                                                double altitude,
                                                double latitude0,
                                                double longitude0,
                                                double altitude0,
                                                GeographicPositions::EarthSpheroidType sphType,
                                                int i)
  : TestCase (Name (latitude, longitude, altitude, latitude0, longitude0, altitude0, sphType)),
    m_latitude (latitude),
    m_longitude (longitude),
    m_altitude (altitude),
    m_latitude0 (latitude0),
    m_longitude0 (longitude0),
    m_altitude0 (altitude0),
    m_sphType (sphType),
    m_i (i)
{
}

CartesianToEnuTestCase::~CartesianToEnuTestCase ()
{
}

void
CartesianToEnuTestCase::DoRun (void)
{  
 Vector cart = GeographicPositions::GeographicToEnu (m_latitude,
                                                     m_longitude,
                                                     m_altitude,
                                                     m_latitude0,
                                                     m_longitude0,
                                                     m_altitude0,
                                                     m_sphType);

  NS_TEST_ASSERT_MSG_EQ_TOL (cart.x,
                             PY_X [m_i],
                             TOLERANCE,
                             "x coordinate (" << cart.x << ") is incorrect "
                             "for WGS84 model in iteration " << m_i);
  NS_TEST_ASSERT_MSG_EQ_TOL (cart.y,
                             PY_Y [m_i],
                             TOLERANCE,
                             "y coordinate (" << cart.y << ") is incorrect "
                             "for WGS84 model in iteration " << m_i);
  NS_TEST_ASSERT_MSG_EQ_TOL (cart.z,
                             PY_Z [m_i],
                             TOLERANCE,
                             "z coordinate (" << cart.z << ") is incorrect "
                             "for WGS84 model in iteration " << m_i);


}

/**
 * \ingroup mobility-test
 * \ingroup tests
 *
 * \brief Cartesian To Enu Test Suite
 */
class CartesianToEnuTestSuite : public TestSuite
{
public:
  CartesianToEnuTestSuite ();
};

CartesianToEnuTestSuite::CartesianToEnuTestSuite ()
  : TestSuite ("cartesian-to-enu", UNIT)
{
  NS_LOG_INFO ("creating CartesianToEnuTestSuite");
  for (int i = 0; i < 20; i += 1)
    {
      std::cout << i <<std::endl;
      AddTestCase (new CartesianToEnuTestCase (PY_LATITUDE [i],
                                               PY_LONGITUDE [i],
                                               PY_ALTITUDE [i],
                                               PY_LATITUDE0 [0],
                                               PY_LONGITUDE0 [0],
                                               PY_ALTITUDE0 [0],
                                               GeographicPositions::WGS84,
                                               i),
                   TestCase::QUICK);
      //++i;
      std::cout << i << std::endl;
    }
}

static CartesianToEnuTestSuite g_CartesianToEnuTestSuite;
