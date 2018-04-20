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
 * Author: Benjamin Cizdziel <ben.cizdziel@gmail.com>
 */

#include <ns3/vector.h>
#include <ns3/random-variable-stream.h>

#ifndef GEOGRAPHIC_POSITIONS_H
#define GEOGRAPHIC_POSITIONS_H

namespace ns3
{

/**
 * \ingroup mobility
 *
 * Consists of methods dealing with Earth geographic coordinates and locations.
 */
class GeographicPositions
{

public:
  /** 
   * Spheroid model to use for earth: perfect sphere (SPHERE), Geodetic 
   * Reference System 1980 (GRS80), or World Geodetic System 1984 (WGS84)
   * 
   * Moritz, H. "Geodetic Reference System 1980." GEODETIC REFERENCE SYSTEM 1980. 
   * <http://www.gfy.ku.dk/~iag/HB2000/part4/grs80_corr.htm>.
   * 
   * "Department of Defense World Geodetic System 1984." National Imagery and 
   * Mapping Agency, 1 Jan. 2000. 
   * <http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf>.
   */
  enum EarthSpheroidType
  {
    SPHERE,
    GRS80,
    WGS84
  };

  /**
   * Converts earth geographic/geodetic coordinates (latitude and longitude in 
   * degrees) with a given altitude above earth's surface (in meters) to Earth 
   * Centered Earth Fixed (ECEF) Cartesian coordinates (x, y, z in meters), 
   * where origin (0, 0, 0) is the center of the earth.
   *
   * @param latitude earth-referenced latitude (in degrees) of the point
   * @param longitude earth-referenced longitude (in degrees) of the point
   * @param altitude height of the point (in meters) above earth's surface
   * @param sphType earth spheroid model to use for conversion
   *
   * @return a vector containing the Cartesian coordinates (x, y, z referenced 
   * in meters) of the point (origin (0, 0, 0) is center of earth)
   */
  static Vector GeographicToCartesianCoordinates (double latitude, 
                                                  double longitude, 
                                                  double altitude,
                                                  EarthSpheroidType sphType);
  
  /** Converts the Earth-Centered Earth-Fixed (ECEF) coordinates
  * or called Cartesian coordinates in ns3 (x, y, z)
  * to East-North-Up coordinates in a Local Tangent Plane that is centered at the
  * Geodetic point (latitude0, longitude0, altitude0)
  * Source: https://gist.github.com/govert/1b373696c9a27ff4c72a
  *
  * @param Cartesian coordinate x of the regared point
  * @param Cartesian coordinate y of the regared point
  * @param Cartesian coordinate z of the regared point
  * @param reference latitude0 (in degrees) of the center-point
  * @param reference longitude0 (in degrees) of the center-point
  * @param reference altitude0 (meter) of the center-point above earth's surface
  * @param sphType earth spheroid model to use for conversion
  *
  * @return a Vector containing the ENU coordinates (x, y, z referenced in meters)
  * of the regared point, centered at the center-point
  *
  */
 static Vector CartesianCoordinatesToEnu (double x,
                                          double y,
                                          double z,
                                          double latitude0,
                                          double longitude0,
                                          double altitude0,
                                          EarthSpheroidType sphType0);


  /** Converts earth geographic/geodetic coordinates (latitude and longitude in
   * degrees) with a given altitude above earth's surface (in meters) to
   * East-North-Up coordinates in a local tangent plane that is centered at the
   * Geodetic point (latitude0, longitude0, altitude0)
   * Source: https://gist.github.com/govert/1b373696c9a27ff4c72a
   *
   * @param latitude (in degrees)
   * @param longitude (in degrees)
   * @param Cartesian coordinate z
   * @param reference latitude0 (in degrees) of the center-point
   * @param reference longitude0 (in degrees) of the center-point
   * @param reference altitude0 (meter) of the center-point, above earth's surface
   * @param sphType earth spheroid model to use for conversion
   *
   *
   * @return a Vector containing the ENU coordinates (x, y, z referenced in meters)
   * of the regared point, centered at the center-point
   */
 static Vector GeographicToEnu (double latitude,
                                double longitude,
                                double altitude,
                                double latitude0,
                                double longitude0,
                                double altitude0,
                                EarthSpheroidType sphType0);

  /**
   * Generates uniformly distributed random points (in ECEF Cartesian 
   * coordinates) within a given altitude above earth's surface centered around 
   * a given origin point (on earth's surface, in geographic/geodetic coordinates) 
   * within a given distance radius (using arc length of earth's surface, not 
   * pythagorean distance).
   * Distance radius is measured as if all generated points are on earth's 
   * surface (with altitude = 0).
   * Assumes earth is a perfect sphere.
   * 
   * @param originLatitude origin point latitude in degrees
   * @param originLongitude origin point longitude in degrees
   * @param maxAltitude maximum altitude in meters above earth's surface with
   * which random points can be generated
   * @param numPoints number of points to generate
   * @param maxDistFromOrigin max distance in meters from origin with which 
   * random transmitters can be generated (all transmitters are less than or 
   * equal to this distance from the origin, relative to points being on earth's 
   * surface)
   * @param uniRand pointer to the uniform random variable to use for random 
   * location and altitude generation
   *
   * @return a list containing the vectors (x, y, z location referenced in 
   * meters from origin at center of earth) of each point generated
   */
  static std::list<Vector> RandCartesianPointsAroundGeographicPoint (double originLatitude, 
                                                                     double originLongitude,
                                                                     double maxAltitude, 
                                                                     int numPoints, 
                                                                     double maxDistFromOrigin,
                                                                     Ptr<UniformRandomVariable> uniRand);

};

} // namespace ns3

#endif /* GEOGRAPHIC_POSITIONS_H */

