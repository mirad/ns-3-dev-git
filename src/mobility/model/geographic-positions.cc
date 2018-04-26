/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2018 Bonn-Rhein-Sieg University of Applied Sciences
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
 * Author: Michael Rademacher <michael.rademacher@h-brs.de>
 */

#include <ns3/log.h>
#include <cmath>
#include "geographic-positions.h"

NS_LOG_COMPONENT_DEFINE ("GeographicPositions");

namespace ns3 {

/// Earth's radius in meters if modeled as a perfect sphere
static const double EARTH_RADIUS = 6371e3; 

/**
 * GRS80 and WGS84 sources
 * 
 * Moritz, H. "Geodetic Reference System 1980." GEODETIC REFERENCE SYSTEM 1980. 
 * <http://www.gfy.ku.dk/~iag/HB2000/part4/grs80_corr.htm>.
 * 
 * "Department of Defense World Geodetic System 1984." National Imagery and 
 * Mapping Agency, 1 Jan. 2000. 
 * <http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf>.
 */

/// Earth's semi-major axis in meters as defined by both GRS80 and WGS84
static const double EARTH_SEMIMAJOR_AXIS = 6378137;

/// Earth's first eccentricity as defined by GRS80
static const double EARTH_GRS80_ECCENTRICITY = 0.0818191910428158;

/// Earth's first eccentricity as defined by WGS84
static const double EARTH_WGS84_ECCENTRICITY = 0.0818191908426215;

Vector
GeographicPositions::GeographicToCartesianCoordinates (double latitude, 
                                                       double longitude, 
                                                       double altitude,
                                                       EarthSpheroidType sphType)
{
  NS_LOG_FUNCTION_NOARGS ();
  double latitudeRadians = 0.01745329 * latitude;
  double longitudeRadians = 0.01745329 * longitude;
  double a; // semi-major axis of earth
  double e; // first eccentricity of earth
  if (sphType == SPHERE)
    {
      a = EARTH_RADIUS;
      e = 0;
    }
  else if (sphType == GRS80)
    {
      a = EARTH_SEMIMAJOR_AXIS;
      e = EARTH_GRS80_ECCENTRICITY;
    }
  else // if sphType == WGS84
    {
      a = EARTH_SEMIMAJOR_AXIS;
      e = EARTH_WGS84_ECCENTRICITY;
    }

  double Rn = a / (sqrt (1 - pow (e, 2) * pow (sin (latitudeRadians), 2))); // radius of curvature
  double x = (Rn + altitude) * cos (latitudeRadians) * cos (longitudeRadians);
  double y = (Rn + altitude) * cos (latitudeRadians) * sin (longitudeRadians);
  double z = ((1 - pow (e, 2)) * Rn + altitude) * sin (latitudeRadians);
  Vector cartesianCoordinates = Vector (x, y, z);
  return cartesianCoordinates;
}

Vector
GeographicPositions::CartesianCoordinatesToEnu (double x,
                                                double y,
                                                double z,
                                                double latitude0,
                                                double longitude0,
                                                double altitude0,
                                                EarthSpheroidType sphType0)
{
  NS_LOG_FUNCTION_NOARGS ();
  double latitudeRadians0 = 0.01745329 * latitude0; // lambda
  double longitudeRadians0 = 0.01745329 * longitude0; // phi
  // Get the cartesianCoordinates (ECEF) from the local Geodetic point
  Vector cartesianCoordinates0 = GeographicToCartesianCoordinates(latitude0, longitude0, altitude0, sphType0);

  // Conduct some difference calculuation for better readability of the matrix operation below
  double xd, yd, zd;
  xd = x - cartesianCoordinates0.x;
  yd = y - cartesianCoordinates0.y;
  zd = z - cartesianCoordinates0.z;

  double sin_phi = sin(longitudeRadians0);
  double cos_phi = cos(longitudeRadians0);
  double sin_lambda = sin(latitudeRadians0);
  double cos_lambda = cos(latitudeRadians0);
  // Codunct the transformation using a matrix multiplication.
  double xEast  = -sin_phi                  * xd   + cos_phi              * yd;
  double yNorth = -cos_phi    * sin_lambda  * xd   - sin_lambda * sin_phi * yd   + cos_lambda * zd;
  double yUp    =  cos_lambda * cos_phi     * xd   + cos_lambda * sin_phi * yd   + sin_lambda * zd;

  Vector enuCoordinates = Vector(xEast, yNorth, yUp);
  return enuCoordinates;
}


Vector
GeographicPositions::GeographicToEnu (double latitude,
                                      double longitude,
                                      double altitude,
                                      double latitude0,
                                      double longitude0,
                                      double altitude0,
                                      EarthSpheroidType sphType)
{
  // First to Cartesian, afterwards to ENU
  Vector cartCoord = GeographicToCartesianCoordinates(latitude, longitude, altitude, sphType);
  Vector enuCoordinates = CartesianCoordinatesToEnu(cartCoord.x,cartCoord.y, cartCoord.z, latitude0, longitude0, altitude0, sphType);

  return enuCoordinates;
}

std::list<Vector>
GeographicPositions::RandCartesianPointsAroundGeographicPoint (double originLatitude, 
                                                               double originLongitude, 
                                                               double maxAltitude,
                                                               int numPoints, 
                                                               double maxDistFromOrigin,
                                                               Ptr<UniformRandomVariable> uniRand)
{
  NS_LOG_FUNCTION_NOARGS ();
  // fixes divide by zero case and limits latitude bounds
  if (originLatitude >= 90)
    {
      NS_LOG_WARN ("origin latitude must be less than 90. setting to 89.999");
      originLatitude = 89.999;
    }
  else if (originLatitude <= -90)
    {
      NS_LOG_WARN ("origin latitude must be greater than -90. setting to -89.999");
      originLatitude = -89.999;
    }

  // restricts maximum altitude from being less than zero (below earth's surface).
  // sets maximum altitude equal to zero if parameter is set to be less than zero.
  if (maxAltitude < 0)
    {
      NS_LOG_WARN ("maximum altitude must be greater than or equal to 0. setting to 0");
      maxAltitude = 0;
    }

  double originLatitudeRadians = originLatitude * (M_PI / 180);
  double originLongitudeRadians = originLongitude * (M_PI / 180);
  double originColatitude = (M_PI / 2) - originLatitudeRadians;

  double a = maxDistFromOrigin / EARTH_RADIUS; // maximum alpha allowed 
                                               // (arc length formula)
  if (a > M_PI)
    {
      a = M_PI; // pi is largest alpha possible (polar angle from origin that 
                // points can be generated within)
    }
  
  std::list<Vector> generatedPoints;
  for (int i = 0; i < numPoints; i++)
    {
      // random distance from North Pole (towards center of earth)
      double d = uniRand->GetValue (0, EARTH_RADIUS - EARTH_RADIUS * cos (a)); 
      // random angle in latitude slice (wrt Prime Meridian), radians
      double phi = uniRand->GetValue (0, M_PI * 2); 
      // random angle from Center of Earth (wrt North Pole), radians
      double alpha = acos((EARTH_RADIUS - d) / EARTH_RADIUS); 

      // shift coordinate system from North Pole referred to origin point referred
      // reference: http://en.wikibooks.org/wiki/General_Astronomy/Coordinate_Systems
      double theta = M_PI / 2 - alpha; // angle of elevation of new point wrt 
                                       // origin point (latitude in coordinate 
                                       // system referred to origin point)
      double randPointLatitude = asin(sin(theta)*cos(originColatitude) + 
                                 cos(theta)*sin(originColatitude)*sin(phi)); 
                                 // declination
      double intermedLong = asin((sin(randPointLatitude)*cos(originColatitude) - 
                            sin(theta)) / (cos(randPointLatitude)*sin(originColatitude))); 
                            // right ascension
      intermedLong = intermedLong + M_PI / 2; // shift to longitude 0

      //flip / mirror point if it has phi in quadrant II or III (wasn't 
      //resolved correctly by arcsin) across longitude 0
      if (phi > (M_PI / 2) && phi <= ((3 * M_PI) / 2))
      intermedLong = -intermedLong;

      // shift longitude to be referenced to origin
      double randPointLongitude = intermedLong + originLongitudeRadians; 

      // random altitude above earth's surface
      double randAltitude = uniRand->GetValue (0, maxAltitude);

      Vector pointPosition = GeographicPositions::GeographicToCartesianCoordinates 
                             (randPointLatitude * (180/M_PI), 
                              randPointLongitude * (180/M_PI),
                              randAltitude,
                              SPHERE);
                              // convert coordinates 
                              // from geographic to cartesian

      generatedPoints.push_back (pointPosition); //add generated coordinate 
                                                      //points to list
    }
  return generatedPoints;
}

} // namespace ns3

