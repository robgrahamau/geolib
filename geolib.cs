//**************************************************************************************
//                              COPYRIGHT AND LICENSE INFORMATION
//**************************************************************************************
//The following Library is released by Robert Graham for public use.
//**************************************************************************************
//                           END OF COPYRIGHT AND LICENSE INFORMATION 
//**************************************************************************************
/* GeoLibrary and other data Class. 
Copyright 2016 Robert Graham & Peter Story
Last Date Modified: 11/03/17.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace geolibrary
{
    /// <summary>
    /// Contains Latitude, Longitude, Altitude
    /// </summary>
    public struct GeoLocation
    {
        public double Latitude { get; set; }
        public double Longitude { get; set; }
        public double altitude { get; set; }
    }
    /// <summary>
    /// Class Contains Converters for converting to and from unix time standards.
    /// </summary>
    public static class timeconvert
    {
        /// <summary>
        /// Converts from Unix Epoch time to standard Date Time.
        /// </summary>
        /// <param name="unixTime"></param>
        /// <returns></returns>
        public static DateTime FromUnixTime(this long unixTime)
        {
            var epoch = new DateTime(1970, 1, 1, 0, 0, 0, DateTimeKind.Utc);
            return epoch.AddSeconds(unixTime);
        }
        /// <summary>
        /// Converts to Unix epoch time
        /// </summary>
        /// <param name="date"></param>
        /// <returns></returns>
        public static long ToUnixTime(this DateTime date)
        {
            var epoch = new DateTime(1970, 1, 1, 0, 0, 0, DateTimeKind.Utc);
            return Convert.ToInt64((date - epoch).TotalSeconds);
        }
    }

    /// <summary>
    /// Geographical Math Library contains the following functions
    /// FindPointatdistancefrom - Returns a Geolocation when you give it data
    /// distance - gives you the distance between two lat/lons or Geolocations
    /// bearing - gives you the bearing from one lat long pair (or Geolocation) to another.
    /// DegreesToRadians - Converts Degrees to Radians
    /// RadiantsToDegrees - Converts Radians to Degrees
    /// KMtoNm - Kilometers to Nautical Miles.
    /// NmtoKM - Nautical Miles to Kilometers.
    /// MtoFeet - Meters to Feet
    /// FeettoM - Feet to Meters
    /// </summary>
    public class geolib
    {

         /// <summary>
        /// Using distance in Kms, finds the point in dgerees at a bearing in radians from starting point in degrees.
        /// </summary>
        /// <param name="startPoint"></param>
        /// <param name="initialBearingRadians"></param>
        /// <param name="distanceKilometres"></param>
        /// <returns></returns>
        public static GeoLocation FindPointAtDistanceFrom(GeoLocation startPoint, double initialBearingRadians, double distanceKilometres)
        {
            const double radiusEarthKilometres = 6371.01;
            var distRatio = distanceKilometres / radiusEarthKilometres;
            var distRatioSine = Math.Sin(distRatio);
            var distRatioCosine = Math.Cos(distRatio);

            var startLatRad = DegreesToRadians(startPoint.Latitude);
            var startLonRad = DegreesToRadians(startPoint.Longitude);

            var startLatCos = Math.Cos(startLatRad);
            var startLatSin = Math.Sin(startLatRad);

            var endLatRads = Math.Asin((startLatSin * distRatioCosine) + (startLatCos * distRatioSine * Math.Cos(initialBearingRadians)));

            var endLonRads = startLonRad
                + Math.Atan2(
                    Math.Sin(initialBearingRadians) * distRatioSine * startLatCos,
                    distRatioCosine - startLatSin * Math.Sin(endLatRads));

            return new GeoLocation
            {
                Latitude = RadiansToDegrees(endLatRads),
                Longitude = RadiansToDegrees(endLonRads)
            };
        }

        /// <summary>
        /// Returns the distance between 2 lat/long's in either Miles (M) Kilometers (K) nm (N) or (R) Radians
        /// </summary>
        /// <param name="lat1">Latitude 1</param>
        /// <param name="lon1">Longitude 1</param>
        /// <param name="lat2">Latitude 2</param>
        /// <param name="lon2">Longitude 2</param>
        /// <param name="unit">M/N or K</param>
        /// <returns>Distance in UNIT</returns>
         public static double getdistance(double lat1, double lon1, double lat2, double lon2, char unit)
        {
            double dist_return;
            double theta = lon2 - lon1; //Changed from lon1 - lon2 as is not correct as per http://www.movable-type.co.uk/scripts/latlong.html
            double dist = Math.Sin(DegreesToRadians(lat1)) * Math.Sin(DegreesToRadians(lat2)) + Math.Cos(DegreesToRadians(lat1)) * Math.Cos(DegreesToRadians(lat2)) * Math.Cos(DegreesToRadians(theta));
            dist = Math.Acos(dist); //radians
            double dist_deg = RadiansToDegrees(dist);
            double dist_miles = dist_deg * 60 * 1.1515;
            switch (unit)
            {
                case 'M':
                case 'm':
                    dist_return = dist_miles;
                    break;
                case 'K':
                case 'k':
                    dist_return = dist_miles * 1.609344;
                    break;
                case 'N':
                case 'n':
                    dist_return = dist_miles * 0.8684;
                    break;
                case 'R':
                case 'r':
                    dist_return = dist;
                    break;

                default:
                    Debug.WriteLine("ERROR: Unknown unit specifier: " + unit, "getdistance");
                    dist_return = 0;
                    break;
            }
            return (dist_return);

        }
        /// <summary>
        /// Returns the distance between 2 Geolocations in either Miles (M) Kilometers (K) nm (N) or (R) Radians
        /// </summary>
        /// <param name="start">Start Coordinate</param>
        /// <param name="end">End Coordinate</param>
        /// <param name="unit">M/N or K</param>
        /// <returns></returns>
        public static double getdistance(GeoLocation start, GeoLocation end, char unit)
        {
            double dist_return;
            double theta = end.Longitude - start.Longitude;
            double dist = Math.Sin(DegreesToRadians(start.Latitude)) * Math.Sin(DegreesToRadians(end.Latitude)) + Math.Cos(DegreesToRadians(start.Latitude)) * Math.Cos(DegreesToRadians(end.Latitude)) * Math.Cos(DegreesToRadians(theta));
            dist = Math.Acos(dist); //radians
            double dist_deg = RadiansToDegrees(dist);
            double dist_miles = dist_deg * 60 * 1.1515;
            switch (unit)
            {
                case 'M':
                case 'm':
                    dist_return = dist_miles;
                    break;
                case 'K':
                case 'k':
                    dist_return = dist_miles * 1.609344;
                    break;
                case 'N':
                case 'n':
                    dist_return = dist_miles * 0.8684;
                    break;
                case 'R':
                case 'r':
                    dist_return = dist;
                    break;

                default:
                    Debug.WriteLine("ERROR: Unknown unit specifier: " + unit, "getdistance");
                    dist_return = 0;
                    break;
            }
            return (dist_return);

        }
        /// <summary>
        /// Returns bearing in degrees from lat1/lon1 to lat2/lon2 coordinates
        /// </summary>
        /// <param name="lat1"></param>
        /// <param name="lon1"></param>
        /// <param name="lat2"></param>
        /// <param name="lon2"></param>
        /// <returns></returns>
        public static double bearing(double lat1, double lon1, double lat2, double lon2)
        {
            lat1 = DegreesToRadians(lat1);
            lat2 = DegreesToRadians(lat2);
            lon1 = DegreesToRadians(lon1);
            lon2 = DegreesToRadians(lon2);

            double x = Math.Cos(lat1) * Math.Sin(lat2) - Math.Sin(lat1) * Math.Cos(lat2) * Math.Cos(lon2 - lon1);
            double y = Math.Sin(lon2 - lon1) * Math.Cos(lat2);
            double abs_brng = Math.Atan2(y, x);
            double brng = d_mod(abs_brng + 2 * Math.PI, 2 * Math.PI);
            brng = RadiansToDegrees(brng);
            return brng;
        }
        /// <summary>
        /// Returns bearing in degrees from start Geolocation to end Geolocation
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <returns></returns>
        public static double bearing(GeoLocation start, GeoLocation end)
        {
            double lat1 = DegreesToRadians(start.Latitude);
            double lat2 = DegreesToRadians(end.Latitude);
            double lon1 = DegreesToRadians(start.Longitude);
            double lon2 = DegreesToRadians(end.Longitude);

            double x = Math.Cos(lat1) * Math.Sin(lat2) - Math.Sin(lat1) * Math.Cos(lat2) * Math.Cos(lon2 - lon1);
            double y = Math.Sin(lon2 - lon1) * Math.Cos(lat2);
            double abs_brng = Math.Atan2(y, x);
            double brng = d_mod(abs_brng + 2 * Math.PI, 2 * Math.PI);
            brng = RadiansToDegrees(brng);
            return brng;
        }
        /// <summary>
        ///  Returns bearing in radians from start lat/lon to end lat/lon coordinates
        /// </summary>
        /// <param name="lat1"></param>
        /// <param name="lon1"></param>
        /// <param name="lat2"></param>
        /// <param name="lon2"></param>
        /// <returns></returns>
        public static double radian_bearing(double lat1, double lon1, double lat2, double lon2)
        {
            lat1 = DegreesToRadians(lat1);
            lat2 = DegreesToRadians(lat2);
            lon1 = DegreesToRadians(lon1);
            lon2 = DegreesToRadians(lon2);

            double x = Math.Cos(lat1) * Math.Sin(lat2) - Math.Sin(lat1) * Math.Cos(lat2) * Math.Cos(lon2 - lon1);
            double y = Math.Sin(lon2 - lon1) * Math.Cos(lat2);
            double abs_brng = Math.Atan2(y, x);
            double brng = d_mod(abs_brng + 2 * Math.PI, 2 * Math.PI);
            return brng;
        }
        /// <summary>
        /// Returns bearing in radians from start Geolocation to end Geolocation
        /// </summary>
        /// <param name="start"></param>
        /// <param name="end"></param>
        /// <returns></returns>
        public static double radian_bearing(GeoLocation start, GeoLocation end)
        {
            double lat1 = DegreesToRadians(start.Latitude);
            double lat2 = DegreesToRadians(end.Latitude);
            double lon1 = DegreesToRadians(start.Longitude);
            double lon2 = DegreesToRadians(end.Longitude);

            double x = Math.Cos(lat1) * Math.Sin(lat2) - Math.Sin(lat1) * Math.Cos(lat2) * Math.Cos(lon2 - lon1);
            double y = Math.Sin(lon2 - lon1) * Math.Cos(lat2);
            double abs_brng = Math.Atan2(y, x);
            double brng = d_mod(abs_brng + 2 * Math.PI, 2 * Math.PI);
            return brng;
        }
        /// <summary>
        /// Calculates the intersection point3 in degrees from point1, crs13, point2, crs23
        /// </summary>
        /// <param name="point1_deg"></param>
        /// <param name="crs13_radians"></param>
        /// <param name="point2_deg"></param>
        /// <param name="crs23_radians"></param>
        /// <returns></returns>
        public static GeoLocation GetIntersection(GeoLocation point1_deg, double crs13_radians, GeoLocation point2_deg, double crs23_radians)
        {
            double lat1 = DegreesToRadians(point1_deg.Latitude);
            double lon1 = DegreesToRadians(point1_deg.Longitude);
            double lat2 = DegreesToRadians(point2_deg.Latitude);
            double lon2 = DegreesToRadians(point2_deg.Longitude);
            double crs12, crs21, ang1, ang2, ang3, dst13, lat3, lon3, dlon;
            GeoLocation point3 = new GeoLocation();

            double delta_lat = lat2 - lat1;
            double delta_lon = lon2 - lon1;


            //double dst12 = 2 * Math.Asin(Math.Pow(Math.Sqrt(Math.Sin((lat1 - lat2) / 2),2) + Math.Cos(lat1) * Math.Cos(lat2) * Math.Pow(Math.Sin((lon1 - lon2) / 2),2)));
            double dst12 = 2 * Math.Asin(Math.Sqrt(Math.Sin(delta_lat / 2) * Math.Sin(delta_lat / 2) + Math.Cos(lat1) * Math.Cos(lat2) * Math.Sin(delta_lon / 2) * Math.Sin(delta_lon / 2)));
            double dist12 = geolib.RadiansToDegrees(dst12) * 60 * 1.1515 * 0.8684;
           // Debug.WriteLine("dst12: " + dist12, "GetIntersection");

            double ang_a = Math.Acos((Math.Sin(lat2) - Math.Sin(lat1) * Math.Cos(dst12)) / (Math.Sin(dst12) * Math.Cos(lat1)));
            double ang_b = Math.Acos((Math.Sin(lat1) - Math.Sin(lat2) * Math.Cos(dst12)) / (Math.Sin(dst12) * Math.Cos(lat2)));

            if (Math.Sin(lon2 - lon1) > 0)
            {
                crs12 = ang_a;
                crs21 = 2 * Math.PI - ang_b;
            }
            else
            {
                crs12 = 2 * Math.PI - ang_a;
                crs21 = ang_b;
            }
            ang1 = d_mod(crs13_radians - crs12 + Math.PI, 2 * Math.PI) - Math.PI;
            ang2 = d_mod(crs21 - crs23_radians + Math.PI, 2 * Math.PI) - Math.PI;
            //Debug.WriteLine("Angle2-1-3: " + geolib.RadiansToDegrees(ang1) + " Angle1-2-3: " + geolib.RadiansToDegrees(ang2), "GetIntersection");

            if (Math.Sin(ang1) == 0 && Math.Sin(ang2) == 0)
            {
                //"infinity of intersections"
                Debug.WriteLine("ERROR: infinity of intersections found.", "GetIntersection");
                point3.Latitude = 99; //return an invalid degree to signal failure - which means lines are parallel
                return point3;
            }
            else if ((Math.Sin(ang1) * Math.Sin(ang2)) < 0)
            {
                //"intersection ambiguous"
                Debug.WriteLine("ERROR: intersection ambiguous found.", "GetIntersection");
                point3.Latitude = 99; //return an invalid degree to signal failure  - which means lines are parallel
                return point3;
            }
            else
            {
                ang1 = Math.Abs(ang1);
                ang2 = Math.Abs(ang2);
                ang3 = Math.Acos(-Math.Cos(ang1) * Math.Cos(ang2) + Math.Sin(ang1) * Math.Sin(ang2) * Math.Cos(dst12));
                dst13 = Math.Atan2(Math.Sin(dst12) * Math.Sin(ang1) * Math.Sin(ang2), Math.Cos(ang2) + Math.Cos(ang1) * Math.Cos(ang3));
                lat3 = Math.Asin(Math.Sin(lat1) * Math.Cos(dst13) + Math.Cos(lat1) * Math.Sin(dst13) * Math.Cos(crs13_radians));
                dlon = Math.Atan2(Math.Sin(crs13_radians) * Math.Sin(dst13) * Math.Cos(lat1), Math.Cos(dst13) - Math.Sin(lat1) * Math.Sin(lat3));
                lon3 = d_mod(lon1 + dlon + Math.PI, 2 * Math.PI) - Math.PI;
                point3.Latitude = RadiansToDegrees(lat3);
                point3.Longitude = RadiansToDegrees(lon3);
                //Debug.WriteLine("Lat: " + point3.Latitude + " Lon: " + point3.Longitude, "GetIntersection");
                return point3;
            }
        }


        public static double DegreesToRadians(double degrees)
         {
             const double degToRadFactor = Math.PI / 180;
             return degrees * degToRadFactor;
         }

         public static double RadiansToDegrees(double radians)
         {
             const double radToDegFactor = 180 / Math.PI;
             return radians * radToDegFactor;
         }
        
        /// <summary>
        /// Converts KM - NM
        /// </summary>
        /// <param name="km"></param>
        /// <returns>double</returns>
        public static double KMToNm(double km)
         {
             double nm = km * 0.539957;
             return nm;
         }

        /// <summary>
        /// Converts Nm - KM
        /// </summary>
        /// <param name="nm"></param>
        /// <returns>double</returns>
        public static double NmToKm(double nm)
        {
            double km = nm * 1.852;
            return km;
        }
        /// <summary>
        /// Converts M to Feet
        /// </summary>
        /// <param name="m"></param>
        /// <returns>Feet double.</returns>
        public static double MtoFeet(double m)
        {
            double ft = m * 3.28084;
            return ft;
        }
        /// <summary>
        /// converts feet to meters
        /// </summary>
        /// <param name="ft"></param>
        /// <returns>double - meters</returns>
        public static double FeettoM(double ft)
        {
            double m = ft * 0.3048;
            return m;
        }
        /// <summary>
        /// converts Meters to Nautical Miles
        /// </summary>
        /// <param name="meters"></param>
        /// <returns></returns>
        public static double MtoNm(double meters)
        {
            double nm = (meters/1000) * 0.539957;
            return nm;
        }
        /// <summary>
        /// Converts distance in nautical miles to radians
        /// </summary>
        /// <param name="nm"></param>
        /// <returns></returns>
        public static double nm2rad(double nm)
        {
            return (Math.PI / (180 * 60)) * nm;
        }
        /// <summary>
        /// Converts distance in radians to nautical miles
        /// </summary>
        /// <param name="rad"></param>
        /// <returns></returns>
        public static double rad2nm(double rad)
        {
            return ((180 * 60) / Math.PI) * rad;
        }
        /// <summary>
        /// Modulus function for double parameters
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double d_mod(double a, double b)
        {
            //equivalent of a%b for doubles
            //Debug.WriteLine("a=" + a + " b=" + b + " Modulus=" + (a - b * Math.Truncate(a / b)),"d_mod");
            return (a - b * Math.Floor(a / b));
        }
    }
    #region Ray Casting
    public struct Rect
    {
        public double min_lat;
        public double max_lat;
        public double min_lon;
        public double max_lon;
    }
    /// <summary>
    /// Creates a new list of boundary points
    /// </summary>
    public class boundary_list
    {
        private List<GeoLocation> data_list;

        public boundary_list()
        {
            data_list = new List<GeoLocation>();
        }
        public void Addto(GeoLocation point)
        {
            data_list.Add(point);
        }
        /// <summary>
        /// Returns a list of boundary points
        /// </summary>
        /// <returns></returns>
        public List<GeoLocation> Getlist()
        {
            return data_list;
        }
        public void Clearall()
        {
            data_list.Clear();
        }
    }

/// <summary>
/// Checks if coordinate passed is within a boundary loaded into Points
/// </summary>
    public class PolyBoundary
    {
        public List<GeoLocation> Points = new List<GeoLocation>();
        /// <summary>
        /// Returns true if point is contained within the boundary saved in Points
        /// </summary>
        /// <param name="point"></param>
        /// <returns></returns>
        public bool Contains(GeoLocation point)
        {
            //  check if it's in the poly.
            int i;
            int j = this.Points.Count - 1;
            bool contained = false;
            for (i = 0; i < this.Points.Count; j = i++)
            {
                if ((((this.Points[i].Latitude <= point.Latitude) && (point.Latitude < this.Points[j].Latitude)) || ((this.Points[j].Latitude <= point.Latitude) && (point.Latitude < this.Points[i].Latitude))) &&
                   (point.Longitude < (this.Points[j].Longitude - this.Points[i].Longitude) * (point.Latitude - this.Points[i].Latitude) / (this.Points[j].Latitude - this.Points[i].Latitude) + this.Points[i].Longitude))
                    contained = !contained;
            }
            return contained;
        }
    }
    #endregion Ray Casting
    #region Random Parameter Generator
    /// <summary>
    /// Generates type double random numbers within specified limits
    /// </summary>
    public class Random_Parms
    {
        private static Random r_gen = new Random();
        /// <summary>
        /// Returns a random double between min and max
        /// </summary>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <returns></returns>
        public double Getparams(double min, double max)
        {
            return r_gen.NextDouble() * (max - min) + min;
        }
        /// <summary>
        /// Returns an random integer between 0 and max-1
        /// </summary>
        /// <param name="max"></param>
        /// <returns></returns>
        public int Getparams(int max)
        {
            return r_gen.Next(max);
        }
    }
    #endregion Random Parameter Generator
}
