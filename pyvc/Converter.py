#!/usr/bin/env python

import math
import calendar
from geographiclib.geodesic import Geodesic

class Converter:
        def __init__(self):
                self.lat0 = 31.5
                self.lon0 = -126.0
                self.earth_radius = self.km_m(6371.0)
        
        def distanceOnEarthsSurface(self, lat1, lon1, lat2, lon2):
            conv = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
            return conv["s12"]
        
        def latlon_xy_2pt(self, lat1, lon1, lat2, lon2):
            conv = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
            return (conv["s12"]*math.sin(self.deg_rad(conv["azi1"])), conv["s12"]*math.cos(self.deg_rad(conv["azi1"])))
        
        def yearDecimalToYearMonthDay(self,yearDecimal,time=False):
                decimal, year = math.modf(yearDecimal)
                
                decimal, month = math.modf(12.0*decimal)
            
                #print float(year)%10000.0, year
            
                decimal, day = math.modf(calendar.monthrange(1972, int(month + 1))[1]*decimal)
                
                '''
                if year <= 9999:
                        decimal, day = math.modf(calendar.monthrange(int(year), int(month + 1))[1]*decimal)
                else:
                        remainder, thousands = math.modf(yearDecimal/10000)
                        
                        year_mod = year - thousands*10000
                    
                        #print year_mod + 1.0, year
                        decimal, day = math.modf(calendar.monthrange(int(year_mod + 1.0), int(month + 1))[1]*decimal)
                '''
                
                decimal, hour = math.modf(24.0*decimal)
                
                decimal, minute = math.modf(60.0*decimal)
                
                decimal, second = math.modf(60.0*decimal)
                
                #print year, month + 1, day + 1, hour + 1 
                
                if time:
                    return int(year), int(month + 1), int(day + 1), int(hour), int(minute), int(second), decimal
                else:
                    return int(year), int(month + 1), int(day + 1)
        
        def year_sec(self, year):
            return 365.0 * 24.0 * 60.0 * 60.0 * year
        
        def sec_year(self, sec):
            return (1.0/365.0) * (1.0/24.0) * (1.0/60.0) * (1.0/60.0) * sec
            
        def km_m(self, km):
                return km * 10.0**(3.0)
        
        def m_km(self, m):
                return m * 10.0**(-3.0)
        
        def deg_rad(self, deg):
                return deg * (math.pi/180.0)
        
        def rad_deg(self, rad):
                return rad * (180.0/math.pi)
    
        def setLatLon0(self, lat0, lon0):
                self.lat0 = lat0
                self.lon0 = lon0
        
        def latlon_xy(self, lat, lon):
            conv = Geodesic.WGS84.Inverse(self.lat0, self.lon0, lat, lon)
            '''
            print conv["lat1"], conv["lon1"]
            print conv["lat2"], conv["lon2"]
            print conv["s12"]
            print conv["azi1"], conv["azi2"]
            
            print Geodesic.WGS84.Direct(self.lat0, self.lon0, conv["azi1"], conv["s12"])
            
            print conv["s12"]*math.sin(self.deg_rad(conv["azi1"])), conv["s12"]*math.cos(self.deg_rad(conv["azi1"]))
            '''
            #print conv
            return (conv["s12"]*math.sin(self.deg_rad(conv["azi1"])), conv["s12"]*math.cos(self.deg_rad(conv["azi1"])))
            
        def latlon_xy_old(self, lat, lon):
                x = self.arclen((self.lat0, self.lon0), (self.lat0, lon))
                y = self.arclen((self.lat0,0),(lat,0))
                
                if (lon < self.lon0):
                    x *= -1
                
                if (lat < self.lat0):
                    y *= -1
                
                return (x,y)
                
        def xy_latlon(self, x, y):
            s12 = (x**2.0 + y**2.0)**(0.5)
            azi1 = self.rad_deg(math.atan2(x,y))
            conv = Geodesic.WGS84.Direct(self.lat0, self.lon0, azi1, s12)
            
            #print conv
            return (conv["lat2"], conv["lon2"])
            
        def xy_latlon_old(self, x, y):
            #print x, y, self.lat0, self.lon0, self.earth_radius
            
                new_lat = self.rad_deg(y/self.earth_radius)+self.lat0;

                new_lon = 2.0 * self.rad_deg(math.asin(math.sin(x/(2.0 * self.earth_radius))/math.cos(self.deg_rad(new_lat))))+self.lon0;

                return new_lat,new_lon
        
        def arclen(self, pt1, pt2):
                # pts are always (lat,lon)

                dlon = self.deg_rad(pt2[1]-pt1[1])
                dlat = self.deg_rad(pt2[0]-pt1[0])
                lat1 = self.deg_rad(pt1[0])
                lat2 = self.deg_rad(pt2[0])

                a = math.sin(dlat/2.0)**2.0 + math.cos(lat1)*math.cos(lat2)*( math.sin(dlon/2.0)**2.0 )
                c = 2.0*math.atan2(math.sqrt(a), math.sqrt(1.0-a))
                d = self.earth_radius*c

                return d