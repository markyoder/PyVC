#!/usr/bin/env python
import math

class LatLonPoint(object):
    def __init__(self, lat = None, lon = None):
        self.lat = lat
        self.lon = lon
    
    def __str__(self):
        return 'lat: %f, lon: %f'%(self.lat, self.lon)

class MyVector(object):
    def __init__(self, x,y,z):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        #self.x = round(x,5)
        #self.y = round(y,5)
        #self.z = round(z,5)
    
    def norm(self):
        norFac = 1.0/self.length()
        
        return MyVector(norFac*self.x, norFac*self.y, norFac*self.z)
    
    def printMe(self, prepend=''):
        print '%s%f %f %f'%(prepend,self.x, self.y, self.z)
        
    def times(self, num):
        #print 'times', num, MyVector(float(num)*self.x, float(num)*self.y, float(num)*self.z)
        return MyVector(float(num)*self.x, float(num)*self.y, float(num)*self.z)
        
    def cross(self, vec):
        x = self.y*vec.z - vec.y*self.z
        y = -1* (self.x*vec.z - vec.x*self.z)
        z = self.x*vec.y - vec.x*self.y
        return MyVector(x,y,z)
    
    def dot(self, vec):
        return self.x * vec.x + self.y * vec.y + self.z * vec.z
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __add__(self, other):
        #print 'add'
        if type(other) == type(self):
            return MyVector(self.x+other.x, self.y+other.y, self.z+other.z)
        else:
            raise TypeError('Only vectors can be added to vectors')
    
    def __rsub__(self, other):
        if type(other) == type(self):
            return MyVector(other.x-self.x, other.y-self.y, other.z-self.z)
        else:
            raise TypeError('Only vectors can be subtracted from vectors')
    
    def __sub__(self, other):
        if type(other) == type(self):
            return MyVector(self.x-other.x, self.y-other.y, self.z-other.z)
        else:
            raise TypeError('Only vectors can be subtracted from vectors')
    
    def __rmul__(self, other):
        return self.__mul__(other)
        
    def __mul__(self, other):
        #print 'mul', other
        if type(other) == type(self):
            return self.cross(other)
        elif type(other) is int or type(other) is float:
            #print self.times(other)
            return self.times(other)
        else:
            raise TypeError('Only vectors or numbers can be multiplied with vectors')
    
    def __rmod__(self, other):
        return self.__mod__(other)
    
    def __mod__(self, other):
        if type(other) == type(self):
            return self.dot(other)
        else:
            raise TypeError('Only vectors can be dotted into vectors')
            
    def __lt__(self, other):
        if type(other) == type(self):
            return math.acos( self.dot(other) / (self.length() * other.length()))
        else:
            raise TypeError('Angles can be found ony for two vectors')
            
        
    def length(self):
        #print math.sqrt( self.x**2.0 + self.y**2.0 + self.z**2.0 )
        return math.sqrt( self.x**2.0 + self.y**2.0 + self.z**2.0 )
    
    def __str__(self):
        return '<%f, %f, %f>'%(self.x , self.y , self.z)
        
            