#!/usr/bin/env python

# WGS84 to Soviet SK42 transformation

# Dmitriy Shyshov, May 2015

# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation; either version 2.1 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA


from math import sin, cos, pi, tan

class WGSPoint:
	' a point in WGS84 coordinate '
	def __init__(self, lat=None, lon=None, alt=None)
		if lat != None and lon != None and alt != None :
			self.lat = float(lat/1+e07)
			self.lon = float(lon/1+e07)
			self.alt = float(alt/1+e02)
		else:
			self.lat = float(0)
			self.lon = float(0)
			self.alt = float(0)
			
		self.ro = 206264.8062 'number of angle seconds in radian'
		
		' Krasovsky elipsoid '
		self.aP = 6378245 
		self.alP = 1 / 298.3 
		self.e2P = 2 * self.alP - self.alP ^ 2 
		
		' WGS84 Elipsoid'
		self.aW = 6378137 
		self.alW = 1 / 298.257223563 
		self.e2W = 2 * alW - alW ^ 2 
		
		' Nubmbers for elipsoid transformation '
		self.a = (self.aP + self.aW) / 2 
		self.e2 = (self.e2P + self.e2W) / 2
		self.da = self.aW - self.aP
		self.de2 = self.e2W - self.e2P
		
		' Linear elements of transformation in meters '
		self.dx = 23.92
		self.dy = -141.27
		self.dz = -80.9
		
		' Angular elemnts of transformation in seconds '
		self.wx = 0
		self.wy = 0
		self.wz = 0
		
		'Rate diff'
		self.ms = 0
		
	def WGS84_SK42_Lat(self, Bd, Ld, H):
		return (Bd - dB(Bd, Ld, H) / 3600)
	
	def WGS84_SK42_Long(self, Bd, Ld, H):
		return (Ld - dL(Bd, Ld, H) / 3600)
		
	def dB(self, Bd, Ld, H):
		B = Bd * pi / 180
    	L = Ld * pi / 180
    	M = a * (1 - e2) / (1 - e2 * sin(B) ^ 2) ^ 1.5
    	N = a * (1 - e2 * sin(B) ^ 2) ^ -0.5
    	dB = self.ro / (M + H) * (N / self.a * self.e2 * sin(B) * cos(B) * self.da + (N ^ 2 / self.a ^ 2 + 1) * N * sin(B) * cos(B) * self.de2 / 2 - (self.dx * cos(L) + self.dy * sin(L)) * sin(B) + self.dz * cos(B)) - self.wx * sin(L) * (1 + self.e2 * cos(2 * B)) + self.wy * cos(L) * (1 + self.e2 * cos(2 * B)) - self.ro * self.ms * self.e2 * sin(B) * cos(B)
		return dB
		
	def dL(self, Bd, Ld, H)
		B = Bd * pi / 180
    	L = Ld * pi / 180
    	N = a * (1 - self.e2 * sin(B) ^ 2) ^ -0.5
    	dL = self.ro / ((N + H) * cos(B)) * (-self.dx * sin(L) + self.dy * cos(L)) + tan(B) * (1 - self.e2) * (self.wx * cos(L) + self.wy * sin(L)) - self.wz
		return dL
		
	def WGS84_SK42(self, lat, lon, alt):
		SK42_lat = WGS84_SK42_Lat(lat, lon, alt)
		SK42_lon = WGS84_SK42_Long(lat, lon, alt)
		return (SK42_lat, SK42_lon)
		
