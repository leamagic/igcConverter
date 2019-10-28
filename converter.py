import xml.etree.ElementTree as ET
import sys
import re
import math
import datetime  


RADIUS = 6378137.0 # in meters on the equator

def y2lat(aY):
   return math.degrees(math.atan(math.exp(aY / RADIUS)) * 2 - math.pi/2)

def x2lon(aX):
   return math.degrees(aX / RADIUS)

def lat2y(a):
  return math.log(math.tan(math.pi / 4 + math.radians(a) / 2)) * RADIUS

def lon2x(a):
  return math.radians(a) * RADIUS

class Diameter:
  # convex hull (Graham scan by x-coordinate) and diameter of a set of points
  # David Eppstein, UC Irvine, 7 Mar 2002
  def diameter(self, Points):
    '''Given a list of 2d points, returns the pair that's farthest apart.'''
    usedPoints = map(lambda p: (lon2x(p[0]), lat2y(p[1])), Points)
    diam,pair = max([((p[0]-q[0])**2 + (p[1]-q[1])**2, (p,q)) for p,q in self.rotatingCalipers(usedPoints)])
    lon1 = x2lon(pair[0][0])
    lat1 = y2lat(pair[0][1])
    lon2 = x2lon(pair[1][0])
    lat2 = y2lat(pair[1][1])
    
    print("Max distance between (%s,%s) and (%s,%s)" % (lon1, lat1, lon2, lat2))
    return Haversine((lon1, lat1), (lon2, lat2)).meters

  def __init__(self,Points):
  	self.meters = self.diameter(Points)

  def orientation(self, p,q,r):
    '''Return positive if p-q-r are clockwise, neg if ccw, zero if colinear.'''
    return (q[1]-p[1])*(r[0]-p[0]) - (q[0]-p[0])*(r[1]-p[1])

  def hulls(self, Points):
    '''Graham scan to find upper and lower convex hulls of a set of 2d points.'''
    U = []
    L = []
    Points.sort()
    for p in Points:
        while len(U) > 1 and self.orientation(U[-2],U[-1],p) <= 0: U.pop()
        while len(L) > 1 and self.orientation(L[-2],L[-1],p) >= 0: L.pop()
        U.append(p)
        L.append(p)
    return U,L

  def rotatingCalipers(self, Points):
    '''Given a list of 2d points, finds all ways of sandwiching the points
    between two parallel lines that touch one point each, and yields the sequence
    of pairs of points touched by each pair of lines.'''
    U,L = self.hulls(Points)
    i = 0
    j = len(L) - 1
    while i < len(U) - 1 or j > 0:
        yield U[i],L[j]
        
        # if all the way through one side of hull, advance the other side
        if i == len(U) - 1: j -= 1
        elif j == 0: i += 1
        
        # still points left on both lists, compare slopes of next hull edges
        # being careful to avoid divide-by-zero in slope calculation
        elif (U[i+1][1]-U[i][1])*(L[j][0]-L[j-1][0]) > \
                (L[j][1]-L[j-1][1])*(U[i+1][0]-U[i][0]):
            i += 1
        else: j -= 1


class Haversine:
    '''
    use the haversine class to calculate the distance between
    two lon/lat coordnate pairs.

    output distance available in kilometers, meters, miles, and feet.

    example usage: Haversine([lon1,lat1],[lon2,lat2]).feet
    
    '''
    def __init__(self,coord1,coord2):
        lon1,lat1=coord1
        lon2,lat2=coord2
        
        R=6371000                               # radius of Earth in meters
        phi_1=math.radians(lat1)
        phi_2=math.radians(lat2)

        delta_phi=math.radians(lat2-lat1)
        delta_lambda=math.radians(lon2-lon1)

        a=math.sin(delta_phi/2.0)**2+\
           math.cos(phi_1)*math.cos(phi_2)*\
           math.sin(delta_lambda/2.0)**2
        c=2*math.atan2(math.sqrt(a),math.sqrt(1-a))
        
        self.meters=R*c                         # output distance in meters
        self.km=self.meters/1000.0              # output distance in kilometers
        self.miles=self.meters*0.000621371      # output distance in miles
        self.feet=self.miles*5280               # output distance in feet

class Pythagoras:
    '''
    use the pythagoras class to calculate the distance between
    two lon/lat coordnate pairs.

    output distance available in kilometers, meters, miles, and feet.

    example usage: Pythagoras([lon1,lat1],[lon2,lat2]).feet
    
    '''
    def __init__(self,coord1,coord2):
        lon1,lat1=coord1
        lon2,lat2=coord2
        
        self.meters = math.sqrt(math.pow(lon2-lon1,2) + math.pow(lat2-lat1,2))

def findMaxHeight(blines):
	maxheight = 0
	return reduce(lambda maxheight, h: max(maxheight, h), map(lambda l: l['presalt'], blines))

def calculateDuration(blines):
	if not blines: return 0
	return blines[-1]['time'] - blines[0]['time'] 

def calculate2distance3(blines):
	Points = map(lambda bl: (bl['longdec'], bl['latdec']), blines)
	return Diameter(Points).meters

filepath = sys.argv[1]

with open(filepath, 'r') as file:
	line = file.readline()
	numofb = 0
	blines = []
	while line:
		char = line[0]
		if line.startswith('HFDTE'):
			# This is the date of the first sample
			found = re.search('HFDTE(\d{2})(\d{2})(\d{2})', line)
			if not found: break
			year = int(found.group(3)) + 2000
			month = int(found.group(2))
			day = int(found.group(1)) 
			startdate = datetime.date(year, month, day)
		if char == 'B':
			numofb += 1
			# 6 digits for the time, 7 latitude, 8 longitude A 5 Press Alt 5 GNSS alt
			found = re.search('B(\d{6})(\d{2})(\d{5})[NS](\d{3})(\d{5})[EW]A(\d{5})(\d{5})', line)
			if not found: break
			time = int(found.group(1))
			latdeg = found.group(2)
			latmin = found.group(3)
			longdeg = found.group(4)
			longmin = found.group(5)
			latdec = int(latdeg) + float(latmin)/60000
			longdec = int(longdeg) + float(longmin)/60000
			presalt = int(found.group(6))
			gnssalt = int(found.group(7))
			# print("Time is %d, latitude: %f, longitude: %s, height: %d" % (time, latdec, longdec, gnssalt))
			blines.append({'time': time, 'latdec': latdec, 'longdec': longdec, 'presalt': presalt, 'gnssalt': gnssalt})
		line = file.readline()
	print('Max heigth: %d' % findMaxHeight(blines))
	print('Duration: %d seconds' % calculateDuration(blines))
	print('Flight date: %s' % startdate)
	print('Distance between 2 points: %sm' % calculate2distance3(blines))
	

print ('Number of Bs: %d' % numofb)



exit(0)








# create the file structure
data = ET.Element('data')  
items = ET.SubElement(data, 'items')  
item1 = ET.SubElement(items, 'item')  
item2 = ET.SubElement(items, 'item')  
item1.set('name','item1')  
item2.set('name','item2')  
item1.text = 'item1abc'  
item2.text = 'item2abc'

# create a new XML file with the results
mydata = ET.tostring(data)  
myfile = open("items2.xml", "w")  
myfile.write(mydata)  