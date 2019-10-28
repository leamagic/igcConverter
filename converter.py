import xml.etree.ElementTree as ET
import sys
import re
import math
import datetime


class StartLongest:
  # Start with longest path, i.e. all points used and got back to n points
  def __init__(self,Points):
  	self.totaldistance = []
  	self.points = list(Points)
  	self.calculateLongestPath(self.points)

  def calculateLongestPath(self, Points):
  	distancepoints = []

  	for i,p in enumerate(Points):
  		if i == 0:
  			distance = 0
  			lastpoint = p
  		else :
  			distance = Haversine(lastpoint, p).meters
  		lastpoint = p
  		distancepoints.append({'distance': distance, 'index': i, 'point': p})
  	self.distancepoints = distancepoints
  	self.totaldistance.append(sum(map(lambda dp: dp['distance'], distancepoints)))

  def getMinIndex(self):

  	# Find the point with minimum distance
  	if len(self.points) < 3: 
  		# We have two or less points, remove the first
  		return 0;
  	# If the first leg is the smalles, remove the first point,
  	# if not remove the one between the smallest
  	smallestd = min(map(lambda dp: dp['distance'], self.distancepoints))
  	if self.distancepoints[1]['distance'] == smallestd:
  		# The first distance is the smallest, remove it
  		print("removing the first point")
  		return 0
  	if self.distancepoints[-1]['distance'] == smallestd:
  		# The last distance is the smallest, remove it
  		print("removing the last point")
  		return len(self.distancepoints) -1

  	# Guess that we need to remove the second element, if not proven otherwise
  	indextoremove = 1
  	mindistance = self.distancepoints[1]['distance'] + self.distancepoints[2]['distance']
  	for i,p in enumerate(self.distancepoints):
  		if i == 0: continue
  		if i == 1: continue
  		candidated = self.distancepoints[i-1]['distance'] + self.distancepoints[i]['distance']
  		if candidated < mindistance:
  			# The index we are at has the shortest path so far with its distance plus
  			# the last point's distance. The last point is hence in the middle of the shortest
  			# path so far
  			indextoremove = i
  			mindistance = candidated
  	return indextoremove

  def calculateAll(self):
  	while self.points:
  		indextoremove = self.getMinIndex()
  		# print("Deleting minindex %s. Left %s" % (indextoremove, len(self.distancepoints)))
  		del self.points[indextoremove]
  		dptoremove = self.distancepoints[indextoremove]
  		# Always make sure the first element has distance 0
  		self.distancepoints[0]['distance'] = 0
  		
  		del self.distancepoints[indextoremove]
  		# If it was the first or last 
  		if indextoremove == 0: 
  			self.totaldistance.append(sum(map(lambda dp: dp['distance'], self.distancepoints)))
  			continue
  		if indextoremove >= len(self.distancepoints): 
  			self.totaldistance.append(sum(map(lambda dp: dp['distance'], self.distancepoints)))
  			continue
  		# We need to recalculate the distance between the two points
  		# surrounding the removed point
  		newdistance = Haversine(self.distancepoints[indextoremove-1]['point'], self.distancepoints[indextoremove]['point']).meters
  		self.distancepoints[indextoremove]['distance'] = newdistance
  		# Then calculate new totaldistance
  		self.totaldistance.append(sum(map(lambda dp: dp['distance'], self.distancepoints)))

  

class Diameter:
  # convex hull (Graham scan by x-coordinate) and diameter of a set of points
  # David Eppstein, UC Irvine, 7 Mar 2002
  def diameter(self, Points):
       '''Given a list of 2d points, returns the pair that's farthest apart.'''
       usedPoints = list(Points)
       usedPoints = map(lambda p: (6371000*math.cos(math.radians(p[1]))*math.cos(math.radians(p[0])),6371000*math.cos(math.radians(p[1]))*math.sin(math.radians(p[0]))), Points)
       diam,pair = max([((p[0]-q[0])**2 + (p[1]-q[1])**2, (p,q))
                        for p,q in self.rotatingCalipers(usedPoints)])
       print("The diameter is %s" % diam)
       found1 = 0
       for i,p in enumerate(Points):
       	if p[0] == pair[1][0] and p[1] == pair[1][1]:
       		print("Found a match in index %d" % i)
       		found1 = i
       found2 = 0
       for i,p in enumerate(Points):
       	if p[0] == pair[0][0] and p[1] == pair[0][1]:
       		print("Found a match in index %d" % i)
       		found2 = i
       #print("Max distance between (%s,%s) and (%s,%s)" % (pair[0][0],pair[0][1],pair[1][0],pair[1][1]))
       print("Hard distance is %d" % Haversine(Points[10], Points[6537]).meters)
       return Haversine(pair[0], pair[1]).meters
       # return diam

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

def calculate2distance4(blines):
	Points = map(lambda bl: (bl['longdec'], bl['latdec']), blines)
	sl = StartLongest(Points)
	sl.calculateAll()
	return sl

def calculate2distance2(blines):
	# Find the points that have max distance
	maxd = 0
	lines = list(blines)
	line = lines.pop()
	start = (line['longdec'], line['latdec'])
	distances = map(lambda bl: Haversine(start, (bl['longdec'], bl['latdec'])).meters, lines)
	indexmax1 = 0
	for index, d in enumerate(distances):
		if d > maxd:
			maxd = d
			indexmax1 = index
	print("Found max in round 1 %d as index %d" % (maxd, indexmax1))
	start = (lines[indexmax1]['longdec'], lines[indexmax1]['latdec'])
	distances2 = map(lambda bl: Haversine(start, (bl['longdec'], bl['latdec'])).meters, lines)
	maxround2  = 0
	indexmax2 = 0
	for index, d in enumerate(distances2):
		if d > maxround2:
			maxround2 = d
			indexmax2 = index
	print("Found max in round 2 %d as index %d" % (maxround2, indexmax2))
	maxround2 = reduce(lambda a, b: max(a, b), distances2)
	print("Found max in round 2 %d" % maxround2)
	return max(maxround2, maxd)


def calculate2distance(blines):
	# Find the points that have max distance
	maxd = 0
	lines = list(blines)
	line = lines.pop()
	while line:
		if not lines: break	
		start = (line['longdec'], line['latdec'])
		#distances = map(lambda bl: Haversine(start, (bl['longdec'], bl['latdec'])).meters, lines)
		distances = []
		for bl in lines:
			maxnow = Pythagoras(start, (bl['longdec'], bl['latdec'])).meters
			#maxnow = Haversine(start, (bl['longdec'], bl['latdec'])).meters
			#maxnow = math.sqrt(math.pow(line['longdec']-bl['longdec'],2) + math.pow(line['latdec']-bl['latdec'],2))
			#maxnow = maxd
			# distances.append(maxnow)
			maxd = max(maxd, maxnow)
		#maxnow = reduce(lambda a, b: max(a, b), distances)
		#maxd = max(maxd, maxnow)
		line = lines.pop()
		# 2min32s with map reduce
		# 2min28swith for-loop then reduce
		# with for-loop and max inlined
		# 2min21s with Pythagoras instead of Haversine
		# 21364m is the distance(?)(open distance)
		# 21170m With convex hull algorithm

	return maxd


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
	#print('Distance between 2 points: %dm' % calculate2distance(blines))
	#print('Distance between 2 points: %dm' % calculate2distance2(blines))
	print('Distance between 2 points: %sm' % calculate2distance3(blines))
	#startlongest = calculate2distance4(blines)
	#startlongest.calculateAll()

	#for i,d in enumerate(startlongest.totaldistance):
	#	print("distance %i: %d" % (i,d))
	
	#print('Total distance: %dm' % startlongest.totaldistance[0])


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