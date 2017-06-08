def line(s, left, top):

    s.Line(point1=(-2 * left, -top), point2=(2 * left, -top))

def rectangle(s, left, top, thickness):
    s.rectangle(point1=(left,top), point2=(-left,top+thickness))


def sketch(s, array):

# Given 2D array of xy coordinates, sketches polygon using shapes
# ....a = mdb.models['standard'].rootAssembly
# ....s = mdb.models['standard'].ConstrainedSketch(name='__profile__', sheetSize=4.0)

    for i in range(0, len(array[0]) - 1):
        s.Line(point1=(array[0][i], array[1][i]), point2=(array[0][i
               + 1], array[1][i + 1]))
def advancedarrayspolygon(inradius, thickness, sides,distance):
    innervertices = vertices(0,0, inradius, sides)
    outervertices = vertices(0,0, inradius+thickness, sides)
    
    innervleft = [innervertices[0][0:len(innervertices[0])//2],innervertices[1][0:len(innervertices[0])//2]]
    innervright = [innervertices[0][len(innervertices[0])//2:len(innervertices[0])],innervertices[1][len(innervertices[0])//2:len(innervertices[0])]]
    outervleft = [outervertices[0][0:len(outervertices[0])//2],outervertices[1][0:len(outervertices[0])//2]]
    outervright = [outervertices[0][len(outervertices[0])//2:len(outervertices[0])],outervertices[1][len(outervertices[0])//2:len(outervertices[0])]]
    
    topnodeouter = max(outervertices[1])
    topnodeinner = max(innervertices[1])
    leftnodeouter = max(outervertices[1])
    
    translatedleft1 = list(map(lambda x: x-leftnodeouter-distance, innervleft[0]))
    translatedright1 = list(map(lambda x: x+leftnodeouter+distance, innervright[0]))
    innerv = [translatedleft1, innervleft[1]]
    innerv[0].extend(translatedright1)
    innerv[1].extend(innervright[1])
    
	
    translatedleft2 = list(map(lambda x: x-leftnodeouter-distance, outervleft[0]))
    translatedright2 = list(map(lambda x: x+leftnodeouter+distance, outervright[0]))
    outerv = [translatedleft2, outervleft[1]]
    outerv[0].extend(translatedright2)
    outerv[1].extend(outervright[1])
    return [innerv,outerv]
def advancedpolygon(
    inradius,
    thickness,
    sides,
    distance,
    ):
	print("advanced polygon started")
	innervertices = vertices(0,0, inradius, sides)
	outervertices = vertices(0,0, inradius+thickness, sides)

	innervleft = [innervertices[0][0:len(innervertices[0])//2],innervertices[1][0:len(innervertices[0])//2]]
	innervright = [innervertices[0][len(innervertices[0])//2:len(innervertices[0])-1],innervertices[1][len(innervertices[0])//2:len(innervertices[0])-1]]

	outervleft = [outervertices[0][0:len(outervertices[0])//2],outervertices[1][0:len(outervertices[0])//2]]
	outervright = [outervertices[0][len(outervertices[0])//2:len(outervertices[0])-1],outervertices[1][len(outervertices[0])//2:len(outervertices[0])-1]]

	topnodeouter = max(outervertices[1])
	topnodeinner = max(innervertices[1])
	leftnodeouter = max(outervertices[1])
	
	arrayleft = [[0,0],[topnodeinner,topnodeouter]]
	arrayleft[0].extend(outervleft[0])
	arrayleft[1].extend(outervleft[1])
	arrayleft[0].extend([0,0])
	arrayleft[1].extend([-topnodeouter,-topnodeinner])
	arrayleft[0].extend(list(reversed(innervleft[0])))
	arrayleft[1].extend(list(reversed(innervleft[1])))
	arrayleft[0].append(0)
	arrayleft[1].append(topnodeinner)
	print(arrayleft)
	print("/n")
	
	translatedleft = list(map(lambda x: x-leftnodeouter-distance, arrayleft[0]))
	
	arrayleftf = [translatedleft, arrayleft[1]]
	
	arrayright = [[0,0],[topnodeouter,topnodeinner]]
	arrayright[0].extend(list(reversed(innervright[0])))
	arrayright[1].extend(list(reversed(innervright[1])))
	arrayright[0].extend([0,0])
	arrayright[1].extend([-topnodeinner,-topnodeouter])
	arrayright[0].extend(outervright[0])
	arrayright[1].extend(outervright[1])
	arrayright[0].append(0)
	arrayright[1].append(topnodeouter)

	translatedright = list(map(lambda x: x+leftnodeouter+distance, arrayright[0]))
        arrayrightf = [translatedright, arrayright[1]]
    
	return [arrayleftf, arrayrightf]




def vertices(
    startx,
    starty,
    radius,
    sides,
    ):

# Returns the coordinates of the vertices of a (sides)-dimensional polygon
# centered at startx, starty of a certain radius (currently the code only works for 0,0)

    mainangle = 2 * 3.14159265 / sides
    xarray = []
    yarray = []
    for i in range(0, sides):
        angle = mainangle * (i + 0.5)

        sincos = trig(angle)
        xval = startx + radius * sincos[0]
        yval = starty + radius * sincos[1]

        xarray.append(xval)
        yarray.append(yval)

   # append starting point so line can return to starting vertice

    xarray.append(xarray[0])
    yarray.append(yarray[0])
    return [xarray, yarray]


def trig(x):
# Returns sin and cos of angle x (rad)
# http://lab.polygonal.de/2007/07/18/fast-and-accurate-sinecosine-approximation/
# Low degree of precision required because CAE should automatically get correct angle once sides
# are constrained

    if x < -3.14159265:
        x += 6.28318531
    elif x > 3.14159265:
        x -= 6.28318531

# compute sine

    if x < 0:
        sin = 1.27323954 * x + .405284735 * x * x
        if sin < 0:
            sin = .225 * (sin * -sin - sin) + sin
        else:
            sin = .225 * (sin * sin - sin) + sin
    else:
        sin = 1.27323954 * x - .405284735 * x * x
        if sin < 0:
            sin = .225 * (sin * -sin - sin) + sin
        else:
            sin = .225 * (sin * sin - sin) + sin

# compute cosine: sin(x + PI/2) = cos(x)

    x += 1.57079632
    if x > 3.14159265:
        x -= 6.28318531
    if x < 0:
        cos = 1.27323954 * x + .405284735 * x * x
        if cos < 0:
            cos = .225 * (cos * -cos - cos) + cos
        else:
            cos = .225 * (cos * cos - cos) + cos
    else:
        cos = 1.27323954 * x - .405284735 * x * x
        if cos < 0:
            cos = .225 * (cos * -cos - cos) + cos
        else:
            cos = .225 * (cos * cos - cos) + cos

    return (sin, cos)  # return tuple


def radiusgen(perimeter, sides):
# Assume for fair comparison that perimeter is equal
# Calculates the radius of the polygon given the perimeter and number of sides

    sincos = trig(3.14159265 / sides)
    radius = perimeter / sides / (2 * sincos[0])
    return radius

    xarray.append(xarray[0])
    yarray.append(yarray[0])
    return [xarray, yarray]

def sidelength(radius, sides):
#Input: Radius of polygon and number of sides
#Calculates the side length of a "sides"-sided polygon of radius 
    	import math 
    	side = 2*radius*math.sin(180/sides)
	return side

def bc_bot(mdb, vertices, thickness):
#Input: Model database, array of vertices of L/R polygon, thickness of polygon
#Applies boundary conditions and rigid body constraints
	from abaqus import *
	from abaqusConstants import *
	import regionToolset

	output = []
	num_vert = len(vertices[0])
	a = mdb.models['standard'].rootAssembly

	bottom = min(vertices[0][1])
	left = min(vertices[0][0])

	#Rigid body definition for top plane
        #a = mdb.models['standard'].rootAssembly
        #e1 = a.instances['planetop-1'].edges
	#coord = e1.getClosest(coordinates=((0,-bottom+thickness,0),))

    	#RP = a.ReferencePoint(point=a.instances['planetop-1'].InterestingPoint(edge=e1.findAt(
        #	coordinates=coord[0][1]), rule=MIDDLE))
	#RP_id = RP.id
        #region2=a.instances['planetop-1'].sets['planar2']
     	#r1 = a.referencePoints
    	#refPoints1 = mdb.models['standard'].rootAssembly.referencePoints[RP_id]		
    	#region1=regionToolset.Region(referencePoints=(refPoints1,))
    	#mdb.models['standard'].RigidBody(name='Constraint-1', refPointRegion=region1, 
        #	bodyRegion=region2)
	

	#Rigid body definition for bot plane
        a = mdb.models['standard'].rootAssembly
        e2 = a.instances['plane-1'].edges
	coord = e2.getClosest(coordinates=((0,bottom+thickness,0),))

    	RP=a.ReferencePoint(point=a.instances['plane-1'].InterestingPoint(edge=e2.findAt(
        	coordinates=coord[0][1]), rule=MIDDLE))
	RP_id = RP.id
        region2=a.instances['plane-1'].sets['planar']
     	r1 = a.referencePoints
    	refPoints1 = mdb.models['standard'].rootAssembly.referencePoints[RP_id]
    	region1=regionToolset.Region(referencePoints=(refPoints1,))
    	mdb.models['standard'].RigidBody(name='Constraint-2', refPointRegion=region1, 
        	bodyRegion=region2)

	mdb.models['standard'].DisplacementBC(name='planelock', 
        	createStepName='Apply load', region=region1, u1=SET, u2=SET, ur3=SET, 
        	amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        	localCsys=None)
  
	s2 = a.instances['plane-1'].edges
	coord = s2.getClosest(coordinates=((-left,bottom-thickness/2,0),))
	side1Edges1 = s2.findAt((coord[0][1],))
	coord = s2.getClosest(coordinates=((left,bottom-thickness/2,0),))
	side1Edges1 += s2.findAt((coord[0][1],))
	
	s2 = a.instances['planetop-1'].edges
	coord = s2.getClosest(coordinates=((-left,-bottom+thickness/2,0),))
	side1Edges1 += s2.findAt((coord[0][1],))
	coord = s2.getClosest(coordinates=((left,-bottom+thickness/2,0),))
	side1Edges1 += s2.findAt((coord[0][1],))

	s2 = a.instances['Frame-1'].edges
	coord = s2.getClosest(coordinates=((-left,bottom+thickness/2,0),))
	side1Edges1 += s2.findAt((coord[0][1],))
	coord = s2.getClosest(coordinates=((left,bottom+thickness/2,0),))
	side1Edges1 += s2.findAt((coord[0][1],))
	coord = s2.getClosest(coordinates=((-left,-bottom-thickness/2,0),))
	side1Edges1 += s2.findAt((coord[0][1],))
	coord = s2.getClosest(coordinates=((left,-bottom-thickness/2,0),))
	side1Edges1 += s2.findAt((coord[0][1],))

	regionbot = a.Set(edges=side1Edges1, name='xsym')

	mdb.models['standard'].XsymmBC(name='xsym', createStepName='Initial', 
        region=regionbot, localCsys=None)

def loader(mdb, radius, thickness, pressure):
#Input: Model database, radius of the polygon, thickness of the top plane, and pressure 
#Applies pressure to top edge of planetop
	from abaqus import *
	from abaqusConstants import *	
	a = mdb.models['standard'].rootAssembly
	
	s2 = a.instances['planetop-1'].edges
	coord = s2.getClosest(coordinates=((0,radius+2*thickness,0),))
	side1Edges1 = s2.findAt((coord[0][1],))
	region = a.Surface(side1Edges=side1Edges1, name='planetop-top')


	data = tuple((0.001 * t / 100, 3.4 * (1 - t / 100.0)) for t in range(100))
   	mdb.models['standard'].TabularAmplitude(name='Blast Amplitude', timeSpan=TOTAL, smooth=SOLVER_DEFAULT, 
        	data=data+((0.001, 0.0), (0.002, 0.0), (100000, 0.0)))
        mdb.models['standard'].Pressure(name='Load 1', createStepName='Apply load', region=region, magnitude=3.4, amplitude='Blast Amplitude')

def radiuschecker(startx, maxrad, x, y):
#Assumed that starty = zero 
#Returns true if x,y is more than maxrad away from (startx,0) 
	return maxrad<((x-startx)**2+y**2)**0.5

def edgecentre(vex):
#Input: 2d array of vertices for the L/R polygon parts
#Returns the list of coordinates for the points between two vertices
#Note: To be used in conjunction with the edgeselector function
	array1x=[]
	array1y=[]
	array2x=[]
	array2y=[]
	length = len(vex[0][0])
	for i in range(1, length):
		array1x.append((vex[0][0][i]+vex[0][0][i-1])/2)
		array1y.append((vex[0][1][i]+vex[0][1][i-1])/2)
		array2x.append((vex[1][0][i]+vex[1][0][i-1])/2)
		array2y.append((vex[1][1][i]+vex[1][1][i-1])/2)
	array1x.append((vex[0][0][length-1]+vex[0][0][0])/2)
	array1y.append((vex[0][1][length-1]+vex[0][1][0])/2)
	array2x.append((vex[1][0][length-1]+vex[1][0][0])/2)
	array2y.append((vex[1][1][length-1]+vex[1][1][0])/2)


	
	return [[array1x,array1y],[array2x,array2y]]

def edgeselector(mdb, vex, radius, thickness):
#Input: Model database, array of vertices for L/R polygon parts, and radius of the polygon 
#Generates the interacting sets for the contacts
#Then creates the contact interaction between them


	vex = edgecentre(vex)
	from abaqus import *
	from abaqusConstants import *
	a = mdb.models['standard'].rootAssembly
	s1 = a.instances['Frame-1'].edges

	mdb.models['standard'].ContactProperty('interior')
    	mdb.models['standard'].interactionProperties['interior'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)	



	startx = min(vex[0][0])
	
	first = False 
	alter = False 
	for i in range(0, len(vex[0][0])):
		x = vex[0][0][i]
		y = vex[0][1][i]
		coord = s1.getClosest(coordinates=((x,y,0),))
		if radiuschecker(startx, radius-thickness*0.4, x,y):
			
			if first==True:
				side1Edges1 += s1.findAt((coord[0][1],))
			else:	#Initializing side1Edges1 for iteration
				side1Edges1 = s1.findAt((coord[0][1],))
				first = True
		else:
			if alter==True:
				side1Edges2 += s1.findAt((coord[0][1],))
			else:	#Initializing side1Edges1 for iteration
				side1Edges2 = s1.findAt((coord[0][1],))
				alter = True

	region1=a.Surface(side1Edges=side1Edges1, name='polygonouterL')
	regionalt1 = a.Surface(side1Edges=side1Edges2, name='polygoninnerl')
	
	first = False 
	alter = False
	for i in range(0, len(vex[0][0])):
		x = vex[1][0][i]
		y = vex[1][1][i]
		coord = s1.getClosest(coordinates=((x,y,0),))
		if radiuschecker(-startx, radius-thickness*0.4, x,y):
			
			if first==True:
				side1Edges1 += s1.findAt((coord[0][1],))
			else: 
				side1Edges1 = s1.findAt((coord[0][1],))
				first=True
		else:
			if alter==True:
				side1Edges2 += s1.findAt((coord[0][1],))
			else:	#Initializing side1Edges1 for iteration
				side1Edges2 = s1.findAt((coord[0][1],))
				alter = True
	



	region2=a.Surface(side1Edges=side1Edges1, name='polygonouterr')
	regionalt2 = a.Surface(side1Edges=side1Edges2, name='polygoninnerr')
	
	mdb.models['standard'].SurfaceToSurfaceContactExp(name='polygon to polygon', 
        	createStepName='Apply load', master=region1, slave=region2, 
        	mechanicalConstraint=KINEMATIC, sliding=FINITE, 
		interactionProperty='interior', initialClearance=OMIT, datumAxis=None, 
		clearanceRegion=None)
	mdb.models['standard'].SelfContactExp(name='selfl', createStepName='Initial', 
        	surface=regionalt1, mechanicalConstraint=KINEMATIC, 
        	interactionProperty='interior')
	mdb.models['standard'].SelfContactExp(name='selfr', createStepName='Initial', 
        	surface=regionalt2, mechanicalConstraint=KINEMATIC, 
        	interactionProperty='interior')
	
	s2 = a.instances['planetop-1'].edges
	coord = s2.getClosest(coordinates=((0,radius,0),))
	side1Edges1 = s2.findAt((coord[0][1],))
	regiontop = a.Surface(side1Edges=side1Edges1, name='planetop-bottom')
	
	mdb.models['standard'].SurfaceToSurfaceContactExp(name='polygonL - planetop', 
        	createStepName='Apply load', master=regiontop, slave=region1, 
        	mechanicalConstraint=KINEMATIC, sliding=FINITE, 
		interactionProperty='interior', initialClearance=OMIT, datumAxis=None, 
		clearanceRegion=None)

	mdb.models['standard'].SurfaceToSurfaceContactExp(name='polygonR - planetop', 
        	createStepName='Apply load', master=regiontop, slave=region2, 
        	mechanicalConstraint=KINEMATIC, sliding=FINITE, 
		interactionProperty='interior', initialClearance=OMIT, datumAxis=None, 
		clearanceRegion=None)

	s2 = a.instances['plane-1'].edges
	coord = s2.getClosest(coordinates=((0,-radius,0),))
	side1Edges1 = s2.findAt((coord[0][1],))
	regionbot = a.Surface(side1Edges=side1Edges1, name='plane-top')

	mdb.models['standard'].SurfaceToSurfaceContactExp(name='polygonL - planebot', 
        	createStepName='Apply load', master=regionbot, slave=region1, 
		mechanicalConstraint=KINEMATIC, sliding=FINITE, 
		interactionProperty='interior', initialClearance=OMIT, datumAxis=None, 
		clearanceRegion=None)

	mdb.models['standard'].SurfaceToSurfaceContactExp(name='polygonR - planebot', 
        	createStepName='Apply load', master=regionbot, slave=region2, 
        	mechanicalConstraint=KINEMATIC, sliding=FINITE, 
		interactionProperty='interior', initialClearance=OMIT, datumAxis=None, 
		clearanceRegion=None)
	
	mdb.models['standard'].SurfaceToSurfaceContactExp(name='planetop-plane',
		createStepName='Apply load', master=regiontop, slave=regionbot,
		mechanicalConstraint=KINEMATIC, sliding=FINITE, 
		interactionProperty='interior', initialClearance=OMIT, datumAxis=None,
		clearanceRegion=None)
def partition(mdb, innerv, outerv):
#Input: model database, array of vertices for the inner polygon, array of vertices for the outer polygon
#Generates partition at each kink of polygon 
    from abaqus import *
    from abaqusConstants import *

    left = min(outerv[0])
    bottom = min(outerv[1])

    p = mdb.models['standard'].parts['Frame']
    f = p.faces
    pickedRegions = f.getByBoundingBox(1.02*left, 1.02*bottom,0,-1.02*left,-1.02*bottom,0)
    p.deleteMesh(regions=pickedRegions)
    p = mdb.models['standard'].parts['Frame']
    f, e, d = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
        0, 0.0, 0.0))
    s = mdb.models['standard'].ConstrainedSketch(name='__profile__', 
        sheetSize=4.16, gridSpacing=0.1, transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models['standard'].parts['Frame']
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
    for i in range(0, len(outerv[0])-1):
    	s.Line(point1=(innerv[0][i],innerv[1][i]), point2=(outerv[0][i],outerv[1][i]))
    p = mdb.models['standard'].parts['Frame']
    f = p.faces
    pickedFaces = f.getByBoundingBox(3*left, 3*bottom,0,-3*left,-3*bottom,0)
    e1, d2 = p.edges, p.datums
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models['standard'].sketches['__profile__']

def const_mass(radius, thickness, sides):
	import math
	
	p_sin = math.sin(180/sides)
	p_cos = math.cos(180/sides) 
	p_tan = p_sin/p_cos

	sidelen = sidelength(radius, sides)
	inner_p = sidelen*sides
	inner_a = sidelen/2/p_tan
	inner_area = inner_p*inner_a/2
    
	sidelen = sidelength(radius+thickness, sides)
	outer_p = sidelen*sides
	outer_a = sidelen/2/p_tan
	outer_area = outer_a*outer_p/2
    
	area = outer_area-inner_area 

	print "Auto-generated list of parameters for constant mass" 
	 
	n = 4
	while n <= 10:
		radicand = area/n/p_sin/p_cos+radius**2
        	if  radicand >0:
			thick = radicand**0.5-radius
              	else:
            		thick = -1 
        	print "%d-sides: Thickness: %f" % (n, thick) 
        	n += 2