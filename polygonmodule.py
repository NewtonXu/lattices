def line(s, left, top):

    s.Line(point1=(-2 * left, -top), point2=(2 * left, -top))



def sketch(s, array):

# Given 2D array of xy coordinates, sketches polygon using shapes
# ....a = mdb.models['standard'].rootAssembly
# ....s = mdb.models['standard'].ConstrainedSketch(name='__profile__', sheetSize=4.0)

    for i in range(0, len(array[0]) - 1):
        s.Line(point1=(array[0][i], array[1][i]), point2=(array[0][i
               + 1], array[1][i + 1]))
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
	sincos = trig(180/sides)
	side = 2*radius*sincos[0]/sincos[1]
	return side

def bc_bot(mdb, vertices, thickness):
#Given array of vertices, applies encastre BC to the bottom most point or points 
	from abaqus import *
	from abaqusConstants import *
	output = []
	num_vert = len(vertices[0])
	a = mdb.models['standard'].rootAssembly
	v = a.instances['Frame-1'].vertices
#	if num_vert%2 != 0:	#Odd number of sides so two points on the bottom need to be fixed
#		region=(v.findAt(((vertices[0][num_vert/2], vertices[1][num_vert/2], 0.0), ), ), None, None, None)
#		mdb.models['standard'].EncastreBC(name='Fixedright', createStepName='Initial', region=region)
		
#		region=(v.findAt(((vertices[0][num_vert/2-1], vertices[1][num_vert/2-1], 0.0), ), ), None, None, None)
#		mdb.models['standard'].EncastreBC(name='Fixedleft', createStepName='Initial', region=region)
#	else:
#		region=(v.findAt(((vertices[0][num_vert/2], vertices[1][num_vert/2], 0.0), ), ), None, None, None)
#		mdb.models['standard'].EncastreBC(name='Fixedbottom', createStepName='Initial', region=region)
	
	bottom = min(vertices[0][1])
	left = min(vertices[0][0])

#Lock down the plane 
	a = mdb.models['standard'].rootAssembly
	region = a.instances['plane-1'].sets['planar']
    	mdb.models['standard'].DisplacementBC(name='planelock', 
        createStepName='Apply load', region=region, u1=SET, u2=SET, ur3=SET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)

#Constrain the edges
	s1 = a.instances['Frame-1'].edges
    	side2Edges1 = s1.getByBoundingSphere((0,0,0),-1*left)
    	region=a.Surface(side2Edges=side2Edges1, name='Surf-2')
	mdb.models['standard'].ContactProperty('interior')
    	mdb.models['standard'].interactionProperties['interior'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)

    	mdb.models['standard'].SelfContactStd(name='interiorcontact', 
        createStepName='Apply load', surface=region, 
        interactionProperty='interior', enforcement=NODE_TO_SURFACE, 
        thickness=OFF, smooth=0.2) 

#Contact between plane and polygon
	a = mdb.models['standard'].rootAssembly
    	s1 = a.instances['plane-1'].edges
    	side1Edges1 = s1.getByBoundingBox(3*left,1.03*bottom,0,-3*left,0.99*bottom,0)
    	region1=a.Surface(side2Edges=side1Edges1, name='plane')
    	a = mdb.models['standard'].rootAssembly
    	region2=a.instances['Frame-1'].sets['everything']
    	mdb.models['standard'].SurfaceToSurfaceContactStd(name='planeandframe', 
        createStepName='Apply load', master=region1, slave=region2, 
        sliding=FINITE, enforcement=NODE_TO_SURFACE, thickness=OFF, 
        interactionProperty='interior', surfaceSmoothing=NONE, 
        adjustMethod=NONE, smooth=0.2, initialClearance=OMIT, datumAxis=None, 
        clearanceRegion=None)

#Contact between planetop and polygon
	a = mdb.models['standard'].rootAssembly
    	s1 = a.instances['planetop-1'].edges
    	side1Edges1 = s1.getByBoundingBox(3*left,-0.99*bottom,0,-3*left,-1.03*bottom,0)
    	region1=a.Surface(side1Edges=side1Edges1, name='planetop')

    	a = mdb.models['standard'].rootAssembly
    	region2=a.instances['Frame-1'].sets['everything']
    	mdb.models['standard'].SurfaceToSurfaceContactStd(name='planetopandframe', 
        createStepName='Apply load', master=region1, slave=region2, 
        sliding=FINITE, enforcement=NODE_TO_SURFACE, thickness=OFF, 
        interactionProperty='interior', surfaceSmoothing=NONE, 
        adjustMethod=NONE, smooth=0.2, initialClearance=OMIT, datumAxis=None, 
        clearanceRegion=None)
#Generate sets for top and bottom edge
	a = mdb.models['standard'].rootAssembly
    	e1 = a.instances['Frame-1'].edges
        edges1 = e1.getByBoundingBox(3*left, -0.99*bottom, 0, -3*left, -1.02*bottom,0)
	a.Set(edges=edges1, name='topedge')

	a = mdb.models['standard'].rootAssembly
    	e1 = a.instances['Frame-1'].edges
        edges1 = e1.getByBoundingBox(3*left, 1.02*bottom, 0, -3*left, 0.99*bottom,0)
	a.Set(edges=edges1, name='botedge')

#Generate overclosure for top and bottom edge 
	regionDef=mdb.models['standard'].rootAssembly.sets['topedge']
        mdb.models['standard'].interactions['planetopandframe'].setValues(
        initialClearance=OMIT, adjustMethod=SET, sliding=FINITE, 
        enforcement=NODE_TO_SURFACE, thickness=OFF, 
        supplementaryContact=SELECTIVE, smooth=0.2, tied=ON, 
        adjustSet=regionDef, bondingSet=None)
        regionDef=mdb.models['standard'].rootAssembly.sets['botedge']
        mdb.models['standard'].interactions['planeandframe'].setValues(
        initialClearance=OMIT, adjustMethod=SET, sliding=FINITE, 
        enforcement=NODE_TO_SURFACE, thickness=OFF, 
        supplementaryContact=SELECTIVE, smooth=0.2, tied=ON, 
        adjustSet=regionDef, bondingSet=None)

#Generate sets for left and right edge 
	a = mdb.models['standard'].rootAssembly
    	e1 = a.instances['Frame-1'].edges
        edges1 = e1.getByBoundingBox(1.02*left, 3*bottom, 0, 0.98*left, -3*bottom,0)
	a.Set(edges=edges1, name='leftedge')

    	e1 = a.instances['Frame-1'].edges
        edges1 = e1.getByBoundingBox(-0.98*left, 3*bottom, 0, -1.02*left, -3*bottom,0)
	a.Set(edges=edges1, name='rightedge')

#Mirror symmetry along x-direction
	region = a.sets['leftedge']
    	mdb.models['standard'].XsymmBC(name='leftedge', createStepName='Apply load', 
        region=region, localCsys=None)

	region=a.sets['rightedge']
	mdb.models['standard'].XsymmBC(name='rightedge', createStepName='Apply load', 
        region=region, localCsys=None)
#Overclosure initial protection
	mdb.models['standard'].StdInitialization(name='CInit-1')


	
 	 
  
def loader(mdb, vex, force=0, velocity=False, velx=0, vely=0, velr3=0, time=2.0, maxinc=10000, initinc=0.01, minimum=2e-05, maximum=0.0125):
#Must give at least one argument for this to work 
	from abaqus import *
	from abaqusConstants import *	
	a = mdb.models['standard'].rootAssembly
#	v = a.instances['Frame-1'].vertices

	left = min(vex[0][0])
	bottom = min(vex[0][1])
	if velocity:
		#Create a step with increments 
		mdb.models['standard'].StaticStep(name='Apply load', previous='Initial', description='Description', timePeriod=time, maxNumInc=maxinc, initialInc=initinc, minInc=minimum, maxInc=maximum)
		session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Apply load')
		region = a.instances['planetop-1'].sets['planar2']
       		mdb.models['standard'].VelocityBC(name='velocity', createStepName='Apply load', 
        	region=region, v1=velx, v2=vely, vr3=velr3, amplitude=UNSET, 
        	localCsys=None, distributionType=UNIFORM, fieldName='')

#	if force != 0: 
#		mdb.models['standard'].StaticLinearPerturbationStep(name='Apply load', 
#    		previous='Initial', description='10kN central load', 
#    		matrixSolver=SOLVER_DEFAULT)
#		session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Apply load')
#		mdb.models['standard'].fieldOutputRequests['F-Output-1'].setValues(
#   		variables=PRESELECT, region=MODEL)
#		
#		region=((v.findAt(((vertices[0][0], vertices[1][0], 0.0), ), ), ), )
#		mdb.models['standard'].ConcentratedForce(name='Force', 
 #   		createStepName='Apply load', 
  #  		region=region, cf2=-force)	
def radiuschecker(startx, maxrad, x, y):
#Assumed that starty = zero 
#Returns true if x,y is more than maxrad away from (startx,0) 
	return maxrad<((x-startx)**2+y**2)**0.5
def edgeselector(mdb, vex, radius):
#Generates the interacting sets for the contact between the two polygons
#Then creates the contact interaction between the two 
	a = mdb.models['standard'].rootAssembly
	s1 = a.instances['Frame-1'].edges
	startx = min(vex[0][0])
	array = [] #temporary
	for i in range(0, len(vex[0][0]-1)):
		x = vex[0][0][i]
		y = vex[0][1][i]
		if radiuschecker(startx, radius, x,y):
			coord = s1.getClosest(coordinates=((x,y,0),))
			array += coord[0][1] #change this to the findat function when possible 
			#Name this set left or right or whatever
	
	array = []
	for i in range(0, len(vex[0][0]-1)):
		x = vex[1][0][i]
		y = vex[1][1][i]
		if radiuschecker(-startx, radius, x,y):	#-startx for centre point on the right region 
			coord = s1.getClosest(coordinates=((x,y,0),))
			array += coord[0][1] #change this to the findat function when possible 
			#Name this set left or right or whatever



