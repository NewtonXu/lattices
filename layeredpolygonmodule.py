def line(s, left, top):
    """Input: s = mdb.models['model-name'].ConstrainedSketch()
              left = (f) min - x of the line
              top = (f) max - y of the line
              Function: Draws a line as specified above
       Output: none
    """
    s.Line(point1=(-2 * left, -top), point2=(2 * left, -top))


def rectangle(s, left, top, thickness):
    """Input: s = mdb.models['model-name'].ConstrainedSketch()
              left = (f) min - x of the rectangle
              top = (f) max - y of the rectangle
              thickness = (f) width desired for rectangle along y-axis
       Function: Draws a rectangle where the vertices are the two points
       Output: none
    """
    s.rectangle(point1=(left, top), point2=(-left, top+thickness))


def sketch(s, array):
    """Input: s = mdb.models['model-name'].ConstrainedSketch()
              array = (list) list[0=x,1=y][#=vertex number]
       Function: Connects all points with straight lines
       Output: none
    """
    for i in range(0, len(array[0]) - 1):
        s.Line(point1=(array[0][i], array[1][i]),
               point2=(array[0][i+ 1], array[1][i + 1]))


def advancedarrayspolygon(inradius, thickness, sides,distance):
    """Input: inradius = (f) outer radius of inner polygon
              thickness = (f) distance between inner and outer polygon
              sides = (int) number of sides of polygon
              distance = (f) distance between the halves of the polygon
       Function: Creates list with coordinates of polygon halves arranged so
                 inner polygon lines up with outer polygon
       Output: (list) list[0=inner,1=outer][0=x,1=y][#=vertex number]
    """
    innervertices = vertices(0,0, inradius, sides)
    outervertices = vertices(0,0, inradius+thickness, sides)
    
    arraylen = len(innervertices[0])
                     
    innervleft = [innervertices[0][0:arraylen//2],
                  innervertices[1][0:arraylen//2]]
    innervright = [innervertices[0][arraylen//2:arraylen-1],
                   innervertices[1][arraylen//2:arraylen-1]]

    outervleft = [outervertices[0][0:arraylen//2],
                  outervertices[1][0:arraylen//2]]
    outervright = [outervertices[0][arraylen//2:arraylen-1],
                   outervertices[1][arraylen//2:arraylen-1]]
    
    topnodeouter = max(outervertices[1])
    topnodeinner = max(innervertices[1])
    leftnodeouter = max(outervertices[1])
    
    innerleft = list(map(lambda x: x-leftnodeouter-distance, innervleft[0]))
    innerright = list(map(lambda x: x+leftnodeouter+distance, innervright[0]))
    innerv = [innerleft, innervleft[1]]
    innerv[0].extend(innerright)
    innerv[1].extend(innervright[1])
    
    outerleft = list(map(lambda x: x-leftnodeouter-distance, outervleft[0]))
    outerright = list(map(lambda x: x+leftnodeouter+distance, outervright[0]))
    outerv = [outerleft, outervleft[1]]
    outerv[0].extend(outerright)
    outerv[1].extend(outervright[1])
    return [innerv,outerv]


def advancedpolygon(inradius,thickness,sides,distance):
    """Input: inradius = (f) outer radius of inner polygon
              thickness = (f) distance between inner and outer polygon
              sides = (int) number of sides of polygon
              distance = (f) distance between the halves of the polygon
       Function: Creates a list of coordinates of the two polygon halves with
                 adjacent vertices in order.
       Output: (list) list[0=left,1=right][0=x,1=y][#=vertex number]
    """
    innervertices = vertices(0,0, inradius, sides)
    outervertices = vertices(0,0, inradius+thickness, sides)

    arraylen = len(innervertices[0])
                     
    innervleft = [innervertices[0][0:arraylen//2],
                  innervertices[1][0:arraylen//2]]
    innervright = [innervertices[0][arraylen//2:arraylen-1],
                   innervertices[1][arraylen//2:arraylen-1]]

    outervleft = [outervertices[0][0:arraylen//2],
                  outervertices[1][0:arraylen//2]]
    outervright = [outervertices[0][arraylen//2:arraylen-1],
                   outervertices[1][arraylen//2:arraylen-1]]

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

    transl_left = list(map(lambda x: x-leftnodeouter-distance, arrayleft[0]))
    arrayleftf = [transl_left, arrayleft[1]]
    
    arrayright = [[0,0],[topnodeouter,topnodeinner]]
    arrayright[0].extend(list(reversed(innervright[0])))
    arrayright[1].extend(list(reversed(innervright[1])))
    arrayright[0].extend([0,0])
    arrayright[1].extend([-topnodeinner,-topnodeouter])
    arrayright[0].extend(outervright[0])
    arrayright[1].extend(outervright[1])
    arrayright[0].append(0)
    arrayright[1].append(topnodeouter)

    transl_right = list(map(lambda x: x+leftnodeouter+distance, arrayright[0]))
    arrayrightf = [transl_right, arrayright[1]]
    
    return [arrayleftf, arrayrightf]


def vertices(startx,starty,radius,sides):
    """Input: startx = (f) x-coord of centre point
              starty = (f) y-coord of centre point
              radius = (f) distance of center to each polygon vertex
              sides = (int) number of polygon sides
       Function: Calculates vertices of the polygon described by the input
       Output: (list) list[0=x,1=y][#=vertex number]
    """
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

    xarray.append(xarray[0])
    yarray.append(yarray[0])
    return [xarray, yarray]


def trig(x):
    """Input: (f) x-radians
       Function: Calculates sin and cos of x
       Output: (tuple) (sin(x), cos(x))
    """
    import math 
    return (math.sin(x), math.cos(x))


def sidelength(radius, sides):
    """Input: radius = (f) radius of polygon
              sides = (int) number of polygon sides
       Funciton: Calculates distance between two adjacent vertices
       Output: (f) side length
    """
    import math 
    side = 2*radius*math.sin(math.pi/sides)
    return side


def rigidbody(mdb, vertices, thickness):
    """Input: mdb = (obj) model database
              vertices = (list) list of polygon vertices
              thickness = (f) thickness of polygon
       Function: Applies rigid body constraint by creating reference point
                 and selecting bottom plane
       Output: none
    """
    from abaqus import *
    from abaqusConstants import *
    import regionToolset
    output = []
    num_vert = len(vertices[0])
    a = mdb.models['standard'].rootAssembly

    bottom = min(vertices[0][1])
    left = min(vertices[0][0])
    
    a = mdb.models['standard'].rootAssembly
    e1 = a.instances['planetop-1'].edges
    coord = e1.getClosest(coordinates=((0,-bottom+thickness,0),))

    RP = a.ReferencePoint(point=a.instances['planetop-1'].InterestingPoint(
        edge=e1.findAt(coordinates=coord[0][1]), rule=MIDDLE))
    RP_id = RP.id
    region2=a.instances['planetop-1'].sets['planar2']
    r1 = a.referencePoints
    refPoints1 = mdb.models['standard'].rootAssembly.referencePoints[RP_id]
    region1=regionToolset.Region(referencePoints=(refPoints1,))
    mdb.models['standard'].RigidBody(
        name='Constraint-1', refPointRegion=region1, bodyRegion=region2)

	#Rigid body definition for bot plane
    a = mdb.models['standard'].rootAssembly
    e2 = a.instances['plane-1'].edges
    coord = e2.getClosest(coordinates=((0,bottom+thickness,0),))

    RP=a.ReferencePoint(point=a.instances['plane-1'].InterestingPoint(
        edge=e2.findAt(coordinates=coord[0][1]), rule=MIDDLE))
    RP_id = RP.id
    region2=a.instances['plane-1'].sets['planar']
    r1 = a.referencePoints
    refPoints1 = mdb.models['standard'].rootAssembly.referencePoints[RP_id]
    region1=regionToolset.Region(referencePoints=(refPoints1,))
    mdb.models['standard'].RigidBody(
        name='Constraint-2', refPointRegion=region1, bodyRegion=region2)

def bc_bot(mdb, vertices, thickness,layer=1, name="Frame-1"):
    """Input: mdb = (obj) model database
              vertices = (list) list of polygon vertices
              thickness = (f) thickness of polygon
       Function: Applies x-symmetry to sides of model and zero-displacement to
                 bottom face
       Output: none
    """
    from abaqus import *
    from abaqusConstants import *
    import regionToolset
    
    output = []
    num_vert = len(vertices[0])
    a = mdb.models['standard'].rootAssembly

    bottom = min(vertices[0][1])
    left = min(vertices[0][0])

    a = mdb.models['standard'].rootAssembly
    s2 = a.instances[name].edges
    s2a = a.instances[name].surfaces
    coord = s2.getClosest(coordinates=((0,bottom-thickness,0),))
    side1Edges1 = s2.findAt((coord[0][1],))
    coord = s2.getClosest(coordinates=((left+thickness,bottom-thickness,0),))
    side1Edges1 += s2.findAt((coord[0][1],))
    coord = s2.getClosest(coordinates=((-left-thickness,bottom-thickness,0),))
    side1Edges1 += s2.findAt((coord[0][1],))

    region2 = a.Set(edges=side1Edges1, name='plane-bot'+str(layer))
    region2a = a.Surface(side1Edges=side1Edges1, name='plane-bot'+str(layer))

    a = mdb.models['standard'].rootAssembly
    s2 = a.instances[name].edges
    coord = s2.getClosest(coordinates=((0,-bottom+thickness,0),))
    side1Edges1 = s2.findAt((coord[0][1],))
    coord = s2.getClosest(coordinates=((left+thickness,-bottom+thickness,0),))
    side1Edges1 += s2.findAt((coord[0][1],))
    coord = s2.getClosest(coordinates=((-left-thickness,-bottom+thickness,0),))
    side1Edges1 += s2.findAt((coord[0][1],))

    region3 = a.Set(edges=side1Edges1, name='plane-top'+str(layer))
    region3a = a.Surface(side1Edges=side1Edges1, name='plane-top'+str(layer))

    if layer ==1:
        mdb.models['standard'].DisplacementBC(name='planelock', 
            createStepName='Apply load', region=region2, u1=SET, u2=SET, ur3=SET, 
            amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
            localCsys=None)

    coord = s2.getClosest(coordinates=((-left,bottom+thickness/2,0),))
    side1Edges1 = s2.findAt((coord[0][1],))
    coord = s2.getClosest(coordinates=((left,bottom+thickness/2,0),))
    side1Edges1 += s2.findAt((coord[0][1],))
    coord = s2.getClosest(coordinates=((-left,-bottom-thickness/2,0),))
    side1Edges1 += s2.findAt((coord[0][1],))
    coord = s2.getClosest(coordinates=((left,-bottom-thickness/2,0),))
    side1Edges1 += s2.findAt((coord[0][1],))

    regionbot = a.Set(edges=side1Edges1, name='xsym'+str(layer))

    mdb.models['standard'].XsymmBC(name='xsym'+str(layer), createStepName='Initial', 
        region=regionbot, localCsys=None)



def loader(mdb, left, bottom, radius, thickness, planethick, pressure, name="Frame-1"):
    """Input: mdb = (obj) model database
              left = (f) min-x coord of polygon
              bottom = (f) min-y coord of polygon
              radius = (f) radius of polygon
              thickness = (f) thickness of polygon
              planethick = (f) thickness of the plane
              pressure = (f) pressure in MPA
       Function: Applies blast pressure to top plane (creates load+amplitude) 
                 Amplitude: (0.001 * t / 100, 3.4 * (1 - t / 100.0)
       Output: none
    """	
    from abaqus import *
    from abaqusConstants import *
    a = mdb.models['standard'].rootAssembly
    
    s2 = a.instances[name].edges
    coord = s2.getClosest(coordinates=((0,radius+thickness+planethick,0),))
    side1Edges1 = s2.findAt((coord[0][1],))
    coord = s2.getClosest(coordinates=((left+thickness, radius+thickness+planethick,0),))
    side1Edges1 += s2.findAt((coord[0][1],))
    coord = s2.getClosest(coordinates=((-left-thickness, radius+thickness+planethick,0),))
    side1Edges1 += s2.findAt((coord[0][1],))
    region = a.Surface(side1Edges=side1Edges1, name='planetop-top')

    #data = tuple((0.00001 * t / 100, 3.4 * (1 - t / 100.0)) for t in range(100))
    mdb.models['standard'].TabularAmplitude(name='Blast Amplitude', 
        timeSpan=TOTAL, smooth=SOLVER_DEFAULT, 
        data=pressure+((0.00101, 0.0), (0.002, 0.0), (100000, 0.0)))
    mdb.models['standard'].Pressure(name='Load 1', createStepName='Apply load', 
        region=region, magnitude=pressure[0][1], amplitude='Blast Amplitude') 

def radiuschecker(startx, maxrad, x, y):
    """Input: startx = (f) x-coordinate of center of polygon
              maxrad = (f) radius to check
              x = (f) x-coordinate of point to check
              y = (f) y-coordinate of point to check
       Function: Returns if (x,y( is more than maxrad away from (startx, 0)
       Output: (bool)
    """
    return maxrad<((x-startx)**2+y**2)**0.5


def edgecentre(vex):
    """Input: vex = (list) list of polygon vertices
       Function: Returns the midpoints of all the sides of the polygon
       Output: (list) list[0=left,1=right][0=x,1=y]
    """
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


def boxselector(mdb, vex, radius, thickness, distance,layer=1,name="Frame-1"):
    """Input: mdb = (obj) model database
              vex = (list) vertices of polygon
              radius = (f) radius of polygon
              thickness = (f) thickness of polygon
              distance = (f) distance between polygon halves
       Function: Use if 4-sided polygon of type 1
                 Selects edges on left halve, right halve, and interior
                 Creates self-contact for left set, interior set, right set
       Output: none
    """
    import math
    startx = min(vex[0][0])
    starty = max(vex[0][1])-math.cos(math.pi/4)*thickness

    vex = edgecentre(vex)
    from abaqus import *
    from abaqusConstants import *
    a = mdb.models['standard'].rootAssembly
    s1 = a.instances[name].edges

    mdb.models['standard'].ContactProperty('interior')
    mdb.models['standard'].interactionProperties['interior'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)

    coord = s1.getClosest(coordinates=((-distance-thickness-radius/2, -starty, 0),))
    regionleft = s1.findAt((coord[0][1],))
    coord = s1.getClosest(coordinates=((-distance-thickness-radius/2, starty,0),))
    regionleft += s1.findAt((coord[0][1],))
    coord = s1.getClosest(coordinates=((-distance-thickness-radius/2, 0,0),))
    regionleft+=s1.findAt((coord[0][1],))

    coord = s1.getClosest(coordinates=((distance+thickness+radius/2, -starty, 0),))
    regionright = s1.findAt((coord[0][1],))
    coord = s1.getClosest(coordinates=((distance+thickness+radius/2, starty,0),))
    regionright += s1.findAt((coord[0][1],))
    coord = s1.getClosest(coordinates=((distance+thickness+radius/2, 0,0),))
    regionright+=s1.findAt((coord[0][1],))

    coord = s1.getClosest(coordinates=((-distance,0,0),))
    regioninner = s1.findAt((coord[0][1],))
    coord = s1.getClosest(coordinates=((distance,0,0),))
    regioninner += s1.findAt((coord[0][1],))
    coord = s1.getClosest(coordinates=((0,starty,0),))
    regioninner += s1.findAt((coord[0][1],))
    coord = s1.getClosest(coordinates=((0,-starty,0),))
    regioninner += s1.findAt((coord[0][1],))

    regionalt1=a.Surface(side1Edges=regionleft, name='polygoninnerl'+str(layer))
    regionalt2 = a.Surface(side1Edges=regionright, name='polygoninnerr'+str(layer))
    regioninner = a.Surface(side1Edges=regioninner, name='polygoninner'+str(layer))
    mdb.models['standard'].SelfContactExp(name='inner'+str(layer), createStepName='Initial', 
        surface=regioninner, mechanicalConstraint=KINEMATIC, 
        interactionProperty='interior')
    mdb.models['standard'].SelfContactExp(name='selfl'+str(layer), createStepName='Initial', 
        surface=regionalt1, mechanicalConstraint=KINEMATIC, 
        interactionProperty='interior')
    mdb.models['standard'].SelfContactExp(name='selfr'+str(layer), createStepName='Initial', 
        surface=regionalt2, mechanicalConstraint=KINEMATIC, 
        interactionProperty='interior')

def edgeselector(mdb, vex, radius, thickness, sides, apothem,layer=1,name="Frame-1"):
    """Input: mdb = (obj) model database
              vex = (list) vertices of polygon
              radius = (f) outer radius of polygon
              thickness = (f) thickness of polygon
              sides = (int) number of polygon sides
              apothem = (f) inner radios of polygon
       Function: Use if not 4-sided polygon, works for type 0 and 1 polygon
                 Finds edges to create interaction sets
                 Creates interactions
       Output: none
    """
    from abaqus import *
    from abaqusConstants import *

    startx = min(vex[0][0])
    starty = apothem
    vex = edgecentre(vex)
    a = mdb.models['standard'].rootAssembly
    s1 = a.instances[name].edges

    mdb.models['standard'].ContactProperty('interior')
    mdb.models['standard'].interactionProperties['interior'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)	
    
    first = False 
    alter = False 
    for i in range(0, len(vex[0][0])):
        x = vex[0][0][i]
        y = vex[0][1][i]
        coord = s1.getClosest(coordinates=((x,y,0),))
        if radiuschecker(startx, apothem+thickness*0.4, x,y) and abs(y)<(apothem+thickness*0.4):
            if first==True:
                side1Edges1 += s1.findAt((coord[0][1],))
            else:	#Initializing side1Edges1 for iteration
                side1Edges1 = s1.findAt((coord[0][1],))
                first = True
        else:
            if abs(y)<(apothem+thickness*0.4):
                if alter==True:
                    side1Edges2 += s1.findAt((coord[0][1],))
                else:	#Initializing side1Edges1 for iteration
                    side1Edges2 = s1.findAt((coord[0][1],))
                    alter = True

    region1=a.Surface(side1Edges=side1Edges1, name='polygonouterL'+str(layer))
    regionalt1 = a.Surface(side1Edges=side1Edges2, name='polygoninnerl'+str(layer))
    regioninner = side1Edges1
    first = False 
    alter = False
    for i in range(0, len(vex[0][0])):
        x = vex[1][0][i]
        y = vex[1][1][i]
        coord = s1.getClosest(coordinates=((x,y,0),))
        if radiuschecker(-startx, apothem+thickness*0.4, x,y) and abs(y)<(apothem+thickness*0.4):
            if first==True:
                side1Edges1 += s1.findAt((coord[0][1],))
            else: 
                side1Edges1 = s1.findAt((coord[0][1],))
                first=True
        else:
            if abs(y)<(apothem+thickness*0.4):
                if alter==True:
                    side1Edges2 += s1.findAt((coord[0][1],))
                else:	#Initializing side1Edges1 for iteration
                    side1Edges2 = s1.findAt((coord[0][1],))
                    alter = True
    coord = s1.getClosest(coordinates=((0,starty,0),))
    regioninner += s1.findAt((coord[0][1],))
    coord = s1.getClosest(coordinates=((0,-starty,0),))
    regioninner += s1.findAt((coord[0][1],))
    regioninner += side1Edges1
    region2=a.Surface(side1Edges=side1Edges1, name='polygonouterr'+str(layer))
    regionalt2 = a.Surface(side1Edges=side1Edges2, name='polygoninnerr'+str(layer))
    regioninner = a.Surface(side1Edges=regioninner, name='polygoninner'+str(layer))
    mdb.models['standard'].SelfContactExp(name='inner'+str(layer), createStepName='Initial', 
        surface=regioninner, mechanicalConstraint=KINEMATIC, 
        interactionProperty='interior')
    mdb.models['standard'].SelfContactExp(name='selfl'+str(layer), createStepName='Initial', 
        surface=regionalt1, mechanicalConstraint=KINEMATIC, 
        interactionProperty='interior')
    mdb.models['standard'].SelfContactExp(name='selfr'+str(layer), createStepName='Initial', 
        surface=regionalt2, mechanicalConstraint=KINEMATIC, 
        interactionProperty='interior')    
def layerinteraction(mdb,layer):
    """Input: mdb = (obj) model database
              layer = (int) layer number
       Function: Creates a contact interaction between layer and layer - 1
       Output: none
    """
    from abaqus import *
    from abaqusConstants import *

    a = mdb.models['standard'].rootAssembly
    regiontop = a.surfaces['plane-top'+str(layer-1)]
    regionbot =  a.surfaces['plane-bot'+str(layer)] 
    mdb.models['standard'].SurfaceToSurfaceContactExp(name='plane-plane'+str(layer),
            createStepName='Apply load', master=regiontop, slave=regionbot,
            mechanicalConstraint=KINEMATIC, sliding=FINITE, 
            interactionProperty='interior', initialClearance=OMIT, datumAxis=None,
            clearanceRegion=None)

def partition(mdb, innerv, outerv, radius, distance, thickness):
    """Input: mdb = (obj) model database
              innerv = (list) list of vertices for inner polygon
              outerv = (list) list of vertices for outer polygon
              radius = (f) outer radius of polygon
              distance = (f) distance between two halves
              thickness = (f) thickness of polygon
       Function: Draws lines between the inner and outer polygon vertices to
                 create partitions for uniform mesh.
                 If polyontype=1, also partitions polygon area from plane.
       Output: none
    """
    from abaqus import *
    from abaqusConstants import *

    left = min(outerv[0])
    bottom = min(outerv[1])
    inbottom = min(innerv[1])
    
    p = mdb.models['standard'].parts['Frame']
    f = p.faces
    pickedRegions = f.getByBoundingBox(3*left, 3*bottom,0,-3*left,-3*bottom,0)
    p.deleteMesh(regions=pickedRegions)
    p= mdb.models['standard'].parts['Frame']
    f, e, d = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(
        0, 0.0, 0.0))
    s = mdb.models['standard'].ConstrainedSketch(name='__profile__', 
        sheetSize=4.16, gridSpacing=0.1, transform=t)
    g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=SUPERIMPOSE)
    p = mdb.models['standard'].parts['Frame']
    p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

    for i in range(0, len(outerv[0])):
    	if abs(outerv[1][i]) < -inbottom:
            	s.Line(point1=(innerv[0][i],innerv[1][i]), point2=(outerv[0][i],outerv[1][i]))
    s.Line(point1=(3*left,-inbottom),point2=(-3*left,-inbottom))
    s.Line(point1=(3*left,inbottom),point2=(-3*left,inbottom))

    p = mdb.models['standard'].parts['Frame']
    f = p.faces
    pickedFaces = f.getByBoundingBox(20*left, 20*bottom,0,-20*left,-20*bottom,0)
    e1, d2 = p.edges, p.datums
    p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del mdb.models['standard'].sketches['__profile__']
def constraint(mdb,i):
    """Input: mdb = (obj) model database
              i = (int) layer number 
       Function: Creates a tie constraint between i and i-1 layer
       Output: none
    """
    from abaqus import *
    from abaqusConstants import *

    a = mdb.models['standard'].rootAssembly
    region1=a.surfaces['plane-top'+str(i-1)]
    a = mdb.models['standard'].rootAssembly
    region2=a.surfaces['plane-bot'+str(i)]
    mdb.models['standard'].Tie(name='Constraint'+str(i), master=region1, slave=region2, 
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, 
        thickness=ON)


def const_distance(radius,thickness,sides,distance,thickness_list):
    #plane_length = 2*(cos(180/n)*(r+t)+d)
    import math
    distance_list = []
    plane_length = 2*(math.cos(math.pi/sides)*(radius+thickness)+distance)
    print(plane_length)
    for n in range(4,14,2):
	index = (n-4)/2
	t = thickness_list[index]
	dist_calc = plane_length/2-math.cos(math.pi/n)*(radius+t)
	distance_list.append(dist_calc)
	print("%d-sides: distance: %f" % (n, dist_calc))
    return distance_list
def iterate(radius,thickness,n,distance,mdb):
    import math
    total_area = polygonarea(mdb)
    length = 2.0*(math.cos(math.pi/n)*(radius+thickness)+distance)
    print("Total area:"+str(total_area)+"Length"+str(length))
    
    theta_1 = math.pi/2.0-math.pi/n
    theta_2 = math.pi/2.0-3.0*math.pi/n
    A = (n-2)*math.sin(math.pi/n)*math.cos(math.pi/n)
    B = 4.0*math.sin(math.pi/2-math.pi/n)
    C = 2*math.sin(theta_1)*(math.cos(theta_1)-(math.sin(theta_1)*(math.cos(theta_2)-math.cos(theta_1)))/(math.sin(theta_2)-math.sin(theta_1)))
    mod = float(math.cos(math.pi/n))
    a = A-C
    b = 2*A*radius+B*length/2.0
    c = -total_area 
    if a>0.001:
        	t = -b+math.sqrt(b**2.0-4.0*a*c)
                t = t/2.0/a
    else:
		t = total_area/(2.0*A*radius+B*length/2.0)
    count = 0 
    while count < 5:
        dist = length/float(2) - mod*(radius+t)
	area = (A+B*mod-C)*t**2+(2.0*A*radius+B*mod*radius+B*dist)*t
        l = 2*(mod*(radius+t)+dist)
        if error(total_area, length, area, l) > 0.00001:
        	a=A+B*mod-C
        	b = 2*A*radius+B*mod*radius+B*dist
        	c = -total_area 

		if a>0.001:
        		temp = -b+math.sqrt(b**2.0-4.0*a*c)
                	temp = temp/2.0/a
		else:
			temp = total_area/(2.0*A*radius+B*length/2.0)

        	dist_2 = length/float(2) - mod*(radius+temp)
	
		area_2 = (A+B*mod-C)*temp**2+(2.0*A*radius+B*mod*radius+B*dist_2)*temp
        	l_2  = 2*(mod*(radius+temp)+dist_2)
       
                if error(total_area, length, area_2, l_2) > 0.00001:
                        t = temp 
			dist = dist_2
		else:
                        count = 6
	else:
            count = 6
	count += 1 
    return t,dist
def parameters_2(radius,thickness,sides,distance,mdb): 
    import math
    total_area = polygonarea(mdb)
    length = 2.0*(math.cos(math.pi/sides)*(radius+thickness)+distance)
    print("Total area:"+str(total_area)+"Length"+str(length))
    for n in range(4,20,2):
        theta_1 = math.pi/2.0-math.pi/n
    	theta_2 = math.pi/2.0-3.0*math.pi/n
    	A = (n-2)*math.sin(math.pi/n)*math.cos(math.pi/n)
    	B = 4.0*math.sin(math.pi/2-math.pi/n)
    	C = 2*math.sin(theta_1)*(math.cos(theta_1)-(math.sin(theta_1)*(math.cos(theta_2)-math.cos(theta_1)))/(math.sin(theta_2)-math.sin(theta_1)))
    	mod = float(math.cos(math.pi/n))
    	a = A-C
    	b = 2*A*radius+B*length/2.0
    	c = -total_area 
    	if a>0.001:
        	t = -b+math.sqrt(b**2.0-4.0*a*c)
                t = t/2.0/a
    	else:
		t = total_area/(2.0*A*radius+B*length/2.0)
    	count = 0 
        while count < 5:
        	dist = length/float(2) - mod*(radius+t)
		area = (A+B*mod-C)*t**2+(2.0*A*radius+B*mod*radius+B*dist)*t
        	l = 2*(mod*(radius+t)+dist)
        	if error(total_area, length, area, l) > 0.00001:
        		a=A+B*mod-C
        		b = 2*A*radius+B*mod*radius+B*dist
        		c = -total_area 

			if a>0.001:
        			temp = -b+math.sqrt(b**2.0-4.0*a*c)
                		temp = temp/2.0/a
			else:
				temp = total_area/(2.0*A*radius+B*length/2.0)

        		dist_2 = length/float(2) - mod*(radius+temp)
	
			area_2 = (A+B*mod-C)*temp**2+(2.0*A*radius+B*mod*radius+B*dist_2)*temp
        		l_2  = 2*(mod*(radius+temp)+dist_2)

       			t = temp 
			dist = dist_2
                	if error(total_area, length, area_2, l_2) < 0.00001:
                        	count = 10
		else:
            		count = 10
		count += 1 
	if count != 10:
            print("%d-sides did not converge, error exceeds preset parameter" % (n))
	out_area = (A+B*mod-C)*t**2+(2*A*radius+B*mod*radius+B*dist)*t
	out_length = 2*(mod*(radius+t)+dist)
        print("%d-sides: Thickness: %f Length: %f Area: %f Length: %f" % (n, t, dist,out_area,out_length)) 

def parameters(radius,thickness,sides,distance,mdb):
    import math
    total_area = polygonarea(mdb)
    length = 2.0*(math.cos(math.pi/sides)*(radius+thickness)+distance)
    print("Total area:"+str(total_area)+"Length"+str(length))
    for n in range(4,20,2):
        theta_1 = math.pi/2.0-math.pi/n
        theta_2 = math.pi/2.0-3.0*math.pi/n
        A = (n-2)*math.sin(math.pi/n)*math.cos(math.pi/n)
	B = 4.0*math.sin(math.pi/2-math.pi/n)
        C = 2*math.sin(theta_1)*(math.cos(theta_1)-(math.sin(theta_1)*(math.cos(theta_2)-math.cos(theta_1)))/(math.sin(theta_2)-math.sin(theta_1)))
        mod = float(math.cos(math.pi/n))

        a = A-C
        b = 2*A*radius+B*length/2.0
        c = -total_area 
	if a>0.001:
        	t = -b+math.sqrt(b**2.0-4.0*a*c)
                t = t/2.0/a
	else:
		t = total_area/(2.0*A*radius+B*length/2.0)
        dist = length/float(2) - mod*(radius+t)
	area = (A+B*mod-C)*t**2+(2.0*A*radius+B*mod*radius+B*dist)*t
        l = 2*(mod*(radius+t)+dist)
	#recheck t 
        
        a=A+B*mod-C
        b = 2*A*radius+B*mod*radius+B*dist
        c = -total_area 

	if a>0.001:
        	temp = -b+math.sqrt(b**2.0-4.0*a*c)
                temp = temp/2.0/a
	else:
		temp = total_area/(2.0*A*radius+B*length/2.0)

        dist_2 = length/float(2) - mod*(radius+temp)
	
	area_2 = (A+B*mod-C)*temp**2+(2.0*A*radius+B*mod*radius+B*dist_2)*temp
        l_2  = 2*(mod*(radius+temp)+dist_2)
	
	print("%d-sides: Thickness: %f Length: %f Area: %f Length: %f T2: %f D2: %f" % (n, t, dist,area,l,temp,dist_2)) 
	out_area = (A+B*mod-C)*temp**2+(2*A*radius+B*mod*radius+B*dist)*temp
	out_length = 2*(mod*(radius+temp)+dist)

	print("%d-sides: Thickness: %f Distance: %f Area: %f Length: %f" % (n,temp,dist,out_area,out_length))
def error(area, length, area_2, length_2):
    error_1 = abs(area_2-area)/area
    error_2 = abs(length_2-length)/length
    return(max(error_1,error_2))
    
def formula(radius,thickness,n,distance,mdb):
    import math
    print("Actual area:"+str(polygonarea(mdb)))
    theta_1 = math.pi/2-math.pi/n
    theta_2 = math.pi/2-3*math.pi/n
    
    A = (n-2)*math.sin(math.pi/n)*math.cos(math.pi/n)
    B = 4*math.sin(math.pi/2-math.pi/n)
    C = 2*math.sin(theta_1)*(math.cos(theta_1)-(math.sin(theta_1)*(math.cos(theta_2)-math.cos(theta_1)))/(math.sin(theta_2)-math.sin(theta_1)))
    m = math.cos(math.pi/n) 

    area = (A+B*m-C)*thickness**2+(2*A*radius+B*m*radius+B*distance)*thickness
    print("Calculated area:"+str(area))

def sidearea(radius,thickness,sides):
    """Input: radius = (f) radius of inner polygon
              thickness = (f) distance between inner and outer polygon
              sides = (int) number of sides of polygon
       Function: Calculates the area of one side of the polygon defined by input
       Output: (f) Area of one side of the polygon
    """
    import math

    p_sin = math.sin(math.pi/sides)
    p_cos = math.cos(math.pi/sides) 
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
    sidearea = area/sides

    return sidearea

def polygonarea(mdb):
    a = mdb.models["standard"].rootAssembly
    return a.getArea(a.instances["Frame-1"].faces)

def const_mass(radius, thickness, sides,distance):
    """Input: radius = (f) outer radius of polygon
              thickness = (f) thickness of polygon
              sides = (int) number of polygon sides
       Function: Calculates for even sided polygons with sides 4 to 10 the
                 thickness value to maintain equal area throughout, and
                 prints the values to console.
                 Note: Uniform area = Uniform mass for the 
		 plane-strain simulation.
		 **This is only for uniform polygon area, for uniform simulation
		   area, use the mass function 
       Output: none  (prints to console) 
    """
    import math
    
    p_sin = math.sin(math.pi/sides)
    p_cos = math.cos(math.pi/sides) 
    p_tan = p_sin/p_cos
    sidelen = sidelength(radius+thickness,sides)
    area_side = sidearea(radius,thickness,sides)
    area = area_side*sides
    p_sin = math.sin(math.pi/2-math.pi/sides)
    p_cos = math.cos(math.pi/2-math.pi/sides)
    if sides%4!=0:
         toplength = radius+thickness
    else:
         toplength = (radius+thickness)*math.cos(math.pi/sides)
    
    if sides ==4:
        modifier = 1
    else:
        modifier = 2 
    #flangearea = thickness*p_sin*(2*toplength+2*distance)-area_side-modifier*thickness*p_cos*thickness*p_sin
    plane_area = thickness*math.sin(math.pi/2-math.pi/sides)*(2*toplength+2*distance) 
    triangles = modifier*thickness*math.sin(math.pi/2-math.pi/sides)*thickness*math.cos(math.pi/2-math.pi/sides)
    total_area = area + 2*(plane_area)-2*area_side-4*triangles 
    print(total_area)
    
    for n in range(4,14,2):
        if n == 4:
            modifier = 1
        else:
            modifier = 2
            
        p_sin = math.sin(math.pi/n)
        p_cos = math.cos(math.pi/n)
        p_tan = p_sin/p_cos

        radicand = area/n/p_sin/p_cos+radius**2    #term under sqrt
	
        if  radicand >0:
            thick = radicand**0.5-radius
        else:
            thick = -1 
        temp_side_len = sidelength(radius+thick,sides)
        p_sin = math.sin(math.pi/2-math.pi/n)
        p_cos = math.cos(math.pi/2-math.pi/n)
        
        if n%4!=0:
            toplength = radius+thick
        else:
            toplength = (radius+thick)*math.cos(math.pi/n)

        
        #dist_calc = (flangearea + modifier*thick*p_cos*thick*p_sin+sidearea(radius,thick,n))/(2*thick*p_sin)-toplength
        #dist_calc = (total_area-n*sidearea(radius,thick,n)+2*sidearea(radius,thick,n)+4*modifier*thick*p_cos*thick*p_sin)/(4*thickness*math.sin(math.pi/2-math.pi/n))-0.5*toplength
	print("%d-sides: Thickness: %f" % (n, thick))  
def getrf(session, layers, odb_name):
#odb_name = "8-sides-smallblast.odb" as str 
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import time
    from abaqus import *
    from abaqusConstants import *

    o1 = session.openOdb(name='/mnt/compute-0-4/newton/'+odb_name)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    odb = session.odbs['/mnt/compute-0-4/newton/'+odb_name]
    for key in session.xyDataObjects.keys():
	del session.xyDataObjects[key]
    
    #Extract reaction force 
    if layers == 1:
    	session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
        	NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=('PLANE-BOT1', ))
    else:
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
        	NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=('PLANE-BOT1', ))
    array = session.xyDataObjects.values()
    array = tuple(array)
    xy121 = sum(array)
    xy121.setValues(
        sourceDescription="same")
    tmpName = xy121.name
    session.xyDataObjects.changeKey(tmpName, 'XYData-1')
    x0 = session.xyDataObjects['XYData-1']
   
#    savetime = str(time.time())
    session.xyReportOptions.setValues(numberFormat=AUTOMATIC)
    session.writeXYReport(fileName='rf'+odb_name+'.dat', xyData=(x0, ))
    print('rf'+odb_name+'.dat')
def getu(session, layers, odb_name, type=2):
#odb_name = "8-sides-smallblast.odb" as str 
#type = 0: min, type = 1:avg,, type = 2:max
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import time
    from abaqus import *
    from abaqusConstants import *
    
    savetime = str(time.time())
    o1 = session.openOdb(name='/mnt/compute-0-4/newton/'+odb_name)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    odb = session.odbs['/mnt/compute-0-4/newton/'+odb_name]
    
    outputset = ()
    for i in range(1,layers+1):
        for key in session.xyDataObjects.keys():
		del session.xyDataObjects[key]
        if layers == 1:
            nodeset=('PLANE-TOP1',)
        else:
            nodeset = ('PLANE-TOP'+str(i),)
        session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
            NODAL, ((COMPONENT, 'U2'), )), ), nodeSets = nodeset)
    	array = session.xyDataObjects.values()
    	array = tuple(array)
	if type == 0 or type=="min":
    	    xy121 = min(array)
	    type = "min"
	elif type==1 or type=="avg": 
            xy121 = avg(array)
	    type = "avg"
	else:
	    xy121 = max(array)
	    type = "max"
    	xy121.setValues(
            sourceDescription="same")
    	tmpName = xy121.name
    	session.xyDataObjects.changeKey(tmpName, 'XYData-'+str(i))
	outputset += (session.xyDataObjects['XYData-'+str(i)],)
    
    
    	session.xyReportOptions.setValues(numberFormat=AUTOMATIC)
    	
    	session.writeXYReport(fileName='U'+type+odb_name+str(i)+'.dat', xyData=outputset)
        #Be advised, i = 1 is the bottom-most layer
    	print('U'+type+odb_name+str(i)+'.dat')

def getse(session, odb_name):
#odb_name = "8-sides-smallblast.odb" as str 
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import time
    from abaqus import *
    from abaqusConstants import *

    o1 = session.openOdb(name='/mnt/compute-0-4/newton/'+odb_name)
    session.viewports['Viewport: 1'].setValues(displayedObject=o1)
    odb = session.odbs['/mnt/compute-0-4/newton/'+odb_name]
    for key in session.xyDataObjects.keys():
	del session.xyDataObjects[key]
    
    #Extract strain energy 
    session.XYDataFromHistory(name='XYData-1', odb=odb, 
        outputVariableName='Strain energy: ALLSE for Whole Model', steps=(
        'Apply load', ), )
    array = session.xyDataObjects.values()
    array = tuple(array)
    xy121 = sum(array)
    xy121.setValues(
        sourceDescription="same")
    #tmpName = xy121.name
    #session.xyDataObjects.changeKey(tmpName, 'XYData-1')
    x0 = session.xyDataObjects['XYData-1']
   
#    savetime = str(time.time())
    session.xyReportOptions.setValues(numberFormat=AUTOMATIC)
    session.writeXYReport(fileName='se'+odb_name+'.dat', xyData=(x0, ))
    print('se'+odb_name+'.dat')
