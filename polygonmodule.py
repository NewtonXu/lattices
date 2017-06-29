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


def bc_bot(polygontype, mdb, vertices, thickness):
    """Input: polygontype = (int) 0 = discrete plane, 1 = integrated plane
              mdb = (obj) model database
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
    if polygontype == 0:
        a = mdb.models['standard'].rootAssembly
        s2 = a.instances['plane-1'].edges
        coord = s2.getClosest(coordinates=((0,bottom-thickness,0),))
        side1Edges1 = s2.findAt((coord[0][1],))
        region2 = a.Set(edges=side1Edges1, name='plane-bot')

        mdb.models['standard'].DisplacementBC(name='planelock', 
            createStepName='Apply load', region=region2, u1=SET, u2=SET, ur3=SET, 
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
    else:
        a = mdb.models['standard'].rootAssembly
        s2 = a.instances['Frame-1'].edges
        coord = s2.getClosest(coordinates=((0,bottom-thickness,0),))
        side1Edges1 = s2.findAt((coord[0][1],))
        coord = s2.getClosest(coordinates=((left+thickness,bottom-thickness,0),))
        side1Edges1 += s2.findAt((coord[0][1],))
        coord = s2.getClosest(coordinates=((-left-thickness,bottom-thickness,0),))
        side1Edges1 += s2.findAt((coord[0][1],))

        region2 = a.Set(edges=side1Edges1, name='plane-bot')

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

        regionbot = a.Set(edges=side1Edges1, name='xsym')

        mdb.models['standard'].XsymmBC(name='xsym', createStepName='Initial', 
            region=regionbot, localCsys=None)


def loader(polygontype, mdb, left, bottom, radius, thickness, planethick, pressure):
    """Input: polygontype = (int) 0 = discrete plane, 1 = integrated plane
              mdb = (obj) model database
              left = (f) min-x coord of polygon
              bottom = (f) min-y coord of polygon
              radius = (f) radius of polygon
              thickness = (f) thickness of polygon
              planethick = (f) thickness of the plane
                           Note if polygontype == 1, this equals thickness
              pressure = (f) pressure in MPA
       Function: Applies blast pressure to top plane (creates load+amplitude) 
                 Amplitude: (0.00001 * t / 100, 3.4 * (1 - t / 100.0)
       Output: none
    """	
    from abaqus import *
    from abaqusConstants import *
    a = mdb.models['standard'].rootAssembly
    if polygontype ==0: 
        s2 = a.instances['planetop-1'].edges
        coord = s2.getClosest(coordinates=((0,radius+thickness+planethick,0),))
        side1Edges1 = s2.findAt((coord[0][1],))
        region = a.Surface(side1Edges=side1Edges1, name='planetop-top')
    else:
        s2 = a.instances['Frame-1'].edges
        coord = s2.getClosest(coordinates=((0,radius+thickness+planethick,0),))
        side1Edges1 = s2.findAt((coord[0][1],))
        coord = s2.getClosest(coordinates=((left+thickness, radius+thickness+planethick,0),))
        side1Edges1 += s2.findAt((coord[0][1],))
        coord = s2.getClosest(coordinates=((-left-thickness, radius+thickness+planethick,0),))
        side1Edges1 += s2.findAt((coord[0][1],))
        region = a.Surface(side1Edges=side1Edges1, name='planetop-top')

    data = tuple((0.001 * t / 100, 3.4 * (1 - t / 100.0)) for t in range(100))
    mdb.models['standard'].TabularAmplitude(name='Blast Amplitude', 
        timeSpan=TOTAL, smooth=SOLVER_DEFAULT, 
        data=data+((0.00101, 0.0), (0.002, 0.0), (100000, 0.0)))
    mdb.models['standard'].Pressure(name='Load 1', createStepName='Apply load', 
        region=region, magnitude=3.4, amplitude='Blast Amplitude') 

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


def boxselector(mdb, vex, radius, thickness, distance):
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
    startx = min(vex[0][0])
    starty = max(vex[0][1])-2*thickness

    vex = edgecentre(vex)
    from abaqus import *
    from abaqusConstants import *
    a = mdb.models['standard'].rootAssembly
    s1 = a.instances['Frame-1'].edges

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

    regionalt1=a.Surface(side1Edges=regionleft, name='polygoninnerl')
    regionalt2 = a.Surface(side1Edges=regionright, name='polygoninnerr')
    regioninner = a.Surface(side1Edges=regioninner, name='polygoninner')
    mdb.models['standard'].SelfContactExp(name='inner', createStepName='Initial', 
        surface=regioninner, mechanicalConstraint=KINEMATIC, 
        interactionProperty='interior')
    mdb.models['standard'].SelfContactExp(name='selfl', createStepName='Initial', 
        surface=regionalt1, mechanicalConstraint=KINEMATIC, 
        interactionProperty='interior')
    mdb.models['standard'].SelfContactExp(name='selfr', createStepName='Initial', 
        surface=regionalt2, mechanicalConstraint=KINEMATIC, 
        interactionProperty='interior')


def edgeselector(polygontype, mdb, vex, radius, thickness, sides, apothem):
    """Input: polygontype = (int) 0=discrete plane, 1=integrated plane
              mdb = (obj) model database
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
    startx = min(vex[0][0])
    starty = apothem

    vex = edgecentre(vex)
    from abaqus import *
    from abaqusConstants import *
    a = mdb.models['standard'].rootAssembly
    s1 = a.instances['Frame-1'].edges

    mdb.models['standard'].ContactProperty('interior')
    mdb.models['standard'].interactionProperties['interior'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)	
    
    if polygontype ==0:
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
                else:#Initializing side1Edges1 for iteration
                    side1Edges2 = s1.findAt((coord[0][1],))
                    alter = True

    else:
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

        region1=a.Surface(side1Edges=side1Edges1, name='polygonouterL')
        regionalt1 = a.Surface(side1Edges=side1Edges2, name='polygoninnerl')
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
        region2=a.Surface(side1Edges=side1Edges1, name='polygonouterr')
        regionalt2 = a.Surface(side1Edges=side1Edges2, name='polygoninnerr')
        regioninner = a.Surface(side1Edges=regioninner, name='polygoninner')
        mdb.models['standard'].SelfContactExp(name='inner', createStepName='Initial', 
            surface=regioninner, mechanicalConstraint=KINEMATIC, 
            interactionProperty='interior')
        mdb.models['standard'].SelfContactExp(name='selfl', createStepName='Initial', 
            surface=regionalt1, mechanicalConstraint=KINEMATIC, 
            interactionProperty='interior')
        mdb.models['standard'].SelfContactExp(name='selfr', createStepName='Initial', 
            surface=regionalt2, mechanicalConstraint=KINEMATIC, 
            interactionProperty='interior')

    if polygontype == 0:
        mdb.models['standard'].SurfaceToSurfaceContactExp(name='polygon to polygon', 
            createStepName='Apply load', master=region1, slave=region2, 
            mechanicalConstraint=KINEMATIC, sliding=FINITE, 
            interactionProperty='interior', initialClearance=OMIT, 
            datumAxis=None, clearanceRegion=None)
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
            interactionProperty='interior', initialClearance=OMIT, 
            datumAxis=None, clearanceRegion=None)

        coord = s2.getClosest(coordinates=((0,-radius,0),))
        side1Edges1 = s2.findAt((coord[0][1],))
        regionbot = a.Surface(side1Edges=side1Edges1, name='plane-top')

        mdb.models['standard'].SurfaceToSurfaceContactExp(name='polygonL - planebot', 
            createStepName='Apply load', master=regionbot, slave=region1, 
            mechanicalConstraint=KINEMATIC, sliding=FINITE, 
            interactionProperty='interior', initialClearance=OMIT, 
            datumAxis=None, clearanceRegion=None)

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


def partition(polygontype, mdb, innerv, outerv, radius, distance, thickness):
    """Input: polygontype = (int) 0=discrete plane, 1=integrated plane
              mdb = (obj) model database
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
    
    if polygontype==0:
        p = mdb.models['standard'].parts['Frame']
        f = p.faces
        pickedRegions = f.getByBoundingBox(1.02*left, 1.02*bottom,0,-1.02*left,-1.02*bottom,0)
        p.deleteMesh(regions=pickedRegions)
        p = mdb.models['standard'].parts['Frame']
        f, e, d = p.faces, p.edges, p.datums
        t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, 
            origin=(0, 0.0, 0.0))
        s = mdb.models['standard'].ConstrainedSketch(name='__profile__', 
            sheetSize=4.16, gridSpacing=0.1, transform=t)
        g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
        s.setPrimaryObject(option=SUPERIMPOSE)
        p = mdb.models['standard'].parts['Frame']
        p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
        for i in range(0, len(outerv[0])):
            s.Line(point1=(innerv[0][i],innerv[1][i]), point2=(outerv[0][i],outerv[1][i]))
        p = mdb.models['standard'].parts['Frame']
        f = p.faces
        pickedFaces = f.getByBoundingBox(3*left, 3*bottom,0,-3*left,-3*bottom,0)
        e1, d2 = p.edges, p.datums
        p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
        s.unsetPrimaryObject()
        del mdb.models['standard'].sketches['__profile__']
    else:
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

def const_mass(radius, thickness, sides):
    """Input: radius = (f) outer radius of polygon
              thickness = (f) thickness of polygon
              sides = (int) number of polygon sides
       Function: Calculates for even sided polygons with sides 4 to 10 the
                 thickness value to maintain equal area throughout, and
                 prints the values to console.
                 Note: Uniform area = Uniform mass for the 
		 plane-strain simulation.
       Output: none
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

    print "Auto-generated list of parameters for constant mass" 

    for n in range(4,14,2):
	p_sin = math.sin(math.pi/n)
	p_cos = math.cos(math.pi/n)
        radicand = area/n/p_sin/p_cos+radius**2
        if  radicand >0:
            thick = radicand**0.5-radius
        else:
            thick = -1 
        print "%d-sides: Thickness: %f" % (n, thick) 


def getdat(session, odb_name):
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
    session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', 
        NODAL, ((COMPONENT, 'RF2'), )), ), nodeSets=('PLANE-BOT', ))
    array = session.xyDataObjects.values()
    array = tuple(array)
    xy121 = sum(array)
    xy121.setValues(
        sourceDescription="same")
    tmpName = xy121.name
    session.xyDataObjects.changeKey(tmpName, 'XYData-1')
    x0 = session.xyDataObjects['XYData-1']
    savetime = str(time.time())
    session.xyReportOptions.setValues(numberFormat=AUTOMATIC)
    session.writeXYReport(fileName='rf'+savetime+'.dat', xyData=(x0, ))
    print('rf'+savetime+'.dat')
