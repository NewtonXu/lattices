#
# Getting Started with Abaqus: Interactive Edition
#
# Script for frame example
#
#

from abaqus import *
from abaqusConstants import *
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
import math
import polygonmodule
import os
os.chdir(r"/mnt/compute-0-4/newton")
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
Mdb()
##
##Define inputs here
##
polygontype = 1
startx = 0
starty = 0
radius = 0.5
sides = 4
thickness = 0.071967

if polygontype == 0:
	planethick = thickness
else:
	planethick = thickness * math.sin(math.pi/2-math.pi/sides)
distance = 0.085
sidelen = polygonmodule.sidelength(radius, sides)
apothem = sidelen/2/math.tan(math.pi/sides)
seeder = 0.01
print("Working as intended")
#vex = polygonmodule.vertices(startx, starty, radius, sides)
vex = polygonmodule.advancedpolygon(radius, thickness, sides, distance)


left = min(vex[0][0])
bottom = min(vex[0][1])
mdb.models.changeKey(fromName='Model-1', toName='standard')

##
##  Sketch polygon	
##
s = mdb.models['standard'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
polygonmodule.sketch(s,vex[0])
polygonmodule.sketch(s,vex[1])	
p = mdb.models['standard'].Part(name='Frame', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['standard'].parts['Frame']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['standard'].parts['Frame']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['standard'].sketches['__profile__']

##
##  Sketch plane
##
s = mdb.models['standard'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
polygonmodule.rectangle(s,left,bottom+planethick, -planethick)
p = mdb.models['standard'].Part(name='plane', dimensionality=TWO_D_PLANAR, 
   		type=DEFORMABLE_BODY)
p = mdb.models['standard'].parts['plane']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['standard'].parts['plane']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['standard'].sketches['__profile__']
##
##  Sketch top plane
##

s = mdb.models['standard'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
polygonmodule.rectangle(s,left,-bottom-planethick, planethick)
p = mdb.models['standard'].Part(name='planetop', dimensionality=TWO_D_PLANAR, 
   		type=DEFORMABLE_BODY)
p = mdb.models['standard'].parts['planetop']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['standard'].parts['planetop']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['standard'].sketches['__profile__']

##
##  Instance the frame and plane
##
a = mdb.models['standard'].rootAssembly
p = mdb.models['standard'].parts['Frame']
a.Instance(name='Frame-1', part=p, dependent=ON)
p = mdb.models['standard'].parts['plane']
a.Instance(name='plane-1', part=p, dependent=ON)
p = mdb.models['standard'].parts['planetop']
a.Instance(name='planetop-1', part=p, dependent=ON)
if polygontype == 1:
	a.InstanceFromBooleanMerge(name='Part-3', instances=(a.instances['Frame-1'], 
        a.instances['plane-1'], ), originalInstances=DELETE, domain=GEOMETRY)
    	a = mdb.models['standard'].rootAssembly
    	a.InstanceFromBooleanMerge(name='Frame', instances=(
        a.instances['planetop-1'], a.instances['Part-3-1'], ), 
        originalInstances=DELETE, domain=GEOMETRY)
##
##  Create material 'Steel'
##
mdb.models['standard'].Material('Steel')
mdb.models['standard'].materials['Steel'].Elastic(table=((206010.0, 0.3), ))
mdb.models['standard'].materials['Steel'].Plastic(table=((260.0, 
        0.0), (262.01, 0.02), (282.695, 0.04), (317.17, 0.06), (330.96, 0.08), 
        (344.75, 0.1), (344.75, 0.12), (330.96, 0.14), (324.065, 0.16), (
        303.38, 0.18), (289.59, 0.2), (275.8, 0.22), (248.22, 0.24)))
mdb.models['standard'].materials['Steel'].plastic.RateDependent(type=JOHNSON_COOK, table=((0.011, 25.0), ))
mdb.models['standard'].materials['Steel'].Density(table=((8e-09, ), ))

##
##  Create material 'Rigidbody'
##
mdb.models['standard'].Material('Rigidbody')
mdb.models['standard'].materials['Rigidbody'].Elastic(table=((20000000,),))
mdb.models['standard'].materials['Rigidbody'].Density(table=((7.85e-09, ), ))
##
##  Create beam sections
##
##mdb.models['standard'].BoxProfile(name='Profile-1', b=1.0, a=1.0, 
##        uniformThickness=ON, t1=0.25)
##mdb.models['standard'].BeamSection(name='plane', 
##        integration=DURING_ANALYSIS, poissonRatio=0.0, profile='Profile-1', 
##        material='Rigidbody', temperatureVar=LINEAR, 
##        consistentMassMatrix=False)

##
##Create homogenous solids
##
mdb.models['standard'].HomogeneousSolidSection(name='frame', 
        material='Steel', thickness=None)
mdb.models['standard'].HomogeneousSolidSection(name='plane', 
        material='Steel', thickness=None)


##
##  Assign FRAME AND PLANE section
##
p = mdb.models['standard'].parts['Frame']
e = p.faces	#e=p.edges 
edges = e.getByBoundingBox(3*left, 3*bottom, 0, -3*left, -3*bottom,0)
region = p.Set(faces=edges, name='polygon')
p = mdb.models['standard'].parts['Frame']
p.SectionAssignment(region=region, sectionName='frame', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
if polygontype == 0:
	p1 = mdb.models['standard'].parts['plane']
	#session.viewports['Viewport: 1'].setValues(displayedObject=p1)
	p = mdb.models['standard'].parts['plane']
	e = p.faces #e=p.edges
	edges = e.getByBoundingBox(3.1*left, 3*bottom,0,-3.1*left,0,0)
	region = p.Set(faces=edges, name='planar')
	p = mdb.models['standard'].parts['plane']
	p.SectionAssignment(region=region, sectionName='plane', offset=0.0, 
        	offsetType=MIDDLE_SURFACE, offsetField='', 
        	thicknessAssignment=FROM_SECTION)

	p1 = mdb.models['standard'].parts['planetop']
	#session.viewports['Viewport: 1'].setValues(displayedObject=p1)
	p = mdb.models['standard'].parts['planetop']
	e = p.faces #e=p.edges
	edges = e.getByBoundingBox(10*left, 0,0,-10*left,-3*bottom,0)
	region = p.Set(faces=edges, name='planar2')
	p = mdb.models['standard'].parts['planetop']
	p.SectionAssignment(region=region, sectionName='plane', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

##
##  Set coordinate system (done by default)
##
a = mdb.models['standard'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)

##
##  Set beam orientation
##
#session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
#        engineeringFeatures=ON)
#session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
#        referenceRepresentation=OFF)
#p = mdb.models['standard'].parts['plane']
#session.viewports['Viewport: 1'].setValues(displayedObject=p)
#p = mdb.models['standard'].parts['plane']
#region=p.sets['planar']
#p = mdb.models['standard'].parts['plane']
#p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0,-1.0))

#p = mdb.models['standard'].parts['planetop']
#session.viewports['Viewport: 1'].setValues(displayedObject=p)
#p = mdb.models['standard'].parts['planetop']
#region=p.sets['planar2']
#p = mdb.models['standard'].parts['planetop']
#p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0,-1.0))

#p = mdb.models['standard'].parts['Frame']
#session.viewports['Viewport: 1'].setValues(displayedObject=p)
#p = mdb.models['standard'].parts['Frame']
#region=p.sets['polygon']
#p = mdb.models['standard'].parts['Frame']
#p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0,-1.0))



#p1 = a.instances['Frame-1']
#p1.translate(vector=(-0.035794, 0.331227, 0.0))

##
##   Create the step 
##
mdb.models['standard'].ExplicitDynamicsStep(name='Apply load', previous='Initial', 
        timeIncrementationMethod=FIXED_USER_DEFINED_INC, userDefinedInc=0.000000005)
mdb.models['standard'].steps['Apply load'].setValues(
        timeIncrementationMethod=AUTOMATIC_GLOBAL, scaleFactor=1.0, 
        maxIncrement=None)

##
##  Apply load to top plane
##
#v = a.instances['Frame-1'].vertices
polygonmodule.loader(polygontype, mdb, left, bottom, radius, thickness, planethick, 3)

##
##  Apply bc and initiate interactions between planes and polygon
##
polygonmodule.bc_bot(polygontype, mdb, vex, planethick)
 
##
## Generate polygon to polygon contact interaction
##
if sides !=4:
	polygonmodule.edgeselector(polygontype, mdb,vex,radius, thickness, sides, apothem)
else:
	polygonmodule.boxselector(mdb, vex, radius, thickness, distance)
##
##  Assign global seed
##
p = mdb.models['standard'].parts['Frame']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
    engineeringFeatures=OFF, mesh=ON)
session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
    meshTechnique=ON)
p = mdb.models['standard'].parts['Frame']
p.seedPart(size=seeder, deviationFactor=0.1, minSizeFactor=0.5)

p = mdb.models['standard'].parts['plane']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['plane']
p.deleteMesh()
p = mdb.models['standard'].parts['plane']
p.seedPart(size=seeder/2, minSizeFactor=0.99)
p = mdb.models['standard'].parts['planetop']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['planetop']
p.deleteMesh()
p = mdb.models['standard'].parts['planetop']
p.seedPart(size=seeder/2, minSizeFactor=0.99)


##
##  Generate partition
##
arrays = polygonmodule.advancedarrayspolygon(radius, thickness,sides, distance)
polygonmodule.partition(polygontype, mdb, arrays[0],arrays[1], radius, distance, thickness)
##
##  Generate mesh
##
p = mdb.models['standard'].parts['Frame']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['Frame']
p.generateMesh()
p = mdb.models['standard'].parts['plane']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['plane']
p.generateMesh()
p = mdb.models['standard'].parts['planetop']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['planetop']
p.generateMesh()
mdb.models['standard'].historyOutputRequests['H-Output-1'].setValues(
        frequency=100)
#mdb.models['standard'].fieldOutputRequests['F-Output-1'].setValues(
#        frequency=100)
mdb.models['standard'].fieldOutputRequests['F-Output-1'].setValues(
        timeInterval=1e-07)

##
##Assign mesh type 
##
elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=EXPLICIT, 
        secondOrderAccuracy=OFF, hourglassControl=DEFAULT, 
        distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=EXPLICIT)
p = mdb.models['standard'].parts['Frame']
f = p.faces
faces = f.getByBoundingBox(3*left, 3*bottom,0, -3*left, -3*bottom, 0) 
pickedRegions =(faces, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))

##
##  Create job
##
mdb.Job(name='%d-sides' % (sides), model='standard', 
    description='Two-dimensional overhead hoist frame')
mdb.jobs['%d-sides' % (sides)].setValues(echoPrint=ON, modelPrint=ON, contactPrint=ON, 
    historyPrint=ON)

session.viewports['Viewport: 1'].view.fitView()

##
##  Save model database
##
#mdb.saveAs('Frame')

a = mdb.models['standard'].rootAssembly
a.regenerate()

#mdb.save()

