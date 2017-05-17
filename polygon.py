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
session.journalOptions.setValues(replayGeometry=COORDINATE, 
    recoverGeometry=COORDINATE)
import polygonmodule
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
Mdb()
startx = 0
starty = 0
radius = 0.67
sides = 8
left = -radius
bottom = -radius
vex = polygonmodule.vertices(startx, starty, radius, sides)
mdb.models.changeKey(fromName='Model-1', toName='standard')

##
##  Sketch polygon
##
s = mdb.models['standard'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
polygonmodule.sketch(s,vex)
p = mdb.models['standard'].Part(name='Frame', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['standard'].parts['Frame']
p.BaseWire(sketch=s)
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
polygonmodule.line(s,vex)
p = mdb.models['standard'].Part(name='plane', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['standard'].parts['plane']
p.BaseWire(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['standard'].parts['plane']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['standard'].sketches['__profile__']

##
##  Create material 'Steel'
##
mdb.models['standard'].Material('Steel')
mdb.models['standard'].materials['Steel'].Elastic(table=((200000, 0.3), ))

##
##  Create material 'Rigidbody'
##
mdb.models['standard'].Material('Rigidbody')
mdb.models['standard'].materials['Rigidbody'].Elastic(table=((20000000,),))

##
##  Create beam sections
##
mdb.models['standard'].BoxProfile(name='Profile-1', b=1.0, a=1.0, 
        uniformThickness=ON, t1=0.25)
mdb.models['standard'].BeamSection(name='plane', 
        integration=DURING_ANALYSIS, poissonRatio=0.0, profile='Profile-1', 
        material='Rigidbody', temperatureVar=LINEAR, 
        consistentMassMatrix=False)
mdb.models['standard'].BeamSection(name='frame', integration=DURING_ANALYSIS, 
        poissonRatio=0.3, profile='Profile-1', material='Steel', 
        temperatureVar=LINEAR, consistentMassMatrix=False)

##
##  Assign FRAME AND PLANE section
##
p = mdb.models['standard'].parts['Frame']
e = p.edges
edges = e.getByBoundingBox(-3*radius, -3*radius, 0, 3*radius, 3*radius,0)
region = p.Set(edges=edges, name='polygon')
p = mdb.models['standard'].parts['Frame']
p.SectionAssignment(region=region, sectionName='frame', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
p1 = mdb.models['standard'].parts['plane']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['standard'].parts['plane']
e = p.edges
edges = e.getByBoundingBox(-3*radius, -1.01*radius,0,3*radius,-0.99*radius,0)
region = p.Set(edges=edges, name='plane')
p = mdb.models['standard'].parts['plane']
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
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
        engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)
p = mdb.models['standard'].parts['plane']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['plane']
region=p.sets['plane']
p = mdb.models['standard'].parts['plane']
p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0,-1.0))
p = mdb.models['standard'].parts['Frame']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['Frame']
region=p.sets['polygon']
p = mdb.models['standard'].parts['Frame']
p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0,-1.0))

##
##  Instance the frame and plane
##
a = mdb.models['standard'].rootAssembly
p = mdb.models['standard'].parts['Frame']
a.Instance(name='Frame-1', part=p, dependent=ON)
p = mdb.models['standard'].parts['plane']
a.Instance(name='plane-1', part=p, dependent=ON)
#p1 = a.instances['Frame-1']
#p1.translate(vector=(-0.035794, 0.331227, 0.0))
##
##
##
##  Apply concentrated force to top node
##
#v = a.instances['Frame-1'].vertices
polygonmodule.loader(mdb, vex, velocity = True, vely=-0.001, time = 50)
##
##  Apply encastre bc to bottom left corner
##
polygonmodule.bc_bot(mdb, vex) 
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
p.seedPart(size=radius/10, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['standard'].parts['plane']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['plane']
p.deleteMesh()
p = mdb.models['standard'].parts['plane']
p.seedPart(size=radius, minSizeFactor=0.99)
##
##  Assign element type
##
#elemType1 = mesh.ElemType(elemCode=B21, elemLibrary=STANDARD)
#p = mdb.models['standard'].parts['Frame']
#e = p.edges
#edges = e.getByBoundingBox(-3*left, 1.01*bottom,0,3*left,-1.01*bottom)
#pickedRegions =(edges, )
#p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))
#p = mdb.models['standard'].parts['plane']
#session.viewports['Viewport: 1'].setValues(displayedObject=p)
#elemType1 = mesh.ElemType(elemCode=B21, elemLibrary=STANDARD)
#p = mdb.models['standard'].parts['plane']
#e = p.edges
#edges = e.getByBoundingBox(-3*left, 1.01*bottom, 0, 3*left, 0.99*bottom)
#pickedRegions =(edges, )
#p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))
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
mdb.models['standard'].historyOutputRequests['H-Output-1'].setValues(
        frequency=100)
mdb.models['standard'].fieldOutputRequests['F-Output-1'].setValues(
        frequency=100)
##
##  Create job
##
mdb.Job(name='Frame', model='standard', 
    description='Two-dimensional overhead hoist frame')
mdb.jobs['Frame'].setValues(echoPrint=ON, modelPrint=ON, contactPrint=ON, 
    historyPrint=ON)

session.viewports['Viewport: 1'].view.fitView()

##
##  Save model database
##
mdb.saveAs('Frame')

a = mdb.models['standard'].rootAssembly
a.regenerate()

mdb.save()
