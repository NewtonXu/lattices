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
import polygonggmodule
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
Mdb()
##
##Define inputs here
##
startx = 0
starty = 0
radius = 0.67
sides = 8
thickness = 0.1
print("Working as intended")

vex = polygonggmodule.vertices(startx, starty, radius, sides)
left = min(vex[0])
bottom = min(vex[1])
mdb.models.changeKey(fromName='Model-1', toName='standard')

##
##  Sketch polygon
##
s = mdb.models['standard'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
polygonggmodule.sketch(s,vex)
vex2 = polygonggmodule.vertices(startx, starty, radius-thickness, sides)
polygonggmodule.sketch(s,vex2) 
p = mdb.models['standard'].Part(name='Frame', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['standard'].parts['Frame']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['standard'].parts['Frame']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['standard'].sketches['__profile__']
y

##
##  Sketch plane
##
s = mdb.models['standard'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
polygonggmodule.line(s,left,-bottom*1.01)
p = mdb.models['standard'].Part(name='plane', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['standard'].parts['plane']
p.BaseWire(sketch=s)
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
polygonggmodule.line(s,left,bottom*1.01)
p = mdb.models['standard'].Part(name='planetop', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['standard'].parts['planetop']
p.BaseWire(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['standard'].parts['planetop']
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
edges = e.getByBoundingBox(-3.1*radius, -3*radius,0,3.1*radius,0,0)
region = p.Set(edges=edges, name='planar')
p = mdb.models['standard'].parts['plane']
p.SectionAssignment(region=region, sectionName='plane', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

p1 = mdb.models['standard'].parts['planetop']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['standard'].parts['planetop']
e = p.edges
edges = e.getByBoundingBox(-3.1*radius, 0,0,3.1*radius,3*radius,0)
region = p.Set(edges=edges, name='planar2')
p = mdb.models['standard'].parts['planetop']
p.SectionAssignment(region=region, sectionName='plane', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)


p = mdb.models['standard'].parts['Frame']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
mdb.models['standard'].HomogeneousSolidSection(name='Section-3', 
    material='Steel', thickness=1.0)
p = mdb.models['standard'].parts['Frame']
f = p.faces
faces = f.findAt(((0.364377, 0.467393, 0.0), ))
region = p.Set(faces=faces, name='everything')
p = mdb.models['standard'].parts['Frame']
p.SectionAssignment(region=region, sectionName='Section-3', offset=0.0, 
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
region=p.sets['planar']
p = mdb.models['standard'].parts['plane']
p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0,-1.0))

p = mdb.models['standard'].parts['planetop']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['planetop']
region=p.sets['planar2']
p = mdb.models['standard'].parts['planetop']
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
p = mdb.models['standard'].parts['planetop']
a.Instance(name='planetop-1', part=p, dependent=ON)
#p1 = a.instances['Frame-1']
#p1.translate(vector=(-0.035794, 0.331227, 0.0))
##
##
##
##  Apply velocity to top plane
##
#v = a.instances['Frame-1'].vertices
polygonggmodule.loader(mdb, radius, vex, velocity = True, vely=-0.001, time = 50)
##
##  Apply encastre bc to bottom left corner
##
polygonggmodule.bc_bot(mdb, vex, thickness) 
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
#p.seedPart(size=radius/10, deviationFactor=0.1, minSizeFactor=0.1)
p.seedPart(size=radius, deviationFactor=0.1, minSizeFactor=0.5)

p = mdb.models['standard'].parts['plane']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['plane']
p.deleteMesh()
p = mdb.models['standard'].parts['plane']
p.seedPart(size=radius*2, minSizeFactor=0.99)
p = mdb.models['standard'].parts['planetop']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['planetop']
p.deleteMesh()
p = mdb.models['standard'].parts['planetop']
p.seedPart(size=radius*2, minSizeFactor=0.99)

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
p = mdb.models['standard'].parts['planetop']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['standard'].parts['planetop']
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
