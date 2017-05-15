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
vex = polygonmodule.vertices(0,0,0.67,6)
mdb.models.changeKey(fromName='Model-1', toName='standard')

##
##  Sketch profile of frame
##
s = mdb.models['standard'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
polygonmodule.sketch(s, vex)
p = mdb.models['standard'].Part(name='Frame', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['standard'].parts['Frame']
p.BaseWire(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['standard'].parts['Frame']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['standard'].sketches['__profile__']
##
##  Create material 'Steel'
##
mdb.models['standard'].Material('Steel')
mdb.models['standard'].materials['Steel'].Elastic(table=((200000, 0.3), ))
##
##  Create truss section
##
mdb.models['standard'].TrussSection(name='FrameSection', material='Steel', 
    area=1.963E-05)
##
##  Assign truss section
##
e = p.edges
edges = e
region = regionToolset.Region(edges=edges)
p.SectionAssignment(region=region, sectionName='FrameSection')

a = mdb.models['standard'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
##
##  Set coordinate system (done by default)
##
a = mdb.models['standard'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
##
##  Instance the frame
##
p = mdb.models['standard'].parts['Frame']
a.Instance(name='Frame-1', part=p, dependent=ON)
#p1 = a.instances['Frame-1']
#p1.translate(vector=(-0.035794, 0.331227, 0.0))
##
##  Create a static linear perturbation step
##
##
##  Apply concentrated force to bottom center
##
#v = a.instances['Frame-1'].vertices
polygonmodule.loader(mdb, vex, velocity = True, vely = -0.001, time = 50)
##
##  Apply encastre bc to bottom left corner
##
polygonmodule.bc_bot(mdb, vex) 
##
##  Assign global seed
##
p.seedPart(size=0.1)
##
##  Assign element type
##
elemType1 = mesh.ElemType(elemCode=T2D2)
e = p.edges
edges = e
pickedRegions =(edges, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))
##
##  Generate mesh
##
p.generateMesh()
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
