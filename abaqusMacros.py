def tplanar():
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
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=4.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.Line(point1=(0.0, 0.0), point2=(0.0, 1.0))
    s.VerticalConstraint(entity=g.findAt((0.0, 0.5)), addUndoState=False)
    s.Line(point1=(0.0, 1.0), point2=(1.0, 1.0))
    s.HorizontalConstraint(entity=g.findAt((0.5, 1.0)), addUndoState=False)
    s.PerpendicularConstraint(entity1=g.findAt((0.0, 0.5)), entity2=g.findAt((0.5, 
        1.0)), addUndoState=False)
    s.Line(point1=(1.0, 1.0), point2=(1.0, 0.0))
    s.VerticalConstraint(entity=g.findAt((1.0, 0.5)), addUndoState=False)
    s.PerpendicularConstraint(entity1=g.findAt((0.5, 1.0)), entity2=g.findAt((1.0, 
        0.5)), addUndoState=False)
    s.Line(point1=(1.0, 0.0), point2=(0.0, 0.0))
    s.HorizontalConstraint(entity=g.findAt((0.5, 0.0)), addUndoState=False)
    s.PerpendicularConstraint(entity1=g.findAt((1.0, 0.5)), entity2=g.findAt((0.5, 
        0.0)), addUndoState=False)
    s.Line(point1=(0.25, 0.25), point2=(0.75, 0.25))
    s.HorizontalConstraint(entity=g.findAt((0.5, 0.25)), addUndoState=False)
    s.Line(point1=(0.75, 0.25), point2=(0.75, 0.75))
    s.VerticalConstraint(entity=g.findAt((0.75, 0.5)), addUndoState=False)
    s.PerpendicularConstraint(entity1=g.findAt((0.5, 0.25)), entity2=g.findAt((
        0.75, 0.5)), addUndoState=False)
    s.Line(point1=(0.75, 0.75), point2=(0.25, 0.75))
    s.HorizontalConstraint(entity=g.findAt((0.5, 0.75)), addUndoState=False)
    s.PerpendicularConstraint(entity1=g.findAt((0.75, 0.5)), entity2=g.findAt((0.5, 
        0.75)), addUndoState=False)
    s.Line(point1=(0.25, 0.75), point2=(0.25, 0.25))
    s.VerticalConstraint(entity=g.findAt((0.25, 0.5)), addUndoState=False)
    s.PerpendicularConstraint(entity1=g.findAt((0.5, 0.75)), entity2=g.findAt((
        0.25, 0.5)), addUndoState=False)
    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=TWO_D_PLANAR, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-1']
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']


def rigidplane():
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
    a = mdb.models['standard'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, connectors=ON, optimizationTasks=OFF, 
        geometricRestrictions=OFF, stopConditions=OFF)
    a = mdb.models['standard'].rootAssembly
    e1 = a.instances['plane-1'].edges
    edges1 = e1.findAt(((-0.619292, -0.619292, 0.0), ))
    region = a.Set(edges=edges1, name='planar')
    mdb.models['standard'].DisplacementBC(name='rigidplane', 
        createStepName='Apply load', region=region, u1=0.0, u2=0.0, ur3=0.0, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)


def orientation():
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
    p = mdb.models['standard'].parts['plane']
    e = p.edges
    edges = e.findAt(((-0.619292, -0.619292, 0.0), ))
    region = p.Set(edges=edges, name='planar')
    p = mdb.models['standard'].parts['plane']
    p.SectionAssignment(region=region, sectionName='plane', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)


def section():
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
    p1 = mdb.models['standard'].parts['plane']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    mdb.models['standard'].BeamSection(name='Section-3', 
        integration=DURING_ANALYSIS, poissonRatio=0.0, profile='Profile-1', 
        material='Rigidbody', temperatureVar=LINEAR, 
        consistentMassMatrix=False)
    mdb.models['standard'].BeamSection(name='Section-4', 
        integration=DURING_ANALYSIS, poissonRatio=0.3, profile='Profile-1', 
        material='Rigidbody', temperatureVar=LINEAR, 
        consistentMassMatrix=False)
    mdb.models['standard'].HomogeneousSolidSection(name='Section-5', 
        material='Rigidbody', thickness=None)
    mdb.models['standard'].SurfaceSection(name='Section-6', useDensity=OFF)
    p = mdb.models['standard'].parts['plane']
    e = p.edges
    edges = e.findAt(((-0.619292, -0.619292, 0.0), ))
    region=p.Set(edges=edges, name='planar')
    mdb.models['standard'].parts['plane'].sectionAssignments[0].setValues(
        region=region)


def property():
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


