


def sketch(s, array):
#Given 2D array of xy coordinates, sketches polygon using shapes 
	for i in range(0, len(array[0])-1):
   		s.Line(point1=(array[0][i], array[1][i]), point2=(array[0][i+1], array[1][i+1]))
def vertices(startx, starty, radius, sides):
#Returns the coordinates of the vertices of a (sides)-dimensional polygon 
#centered at startx, starty of a certain radius

   xarray = []
   yarray = []
   for i in range(0, sides):
	angle = 2*3.14159265*i/sides
	sincos = trig(angle) 
	xval = startx + radius*sincos[0]
	yval = starty + radius*sincos[1]

	xarray.append(xval)
	yarray.append(yval)
   #append starting point so line can return to starting vertice 
   xarray.append(xarray[0])
   yarray.append(yarray[0])
   return [xarray, yarray]
def trig(x):
#Returns sin and cos of angle x (rad) 
#http://lab.polygonal.de/2007/07/18/fast-and-accurate-sinecosine-approximation/
#Low degree of precision required because CAE should automatically get correct angle once sides 
#are constrained 
	if x < -3.14159265: 
    	   x += 6.28318531
	else:
	   if (x >  3.14159265):
    	      x -= 6.28318531
#compute sine
	if x<0:
		sin = 1.27323954 * x + .405284735 * x * x
		if sin<0: 
    	   		sin = .225*(sin*(-sin)-sin)+sin
		else:
    	   		sin = .225*(sin*sin-sin)+sin
	else:
		sin = 1.27323954 * x - .405284735 * x * x
		if sin < 0:
			sin = .225*(sin*(-sin)-sin)+sin
		else:
			sin = .225*(sin*sin-sin)+sin
#compute cosine: sin(x + PI/2) = cos(x)
	x += 1.57079632;
	if x >  3.14159265:
    	  	x -= 6.28318531;
	if x<0:
    		cos = 1.27323954 * x + 0.405284735 * x * x
		if cos < 0:
			cos = .225*(cos*(-cos)-cos)+cos
		else:
			cos = .225*(cos*cos-cos) + cos
	else:
      		cos = 1.27323954 * x - 0.405284735 * x * x
		if cos < 0:
			cos = .225*(cos*(-cos)-cos)+cos
		else:
			cos = .225*(cos*cos-cos)+cos

	return sin,cos #return tuple

def radiusgen(perimeter, sides):
#Assume for fair comparison that perimeter is equal
#Calculates the radius of the polygon given the perimeter and number of sides
	sincos = trig(3.14159265/sides)
	radius = (perimeter/sides)/(2*sincos[0])
	return radius

def bc_bot(mdb, vertices):
#a = mdb.models['standard'].rootAssembly
#Given array of vertices, applies encastre BC to the bottom most point or points 
#Max is always at vertices[0][0] but only apply velocity if specified and force if specified
	output = []
	num_vert = len(vertices[0])
	a = mdb.models['standard'].rootAssembly
	v = a.instances['Frame-1'].vertices
	if num_vert%2 != 1:	#Odd number of sides so two points on the bottom need to be fixed
		region=(v.findAt(((vertices[0][num_vert/2], vertices[1][num_vert/2], 0.0), ), ), None, None, None)
		mdb.models['standard'].EncastreBC(name='Fixedright', createStepName='Initial', region=region)
		
		region=(v.findAt(((vertices[0][num_vert/2+1], vertices[1][num_vert/2+1], 0.0), ), ), None, None, None)
		mdb.models['standard'].EncastreBC(name='Fixedleft', createStepName='Initial', region=region)
	else:
		region=(v.findAt(((vertices[0][num_vert/2], vertices[1][num_vert/2], 0.0), ), ), None, None, None)
		mdb.models['standard'].EncastreBC(name='Fixedbottom', createStepName='Initial', region=region)
		
def loader(mdb, vertices, force=0, velocity=False, velx=0, vely=0, velr3=0, time=2.0, maxinc=10000, initinc=0.01, minimum=2e-05, maximum=0.025):
#Must give at least one argument for this to work 
	from abaqus import *
	from abaqusConstants import *
	a = mdb.models['standard'].rootAssembly
	v = a.instances['Frame-1'].vertices

	if velocity:
		#Create a step with increments 
		mdb.models['standard'].StaticStep(name='Apply load', previous='Initial', description='Description', timePeriod=time, maxNumInc=maxinc, 
        	initialInc=initinc, minInc=minimum, maxInc=maximum)
    		session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Apply load')

		verts1 = v.findAt(((vertices[0][0], vertices[1][0], 0.0), ))
    		region = a.Set(vertices=verts1, name='topnode')
    		mdb.models['standard'].VelocityBC(name='velocity', createStepName='Apply load', 
        	region=region, v1=velx, v2=vely, vr3=velr3, amplitude=UNSET, 
        	localCsys=None, distributionType=UNIFORM, fieldName='')

	if force != 0: 
		mdb.models['standard'].StaticLinearPerturbationStep(name='Apply load', 
    		previous='Initial', description='10kN central load', 
    		matrixSolver=SOLVER_DEFAULT)
		session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Apply load')
		mdb.models['standard'].fieldOutputRequests['F-Output-1'].setValues(
    		variables=PRESELECT, region=MODEL)
		
		region=((v.findAt(((vertices[0][0], vertices[1][0], 0.0), ), ), ), )
		mdb.models['standard'].ConcentratedForce(name='Force', 
    		createStepName='Apply load', 
    		region=region, cf2=-force)		
