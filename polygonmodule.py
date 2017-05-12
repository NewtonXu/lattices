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
	xarray.append(startx + radius*sincos[0])
	yarray.append(starty + radius*sincos[1])
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
    	   sin = 1.27323954 * x + .405284735 * x * x;
	else:
    	   sin = 1.27323954 * x - 0.405284735 * x * x;
#compute cosine: sin(x + PI/2) = cos(x)
	x += 1.57079632;
	if x >  3.14159265:
    	   x -= 6.28318531;
	if x<0:
    	   cos = 1.27323954 * x + 0.405284735 * x * x
	else:
    	   cos = 1.27323954 * x - 0.405284735 * x * x;
	return sin,cos #return tuple

def radiusgen(perimeter, sides):
#Assume for fair comparison that perimeter is equal
#Calculates the radius of the polygon given the perimeter and number of sides
	sincos = trig(3.14159265/sides)
	radius = (perimeter/sides)/(2*sincos[0])
	return radius

