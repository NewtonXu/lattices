# lattices

Scripts intended for use with ABAQUS Finite Element Analysis to generate n-sided polygons for simulation.

polygonmodule.py is a python module with functions that abstract the generation of polygons and the application of load to relevant nodes.

polygon.py is a script that goes through the entire ABAQUS workflow from shape generation to creation of job. It is used to generate input files for the polygon simulations using the polygonmodule. 
