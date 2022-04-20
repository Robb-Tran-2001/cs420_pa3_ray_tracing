Assignment #3: Ray tracing

FULL NAME: Robb Tran

Developed on Windows 10 and VS2019, run in release mode.

MANDATORY FEATURES
------------------

<Under "Status" please indicate whether it has been implemented and is
functioning correctly.  If not, please explain the current status.>

Feature:                                 Status: finish? (yes/no)
-------------------------------------    -------------------------
1) Ray tracing triangles                 yes

2) Ray tracing sphere                    yes

3) Triangle Phong Shading                yes

4) Sphere Phong Shading                  yes

5) Shadows rays                          yes

6) Still images                          yes
   
7) Extra Credit (up to 20 points)
	Recursive Reflection: when doing radiance calculation, I use recursive reflection to
determine the final light and color after the initial ray is reflected by as many times
as the macro MAX_REFLECTION.
	Antialiasing: instead of initializing one ray at a time, I initialize 4 rays in a grid
around the initial ray in the North, South, East, West direction by a little distance.
The pixel color is the average of those 4.
	Soft shadow: instead of one dark direct shadow, I create soft shadows by reducing
the initial light color intensity by the SUB_LIGHT macro, as well as creating the same 
number of lights for each initial light source, each at a random position.

Other:
	MACRO values to change for different effects.
	Moller-Trumbore intersection: ray-triangle intersection, fast method to calculate
ray-triangle intersection without precomputing plane equation the triangle is in.

Files to show off features all with antialiasing and softshadow:
- test1.jpg (30 sublights)
- test2.jpg (30 sublights)
- spheres.jpg (30 sublights)
- snow.jpg (10 sublights, took too long with 30)
- table.jpg (30 sublights)
- toy.jpg (10 sublights, took too long with 30)
- SIGGRAPH.jpg (5 sublights, took too long with 30)