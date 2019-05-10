# Ray Tracing

The following code implements a raytracer in C++ which is designed to render a static scene consisting basic primitives such as Quadric (like sphere, cylinder etc), Box, Planes. Multiple point lights can be set up by the discretion of user to give the desired lighting to the scene.

## Scene Creation

All the objects present in the scene are taken as input from the user in prescribed json input format (please see any of the input files already provided - "objectInput.json"). These contain all the relevant parameters such as refractive index, coefficient of diffusion, affine transformation matrix and texture.

Also, the relevant information for worldview is also taken as input in another json file ("viewerInput.json"). The angle in the file specifies the relative size of screen, smaller it will be, smaller the size of the screen. The width - height defines the number of pixels along lenght and breadth respectively. Eyepoint is fixed at distance of 1 from the screen. Please note : by varying the angle, we are considering the effect of displacing the eye point for a fixed size screen. Any other transformations like changing the viewing direction itself is considered via the transformation matrix.

Input for lightsources is taken from json file named "lightInput.json". This file notes the location and intensity of each light source of the scene.

## Scene Rendering

Scene is rendered using the technique of backward raytracing. Phong local Illumination model along with secondary illumination i.e. reflection and refraction are considered at each point of intersaction of viewing ray and the objects present in the scene. Secondary rays were only considered upto to user defined depth of recursion.
(Please note that simulation time increases manyfolds when the texture is applied to the objects, hence better avoid using it for just simulation purposes)

## OpenGL Visualisation of the Process of Ray Tracing

The whole process of ray tracing is being rendered using opengl. All of the rays - viewing ray, shadow rays, reflection rays and refraction rays are shown with different color coding (Red, white, green, blue respectively). User mentions the pixel coordinates of viewing screen for which he wants to see the ray tracing process.

## Note

All parts except the OpenGL rendering is highly generalized. Since the opengl rendering requires the objects in mesh form, while our ray tracer only required parameters, our code can't render all the scenes of ray tracer in opengl rendering without explicit changes in the code of opengl renderer. For handling this issue, we further propose to use meshes in the raytracer too so that both can be easily linked. For demonstration purposes, one input file is already provided for this simulation purposes.  

To run this code, just clone the repo or download the zip file. Build the solution in x64/Release mode and then run the simulate.sh. File paths for various input files are not currently provided as cmd arguments, instead you may need to change the path in code itself.