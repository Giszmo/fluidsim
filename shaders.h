//==============================================================================
// Shaders for the fluid surface.
//
// TRUFORM, which this program used in 2005, was ATI's fixed-function PN
// triangles: hand the card a triangle with a normal at each corner and it
// builds a cubic Bezier patch for the positions, a quadratic one for the
// normals, and subdivides. The extension is gone and no current driver has it.
// Its descendant is the tessellation stage of OpenGL 4.0, which does the same
// thing except that the patch is evaluated by shaders you write yourself -
// surface_pn_* below is that shader, and it is the same maths Vlachos et al.
// published in 2001.
//
// Two things it can do that TRUFORM could not: the subdivision level is chosen
// per edge, so it is set here from how long the edge is on screen, and a patch
// far away costs nothing. Crack-free because the level and the curve of a
// shared edge depend only on that edge's two endpoints, so both triangles
// sitting on it compute the same numbers.
//
// Everything is a compatibility-profile shader on purpose: it reads
// gl_ModelViewProjectionMatrix, gl_LightSource[] and gl_FrontMaterial, so the
// camera, the two lights and the material stay exactly where the rest of this
// program sets them and the tessellated surface matches the untessellated one.
//==============================================================================
#ifndef SHADERS_H
#define SHADERS_H

//Compile the surface programs. Call once for each GL context, after the window
//exists. Returns false if this machine has no OpenGL 4 tessellation, in which
//case the caller keeps drawing plain triangles.
bool surface_shaders_init();

//False when there is no tessellation on this machine, or the shaders did not
//compile. surface_shaders_error() then says why, for the status line.
bool surface_tessellation_available();
const char * surface_shaders_error();

//Bind the PN triangle program. pixels_per_segment is the edge length in pixels
//that one tessellated segment aims for - smaller is smoother and slower.
//max_level caps the subdivision. Draw with GL_PATCHES and 3 vertices a patch.
void surface_tessellation_bind ( float pixels_per_segment, float max_level );
void surface_tessellation_unbind();

//------------------------------------------------------------------------------
// Screen-space fluid, the other way to draw the same particles
//------------------------------------------------------------------------------
//
// No mesh at all. The particles are splatted as spheres into a depth buffer,
// the depth is smoothed, and the surface is shaded from the normals of the
// smoothed depth - van der Laan, Green and Sainz, "Screen Space Fluid Rendering
// with Curvature Flow", I3D 2009. It is what the real-time fluid renderers do,
// because there is no marching cubes pass, nothing to re-upload per frame, and
// the cost is per pixel rather than per particle.
//
// What it cannot give you is geometry: no reflections off the water, no shadow
// casting, nothing to export. That is what the tessellated marching-cubes path
// above is for.

//Compiled by surface_shaders_init() along with the rest.
bool screenspace_available();
const char * screenspace_error();

//Draw count particles, three floats each, as a fluid surface. radius is in
//world units and smoothing is how many times the depth is blurred - more is
//smoother and slower. Uses the projection and lights already set, writes depth,
//and leaves the GL state as it found it. False if it could not draw.
bool screenspace_render ( const float * positions, unsigned long count,
                          float radius, int smoothing );

//Give the render targets back. Safe on a context that is already gone.
void screenspace_free();

#endif
