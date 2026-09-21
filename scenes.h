#ifndef SCENES_H
#define SCENES_H

#include "fluid.h"

//A scene is everything about a run that is not the solver: the terrain, the
//liquid, the control geometry that pushes it around, and where to point the
//camera. main() used to hold exactly one of these with the others commented
//out in place, so choosing between them meant an edit and a rebuild; -X or
//--scene picks one at run time.
struct Scene
{
	const char * name;
	const char * summary;

	//Terrain. Plain function pointers, because that is what
	//Fluid::set_height_function takes.
	float  ( *height ) ( float x, float y );
	Vektor ( *height_normal ) ( float x, float y );

	//Put the liquid and the control geometry in. wanted_particles is -n.
	void ( *build ) ( Fluid & f, unsigned long wanted_particles );

	//What -n means here: the particle count the scene is written for, so that
	//-n scales its liquid the way it scales the fountain's.
	unsigned long natural_particles;

	//The patch of ground the viewer draws, and at what resolution. The terrain
	//is a function and not geometry, so this is purely a drawing choice - as
	//is whether to draw it at all, which the fountain does not, because its
	//bowl climbs to z=180 within the patch and the camera starts ten units
	//from the origin, inside it. b turns it on.
	bool  ground_shown;
	float ground_x0, ground_x1;
	int   ground_nx;
	float ground_y0, ground_y1;
	int   ground_ny;

	//Whether the viewer should draw this scene's boundary funnel. The funnel is
	//the only piece of control geometry in any scene with a shape worth a
	//wireframe of its own; everything else reads well enough as the dim dots
	//its control particles are drawn as.
	bool draw_funnel;

	//Where to put the camera before the first frame: how far back, and how far
	//to tip the scene towards the viewer. Straight down is fine for a fountain
	//and useless for a valley.
	float camera_distance;
	float camera_pitch;

	//How much of its speed a particle keeps across a bounce off the ground.
	//0.8 is what the solver has always used.
	float ground_restitution;
};

//0 if there is no scene by that name.
const Scene * scene_by_name ( const char * name );
//Every scene, one line each, on stdout.
void scene_list();
//What --scene defaults to.
const Scene * scene_default();
//One ring vertex of the funnel scene, so that the viewer can draw the same
//funnel the solver was given. j is the ring, 1..5; i is the step, 0..10.
void funnel_ring ( int j, int i, float * out );
//Terrain scene tunables, applied before build() runs.
void terrain_set ( float slope, float waviness );

#endif /*SCENES_H*/
