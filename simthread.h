//==============================================================================
// The solver on its own thread.
//
// Until now the simulation was stepped from inside the GLUT display callback:
// one call to Fluid::progress() per drawn frame. That callback runs after
// glutSwapBuffers(), which blocks until the next vertical retrace, so the
// simulation could never take more steps per second than the monitor took
// frames - 60, on any ordinary screen, however many cores were idle.
//
// Here the solver runs in a thread of its own and steps as fast as it can. It
// publishes a snapshot - the particle positions, and the marching cubes
// surface when the viewer is drawing one - whenever the renderer asks for the
// next one, and only then, so that building the mesh costs one build per drawn
// frame rather than one per step. Three snapshots rotate between the two
// threads, so neither ever waits for the other: the solver fills one while the
// renderer draws another.
//
// Everything that touches the Fluid happens on the solver thread. The renderer
// only ever reads a published snapshot, which nothing writes to while it holds
// it.
//==============================================================================
#ifndef SIMTHREAD_H
#define SIMTHREAD_H

class Fluid;

//One consistent view of the fluid, taken between two solver steps. The arrays
//belong to the snapshot and are reused as it comes round again; they grow the
//same way the viewer's arrays used to.
struct FluidSnapshot
{
	//Particle positions, three floats each, boundary particles included - the
	//same layout Fluid::get_particlearray() has always produced.
	float         * particles;
	unsigned long   particles_alloc;
	unsigned long   particlecount;        //all of them, which is what the dots draw
	unsigned long   movingparticlecount;  //the liquid ones, which is what the water is

	//The marching cubes surface, empty unless this snapshot was asked for one.
	float         * vertices;
	float         * normals;
	unsigned long   vn_alloc;             //floats, shared by both arrays
	unsigned long   vn_count;             //vertices
	unsigned int  * indices;
	unsigned long   indices_alloc;
	unsigned long   trianglecount;
	bool            has_surface;

	unsigned long long step;              //solver steps taken when this was built
};

//Start stepping f. Nothing else may touch f afterwards until sim_stop().
//omp_threads is what -t asked for, 0 for the machine's default: it has to be
//set from inside the solver thread, because the OpenMP thread count is a
//property of the thread that asks for a parallel region, not of the program.
void sim_start ( Fluid & f, int omp_threads );
//Stop the thread and wait for it. Safe to call twice, and called at exit.
void sim_stop();

//Simulation controls. Read and written from either thread.
void  sim_set_paused ( bool paused );
bool  sim_paused();
void  sim_set_dt ( float dt );
float sim_dt();
//Steps per second the solver may take, 0 for as many as it can.
void  sim_set_rate_limit ( float steps_per_second );

//Hand back the newest published snapshot and ask for the next one, which will
//carry the surface if with_surface. Null until the solver has published its
//first. The pointer stays valid until the next call.
const FluidSnapshot * sim_acquire_snapshot ( bool with_surface );

//Steps taken since this was last called: the counterpart of the frame counter.
unsigned long long sim_steps_since_last_call();

#endif
