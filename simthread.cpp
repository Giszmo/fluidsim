#include "simthread.h"
#include "fluid.h"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <mutex>
#include <thread>

namespace
{

Fluid * g_fluid = 0;

std::thread             g_thread;
std::mutex              g_mutex;      //guards the three pointers below
std::condition_variable g_wake;       //solver waits on it when it has nothing to do

//The three snapshots. The solver owns g_build, the renderer owns g_shown, and
//g_ready is the hand-over: both swap with it under g_mutex, neither ever waits
//for the other to finish with a buffer.
FluidSnapshot * g_build = 0;
FluidSnapshot * g_ready = 0;
FluidSnapshot * g_shown = 0;
bool            g_ready_fresh = false;

std::atomic<bool>  g_quit ( false );
std::atomic<bool>  g_paused ( true );
std::atomic<bool>  g_want_snapshot ( false );
std::atomic<bool>  g_want_surface ( false );
std::atomic<float> g_dt ( 0.004f );
std::atomic<float> g_rate_limit ( 0.0f );
std::atomic<unsigned long long> g_steps ( 0 );
int g_omp_threads = 0;
unsigned long long g_steps_reported = 0;

FluidSnapshot * new_snapshot()
{
	FluidSnapshot * s = new FluidSnapshot;
	//The arrays start empty; Fluid's resize helpers grow them on first use, the
	//same way the viewer's own arrays used to grow.
	s->particles = new float[3];
	s->particles_alloc = 0;
	s->particlecount = 0;
	s->movingparticlecount = 0;
	s->vertices = new float[3];
	s->normals  = new float[3];
	s->vn_alloc = 0;
	s->vn_count = 0;
	s->indices = new unsigned int[3];
	s->indices_alloc = 0;
	s->trianglecount = 0;
	s->has_surface = false;
	s->step = 0;
	return s;
}

void fill_snapshot ( FluidSnapshot * s, bool with_surface )
{
	Fluid & f = *g_fluid;
	//The surface first: it reads the cell lists the last step left behind, and
	//get_particlearray() does not disturb them either way.
	if ( with_surface )
	{
		s->trianglecount = f.get_surfacegrid ( &s->vertices, &s->normals,
		                                       s->vn_alloc, s->vn_count,
		                                       &s->indices, s->indices_alloc );
		s->has_surface = true;
	}
	else
	{
		s->trianglecount = 0;
		s->vn_count = 0;
		s->has_surface = false;
	}
	f.get_particlearray ( &s->particles, s->particles_alloc );
	s->particlecount = f.particlecount();
	s->movingparticlecount = f.movingparticlecount();
	s->step = g_steps.load();
}

void publish ( )
{
	std::lock_guard<std::mutex> lock ( g_mutex );
	std::swap ( g_build, g_ready );
	g_ready_fresh = true;
}

void solver_loop()
{
	//Has to happen here and not in main(): omp_set_num_threads() sets the count
	//for the calling thread, and the thread that calls the parallel regions is
	//this one.
#ifdef _OPENMP
	if ( g_omp_threads > 0 )
		omp_set_num_threads ( g_omp_threads );
#endif

	//When a step rate is asked for, this is when the next step is due. Sleeping
	//to a deadline rather than for a fixed interval keeps the average right even
	//though a step never takes exactly as long as the one before it.
	std::chrono::steady_clock::time_point next_step = std::chrono::steady_clock::now();

	while ( !g_quit.load() )
	{
		//Paused, and nobody is waiting for a picture of it: wait rather than spin.
		if ( g_paused.load() && !g_want_snapshot.load() )
		{
			std::unique_lock<std::mutex> lock ( g_mutex );
			g_wake.wait_for ( lock, std::chrono::milliseconds ( 2 ) );
			continue;
		}

		if ( !g_paused.load() )
		{
			const float limit = g_rate_limit.load();
			if ( limit > 0.0f )
			{
				std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
				if ( now < next_step )
					std::this_thread::sleep_for ( next_step - now );
				else if ( next_step + std::chrono::seconds ( 1 ) < now )
					next_step = now;  //fell far behind; do not try to catch up
				next_step += std::chrono::nanoseconds (
				                 ( long long ) ( 1.0e9 / limit ) );
			}
			g_fluid->progress ( g_dt.load() );
			g_steps.fetch_add ( 1 );
		}

		//Only build a snapshot when one was asked for, so that the mesh is built
		//once per drawn frame and not once per step. A request that arrives while
		//we are building is served on the next turn round the loop.
		if ( g_want_snapshot.exchange ( false ) )
		{
			fill_snapshot ( g_build, g_want_surface.load() );
			publish();
		}
	}
}

} //namespace

void sim_start ( Fluid & f, int omp_threads )
{
	if ( g_thread.joinable() )
		return;
	g_fluid = &f;
	g_omp_threads = omp_threads;
	g_build = new_snapshot();
	g_ready = new_snapshot();
	g_shown = new_snapshot();
	g_quit.store ( false );
	//Leaving the thread running into static destruction is the one way this can
	//crash on the way out, and the viewer exits from three different places.
	atexit ( sim_stop );
	g_thread = std::thread ( solver_loop );
}

void sim_stop()
{
	if ( !g_thread.joinable() )
		return;
	g_quit.store ( true );
	g_wake.notify_all();
	g_thread.join();
}

void sim_set_paused ( bool paused )
{
	g_paused.store ( paused );
	g_wake.notify_all();
}

bool sim_paused()
{
	return g_paused.load();
}

void sim_set_dt ( float dt )
{
	g_dt.store ( dt );
}

float sim_dt()
{
	return g_dt.load();
}

void sim_set_rate_limit ( float steps_per_second )
{
	g_rate_limit.store ( steps_per_second );
}

const FluidSnapshot * sim_acquire_snapshot ( bool with_surface )
{
	g_want_surface.store ( with_surface );
	{
		std::lock_guard<std::mutex> lock ( g_mutex );
		if ( g_ready_fresh )
		{
			std::swap ( g_shown, g_ready );
			g_ready_fresh = false;
		}
	}
	g_want_snapshot.store ( true );
	g_wake.notify_all();
	//step is 0 only before the solver has ever published: an unfilled snapshot
	//has no particles in it and must not be drawn.
	return g_shown->particlecount ? g_shown : 0;
}

unsigned long long sim_steps_since_last_call()
{
	unsigned long long now = g_steps.load();
	unsigned long long d = now - g_steps_reported;
	g_steps_reported = now;
	return d;
}
