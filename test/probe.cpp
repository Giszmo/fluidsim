// What actually stops a particle in a scene.
//
// Builds a scene exactly as the viewer does, runs it without a window, and
// counts every correction the solver applies to a particle's velocity, split
// by whether the particle was moving into the surface or already moving away
// from it. The second number is the interesting one: a correction applied to
// a particle that is on its way out takes speed off it for nothing, and with
// a low restitution it is what makes water stick to a wall.
//
// It also counts the patches that move a particle rather than turn it - the
// shifts and teleports a scene closes its loop with - and says where the water
// ended up, which is how you tell a loop that is turning from one that is
// quietly emptying the scene into the distance.
//
//     make test/probe
//     ./test/probe funnel 4000
//     ./test/probe funnel 25000 -n 833          #a longer run, or more water
//     ./test/probe funnel 3000 /tmp/floor.txt    #and every near-floor position
//
// It runs on one thread so that a run is reproducible; FLUIDSIM_PROBE_THREADS
// sets another count, which is what the viewer uses and what a difference
// between the two would show up in.
//
// The solver only calls the hooks below when it is built with -DFLUID_PROBE,
// which nothing but this program does.
#include "../fluid.h"
#include "../scenes.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>

struct Bucket
{
	long in, out, flat;
	double dv_in, dv_out;
	double vn_out_worst;
};
static Bucket g_ground, g_wall, g_speed;
static double g_speed_vz = 0.0;
static long g_shifts = 0, g_teleports = 0;

static void tally ( Bucket & b, float vn, float dv )
{
	if ( vn < -1e-4f )      { b.in++;  b.dv_in  += dv; }
	else if ( vn > 1e-4f )
	{
		b.out++; b.dv_out += dv;
		if ( vn > b.vn_out_worst ) b.vn_out_worst = vn;
	}
	else                    b.flat++;
}

// The funnel's floor: a setspeed patch at z=1 over [-6,6]^2 that pushes what
// lands on it east, with the flat ground at z=0 underneath and around it. The
// two treat an arriving particle completely differently, so count arrivals
// separately for each.
static const float FLOOR_Z    = 5.0f;    //"arriving" is below this
static const float FLOOR_V    = 3.0f;    //...and still falling faster than this
static const float PLATE_HALF = 7.0f;    //the patch, plus a cell
static long arr_ground = 0, arr_plate = 0;
static double arr_ground_vz = 0.0, arr_plate_vz = 0.0;
static long arr_ground_in = 0;
static const double BLOWUP         = 300.0;   //units per second; see below
static const long   BLOWUP_REPORTS = 12;
static long xhist[7];

// A step may not move a particle further than its own interaction radius, or
// it crosses a neighbour's kernel unseen and nothing pushes back. Anything that
// comes out of a step above BLOWUP is past that, and the interesting question
// is what touched it while it was being taken - so the hooks tally per particle
// for the step in flight, and the first few offenders are named. A shift or a
// teleport moves a particle on purpose and is not an offence, hence ev_move.
static std::vector<int> ev_ground, ev_wall, ev_speed, ev_move;
static inline void bump ( std::vector<int> & v, unsigned long i )
{
	if ( i < v.size() ) ++v[i];
}

static void arrival ( long & n, double & vz, Vektor x, Vektor v )
{
	if ( x.z() < FLOOR_Z && v.z() < -FLOOR_V ) { n++; vz += -v.z(); }
}

void fluid_probe_ground ( unsigned long pi, Vektor x, Vektor v, Vektor n )
{
	bump ( ev_ground, pi );
	Vektor nn = n.normed();
	float vn = v*nn;
	tally ( g_ground, vn, fabsf ( 1.8f*vn ) );
	arrival ( arr_ground,arr_ground_vz,x,v );
	{
		float a = x.x();
		int b = a < -6 ? 0 : a < -2 ? 1 : a < 2 ? 2 : a < 6 ? 3 : a < 10 ? 4 : a < 15 ? 5 : 6;
		xhist[b]++;
	}
	if ( x.z() < FLOOR_Z && v.z() < -FLOOR_V
	     && fabsf ( x.x() ) <= PLATE_HALF && fabsf ( x.y() ) <= PLATE_HALF )
		arr_ground_in++;
}

void fluid_probe_control ( unsigned long pi, int kind, Vektor x, Vektor v, Vektor n )
{
	bump ( kind==2 ? ev_wall : ( kind>=3 ? ev_move : ev_speed ), pi );
	Vektor nn = n.normed();
	float vn = v*nn;
	if ( kind == 2 )
	{
		tally ( g_wall, vn, fabsf ( 1.01f*vn ) );
		return;
	}
	//The two that move a particle instead of turning it. Neither takes any
	//speed off, so there is nothing to tally - what matters is how often they
	//fire, against how often the loop they belong to ought to be turning.
	if ( kind == 3 ) { ++g_teleports; return; }
	if ( kind == 4 ) { ++g_shifts;    return; }
	tally ( g_speed, vn, fabsf ( vn ) );
	g_speed_vz += fabs ( v.z() );
	if ( x.z() < 50.0f )   //the floor patch, not the one at z=100
		arrival ( arr_plate,arr_plate_vz,x,v );
}

static void report ( const char * what, const Bucket & b, long steps )
{
	long tot = b.in+b.out+b.flat;
	if ( !tot ) { printf ( "%-22s never fired\n", what ); return; }
	printf ( "%-22s %9ld events (%7.1f/step)  into %5.1f%%  away %5.1f%%  grazing %5.1f%%\n",
	         what, tot, ( double ) tot/steps,
	         100.0*b.in/tot, 100.0*b.out/tot, 100.0*b.flat/tot );
	printf ( "%-22s   normal speed taken off: %.4g coming in, %.4g going out; "
	         "worst outgoing normal speed turned back: %.4g\n",
	         "", b.dv_in, b.dv_out, b.vn_out_worst );
}

int main ( int argc, char ** argv )
{
	const char * name = argc>1 ? argv[1] : "funnel";
	long steps = argc>2 ? atol ( argv[2] ) : 4000;
	const char * dumpname = 0;
	unsigned long wanted = 0;
	for ( int i=3;i<argc;++i )
	{
		if ( !strcmp ( argv[i],"-n" ) && i+1<argc ) wanted = strtoul ( argv[++i],0,10 );
		else dumpname = argv[i];
	}
	const Scene * sc = scene_by_name ( name );
	if ( !sc ) { printf ( "no scene %s\n",name ); return 1; }
	if ( !wanted ) wanted = sc->natural_particles;
#ifdef _OPENMP
	omp_set_num_threads ( getenv ( "FLUIDSIM_PROBE_THREADS" ) ? atoi ( getenv ( "FLUIDSIM_PROBE_THREADS" ) ) : 1 );
#endif
	//Sized from the scene, the way main() does it: the control particles share
	//the array with the water and a scene can hold more of them than of it.
	unsigned long need;
	{
		Fluid counter ( 2, 0.5f, 0,0,-9.81f, -200,-200,-200, 200,200,200 );
		counter.set_count_only ( true );
		sc->build ( counter, wanted );
		need = counter.particlecount() + 64;
	}
	Fluid f ( need, 0.5f, 0,0,-9.81f, -200,-200,-200, 200,200,200 );
	f.set_height_function ( sc->height );
	f.set_height_function_normal ( sc->height_normal );
	f.set_ground_restitution ( sc->ground_restitution );
	sc->build ( f, wanted );
	printf ( "scene %s: %lu moving particles, %lu control particles, %ld steps of dt=0.004\n",
	         sc->name, f.movingparticlecount(), f.boundaryparticlecount(), steps );

	FILE * dump = dumpname ? fopen ( dumpname,"w" ) : 0;
	//Where the water goes. A loop that is turning keeps the same water going
	//round; one that is not empties the scene into the distance, and the two
	//look the same for the first few seconds. So watch the extremes, and say at
	//the end how far out anything ever got.
	double run_rmax = 0, run_zmax = -1e30, run_vmax = 0;
	long nan_steps = 0;
	std::vector<float> px, py, pz;

	for ( long s=0; s<steps; ++s )
	{
		if ( dump )
		{
			float * varr = new float[16]; float ** va = &varr; unsigned long len = 5;
			f.get_particlearray ( va,len );
			unsigned long n = f.movingparticlecount();
			for ( unsigned long k=0;k<n;++k )
				if ( ( *va ) [3*k+2] < 8.0f )
					fprintf ( dump,"%ld %lu %.4f %.4f %.4f\n",
					          s,k, ( *va ) [3*k], ( *va ) [3*k+1], ( *va ) [3*k+2] );
			delete [] varr;
		}
		{
			const unsigned long nn = f.movingparticlecount();
			if ( ev_ground.size() != nn )
			{ ev_ground.assign ( nn,0 ); ev_wall.assign ( nn,0 ); ev_speed.assign ( nn,0 ); ev_move.assign ( nn,0 ); }
			else
			{
				std::fill ( ev_ground.begin(),ev_ground.end(),0 );
				std::fill ( ev_wall.begin(),ev_wall.end(),0 );
				std::fill ( ev_speed.begin(),ev_speed.end(),0 );
				std::fill ( ev_move.begin(),ev_move.end(),0 );
			}
		}
		f.progress ( 0.004f );
		{
			float * varr = new float[16]; float ** va = &varr; unsigned long len = 5;
			f.get_particlearray ( va,len );
			const unsigned long n = f.movingparticlecount();
			if ( px.size() != n ) { px.assign ( n,0 ); py.assign ( n,0 ); pz.assign ( n,0 ); }
			static long blown = 0;
			for ( unsigned long k=0;k<n;++k )
			{
				const float x= ( *va ) [3*k], y= ( *va ) [3*k+1], z= ( *va ) [3*k+2];
				if ( x!=x || y!=y || z!=z ) { ++nan_steps; continue; }
				const double r = sqrt ( ( double ) x*x + ( double ) y*y );
				if ( r > run_rmax ) run_rmax = r;
				if ( z > run_zmax ) run_zmax = z;
				//speed from the step just taken, ignoring the jumps a shift or a
				//teleport makes - those are not the particle moving.
				const double dx=x-px[k], dy=y-py[k], dz=z-pz[k];
				const double sp = sqrt ( dx*dx+dy*dy+dz*dz ) /0.004;
				if ( s && sp < 2000 && sp > run_vmax ) run_vmax = sp;
				if ( s && sp > BLOWUP && blown < BLOWUP_REPORTS && ev_move[k]==0 )
				{
					++blown;
					printf ( "  blow-up: step %ld particle %lu  %.2f %.2f %.2f -> %.2f %.2f %.2f"
					         "  speed %.4g   ground %d wall %d setspeed %d\n",
					         s,k, px[k],py[k],pz[k], x,y,z, sp,
					         ev_ground[k], ev_wall[k], ev_speed[k] );
				}
				px[k]=x; py[k]=y; pz[k]=z;
			}
			delete [] varr;
		}
	}
	if ( dump ) fclose ( dump );

	printf ( "\n" );
	report ( "ground", g_ground, steps );
	report ( "boundary walls", g_wall, steps );
	if ( g_wall.out )
		printf ( "%-22s   (the hook runs before the wall's own test, so \"away\" is what it declines to turn)\n", "" );
	report ( "setspeed patches", g_speed, steps );

	printf ( "\nparticles arriving at the floor (z<%g, still falling faster than %g):\n",
	         FLOOR_Z,FLOOR_V );
	printf ( "  bounced by the ground:        %ld, mean arrival speed %.3g\n",
	         arr_ground,arr_ground? arr_ground_vz/arr_ground : 0.0 );
	printf ( "    of those, %ld within the funnel's floor patch and %ld outside it\n",
	         arr_ground_in,arr_ground-arr_ground_in );
	printf ( "  caught by a setspeed patch:   %ld, mean arrival speed %.3g, all of it discarded\n",
	         arr_plate,arr_plate? arr_plate_vz/arr_plate : 0.0 );
	if ( g_speed.in+g_speed.out+g_speed.flat )
		printf ( "  vertical speed overwritten by setspeed, summed over the run: %.4g\n",
		         g_speed_vz );

	printf ( "\npatches that move a particle rather than turn it:\n" );
	printf ( "  shift    %8ld  (%.2f/step)\n", g_shifts,    ( double ) g_shifts/steps );
	printf ( "  teleport %8ld  (%.2f/step)\n", g_teleports, ( double ) g_teleports/steps );

	printf ( "\nover the whole run: furthest from the middle %.4g, highest %.4g, "
	         "fastest %.4g%s\n", run_rmax, run_zmax, run_vmax,
	         nan_steps ? "  *** NaN seen ***" : "" );

	printf ( "\nground bounces by x (the funnel's conveyor and its wall are at x=15):\n" );
	{
		static const char * lab[7] = { "x<-6","-6..-2","-2..2","2..6","6..10","10..15","x>15" };
		for ( int b=0;b<7;++b ) printf ( "  %-8s %8ld\n", lab[b], xhist[b] );
	}
	return 0;
}
