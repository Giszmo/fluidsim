// What actually stops a particle in a scene.
//
// Builds a scene exactly as the viewer does, runs it without a window, and
// counts every correction the solver applies to a particle's velocity, split
// by whether the particle was moving into the surface or already moving away
// from it. The second number is the interesting one: a correction applied to
// a particle that is on its way out takes speed off it for nothing, and with
// a low restitution it is what makes water stick to a wall.
//
//     make test/probe
//     ./test/probe funnel 4000
//     ./test/probe funnel 3000 /tmp/floor.txt    #and every near-floor position
//
// The solver only calls the hooks below when it is built with -DFLUID_PROBE,
// which nothing but this program does.
#include "../fluid.h"
#include "../scenes.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

struct Bucket
{
	long in, out, flat;
	double dv_in, dv_out;
	double vn_out_worst;
};
static Bucket g_ground, g_wall, g_speed;
static double g_speed_vz = 0.0;

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
static long xhist[7];

static void arrival ( long & n, double & vz, Vektor x, Vektor v )
{
	if ( x.z() < FLOOR_Z && v.z() < -FLOOR_V ) { n++; vz += -v.z(); }
}

void fluid_probe_ground ( unsigned long, Vektor x, Vektor v, Vektor n )
{
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

void fluid_probe_control ( unsigned long, int kind, Vektor x, Vektor v, Vektor n )
{
	Vektor nn = n.normed();
	float vn = v*nn;
	if ( kind == 2 )
	{
		tally ( g_wall, vn, fabsf ( 1.01f*vn ) );
		return;
	}
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
	const Scene * sc = scene_by_name ( name );
	if ( !sc ) { printf ( "no scene %s\n",name ); return 1; }
#ifdef _OPENMP
	omp_set_num_threads ( 1 );
#endif
	Fluid f ( 150000, 0.5f, 0,0,-9.81f, -200,-200,-200, 200,200,200 );
	f.set_height_function ( sc->height );
	f.set_height_function_normal ( sc->height_normal );
	f.set_ground_restitution ( sc->ground_restitution );
	sc->build ( f, sc->natural_particles );
	printf ( "scene %s: %lu moving particles, %lu control particles, %ld steps of dt=0.004\n",
	         sc->name, f.movingparticlecount(), f.boundaryparticlecount(), steps );

	FILE * dump = argc>3 ? fopen ( argv[3],"w" ) : 0;
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
		f.progress ( 0.004f );
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

	printf ( "\nground bounces by x (the funnel's floor patch ends at x=6, its teleport wall is at x=15):\n" );
	{
		static const char * lab[7] = { "x<-6","-6..-2","-2..2","2..6","6..10","10..15","x>15" };
		for ( int b=0;b<7;++b ) printf ( "  %-8s %8ld\n", lab[b], xhist[b] );
	}
	return 0;
}
