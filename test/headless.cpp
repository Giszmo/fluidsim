// Headless regression probe for the solver in fluid.h.
//
// Runs the fluid without a window and prints an FNV-1a checksum over every
// moving particle's position, so two builds can be compared bit for bit:
//
//     make check
//
// Scenario 0 drops a tetrahedron of liquid onto the centre of a square plate,
// scenario 1 drops it beside the plate. Both are deterministic: same build,
// same arguments, same checksum.
//
// Scenario 2 is a block big enough to be spread over every core. With it,
//
//     ./headless cmp
//
// runs the solver once on one thread and once on all of them and reports how
// far the two end states differ, which is the check that matters for the
// parallel passes: they may not reproduce the serial sum bit for bit, because
// the per-thread buffers are added up in a different order, but they must not
// differ by more than float rounding.
#include "../fluid.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

// No terrain: the height function sits far below everything.
static float  flat_z ( float, float ) { return -1.0e6f; }
static Vektor flat_n ( float, float ) { return Vektor ( 0,0,1 ); }

static const float PLATE_Z = 0.0f;
static const float HALF    = 6.0f;

// A square plate at z=PLATE_Z spanning [-HALF,HALF]^2, outward normal +z.
static void add_plate ( Fluid & f )
{
	float t[18];
	int c = 0;
	t[c++]=-HALF; t[c++]=-HALF; t[c++]=PLATE_Z;
	t[c++]= HALF; t[c++]=-HALF; t[c++]=PLATE_Z;
	t[c++]=-HALF; t[c++]= HALF; t[c++]=PLATE_Z;
	t[c++]= HALF; t[c++]=-HALF; t[c++]=PLATE_Z;
	t[c++]= HALF; t[c++]= HALF; t[c++]=PLATE_Z;
	t[c++]=-HALF; t[c++]= HALF; t[c++]=PLATE_Z;
	f.trianglelist2boundary ( t, 0, 1 );
}

// A tetrahedron of liquid with its low corner at (cx,cy,cz) and edge length s.
static void add_liquid ( Fluid & f, float cx, float cy, float cz, float s )
{
	float t[12] = { cx,     cy,     cz,
	                cx + s, cy,     cz,
	                cx,     cy + s, cz,
	                cx,     cy,     cz + s };
	f.tetraederlist2liquid ( t, 0, 0 );
}

static void build ( Fluid & f, int scenario )
{
	f.set_height_function ( flat_z );
	f.set_height_function_normal ( flat_n );
	if ( scenario != 2 )
	{
		add_plate ( f );
		add_liquid ( f, ( scenario == 0 ) ? -2.0f : HALF + 2.0f, -2.0f, 8.0f, 4.0f );
	}
	else
	{
		add_liquid ( f, -14.0f, -14.0f, 8.0f, 28.0f );
	}
}

// Run one scenario to completion and hand back every particle position.
static unsigned long run ( int scenario, float simtime, float dt, int threads, float ** out )
{
#ifdef _OPENMP
	if ( threads > 0 )
		omp_set_num_threads ( threads );
#else
	( void ) threads;
#endif
	float half = ( scenario == 2 ) ? 192.0f : 64.0f;
	Fluid f ( 200000, 1.0f, 0, 0, -9.81f, -half,-half,-half, half,half,half );
	build ( f, scenario );

	int steps = ( int ) ( simtime/dt + 0.5f );
	for ( int s = 0; s < steps; ++s )
		f.progress ( dt );

	float * varr = new float[16];
	float ** va  = &varr;
	unsigned long len = 5;
	f.get_particlearray ( va, len );
	unsigned long n = f.movingparticlecount();
	*out = new float[3*n];
	memcpy ( *out, *va, 3*n*sizeof ( float ) );
	delete [] varr;
	return n;
}

// One thread against all of them, on the same scenario.
static int compare_thread_counts()
{
	float * a; float * b;
	int nt = 1;
#ifdef _OPENMP
	nt = omp_get_max_threads();
#endif
	unsigned long na = run ( 2, 0.4f, 0.004f, 1, &a );
	unsigned long nb = run ( 2, 0.4f, 0.004f, nt, &b );
	if ( na != nb )
	{
		printf ( "FAIL     particle counts differ: %lu vs %lu\n", na, nb );
		return 1;
	}
	double worst = 0.0;
	for ( unsigned long i = 0; i < na; ++i )
	{
		double dx = a[3*i]-b[3*i], dy = a[3*i+1]-b[3*i+1], dz = a[3*i+2]-b[3*i+2];
		double d = sqrt ( dx*dx + dy*dy + dz*dz );
		if ( d > worst ) worst = d;
	}
	delete [] a; delete [] b;
	printf ( "CMP      particles=%lu threads=1 vs %d  max displacement %.3e "
	         "(particle diameter 1.0)\n", na, nt, worst );
	if ( worst > 1.0e-3 )
	{
		printf ( "FAIL     the two thread counts do not agree\n" );
		return 1;
	}
	printf ( "OK       within float rounding\n" );
	return 0;
}

int main ( int argc, char ** argv )
{
	if ( argc > 1 && strcmp ( argv[1], "cmp" ) == 0 )
		return compare_thread_counts();

	int   scenario = ( argc > 1 ) ? atoi ( argv[1] )         : 0;
	float simtime  = ( argc > 2 ) ? ( float ) atof ( argv[2] ) : 3.0f;
	float dt       = ( argc > 3 ) ? ( float ) atof ( argv[3] ) : 0.004f;
	int   steps    = ( int ) ( simtime/dt + 0.5f );

	Fluid f ( 200000, 1.0f, 0, 0, -9.81f, -64, -64, -64, 64, 64, 64 );
	build ( f, scenario );

	float * varr = new float[16];
	float ** va  = &varr;
	unsigned long len = 5;

	printf ( "scenario=%d boundary=%lu moving=%lu dt=%g steps=%d simtime=%gs\n",
	         scenario, f.boundaryparticlecount(), f.movingparticlecount(), dt, steps, simtime );

	for ( int s = 1; s <= steps; ++s )
		f.progress ( dt );

	f.get_particlearray ( va, len );
	unsigned long n = f.movingparticlecount();

	float minx=1e9f, maxx=-1e9f, minz=1e9f, maxz=-1e9f;
	int below = 0;
	for ( unsigned long i = 0; i < n; ++i )
	{
		float x = ( *va ) [3*i], z = ( *va ) [3*i+2];
		if ( x < minx ) minx = x;
		if ( x > maxx ) maxx = x;
		if ( z < minz ) minz = z;
		if ( z > maxz ) maxz = z;
		if ( z < PLATE_Z - 0.5f ) ++below;
	}

	unsigned long long h = 1469598103934665603ULL;   // FNV-1a over the raw float bits
	for ( unsigned long i = 0; i < 3*n; ++i )
	{
		unsigned int bits;
		float v = ( *va ) [i];
		memcpy ( &bits, &v, 4 );
		h ^= bits;
		h *= 1099511628211ULL;
	}

	printf ( "RESULT   particles=%lu through_plate=%d x=[%.2f,%.2f] z=[%.2f,%.2f] (plate x=[%.1f,%.1f] z=%.1f)\n",
	         n, below, minx, maxx, minz, maxz, -HALF, HALF, PLATE_Z );
	printf ( "CHECKSUM %016llx\n", h );

	delete [] varr;
	return 0;
}
