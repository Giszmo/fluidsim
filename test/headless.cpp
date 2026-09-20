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
#include "../fluid.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

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

int main ( int argc, char ** argv )
{
	int   scenario = ( argc > 1 ) ? atoi ( argv[1] )         : 0;
	float simtime  = ( argc > 2 ) ? ( float ) atof ( argv[2] ) : 3.0f;
	float dt       = ( argc > 3 ) ? ( float ) atof ( argv[3] ) : 0.004f;
	int   steps    = ( int ) ( simtime/dt + 0.5f );

	Fluid f ( 200000, 1.0f, 0, 0, -9.81f, -64, -64, -64, 64, 64, 64 );
	f.set_height_function ( flat_z );
	f.set_height_function_normal ( flat_n );
	add_plate ( f );
	add_liquid ( f, ( scenario == 0 ) ? -2.0f : HALF + 2.0f, -2.0f, 8.0f, 4.0f );

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
