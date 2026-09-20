// Headless check on the marching cubes surface, run by
//
//     make check
//
// on the scene the viewer actually starts with: the parabolic ground, the jet
// in the middle and a drop of liquid falling into it. Every step it pulls the
// surface out of the solver and asserts the three things that were broken
// before:
//
//  - no triangle is longer than a cube. A vertex sits in the cell it belongs
//    to, stepped at most one cell towards a neighbour, so two corners of one
//    cube are at most two cells apart: 2*sqrt(3)*cellsize. Anything longer is
//    the spike that used to shoot across the whole domain when a stray
//    particle left the box and its cell index wrapped into the middle of the
//    pool.
//  - every index points at a vertex that exists.
//  - triangles come to about twice the vertices, which is what a closed
//    surface looks like; a mesh full of dropped triangles falls well below it.
#include "../fluid.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

static float  sink       ( float a, float b ) { return .05f* ( a*a+b*b ); }
static Vektor sinknormal ( float a, float b ) { return Vektor ( -.1f*a,-.1f*b,1 ); }

static const float PARTICLESIZE = 0.5f;
static const float CELLSIZE     = 4.0f*PARTICLESIZE;
// Two cells, corner to corner: the longest edge a single cube can produce.
static const float MAXEDGE      = 2.0f*1.7320508f*CELLSIZE;

int main ( int argc, char ** argv )
{
	const int  steps  = argc>1 ? atoi ( argv[1] ) : 900;
	const long wanted = argc>2 ? atol ( argv[2] ) : 20000;

	Fluid f ( 200000, PARTICLESIZE, 0,0,-9.81f, -200,-200,-200, 200,200,200 );
	f.set_height_function ( sink );
	f.set_height_function_normal ( sinknormal );

	// The 4x4 plate at the origin that throws fluid straight up.
	float t[18]; int c=0;
	t[c++]=-2;t[c++]=-2;t[c++]=0;  t[c++]=-2;t[c++]= 2;t[c++]=0;  t[c++]= 2;t[c++]=-2;t[c++]=0;
	t[c++]= 2;t[c++]=-2;t[c++]=0;  t[c++]=-2;t[c++]= 2;t[c++]=0;  t[c++]= 2;t[c++]= 2;t[c++]=0;
	f.trianglelist2setspeed ( 0,0,70,t,0,1 );

	// The drop above it, scaled to the particle count asked for.
	{
		const float v0[3] = { -5,-5,55 };
		const float e[3][3] = { { 15,0,5 }, { 0,15,5 }, { 0,0,10 } };
		const float s = powf ( ( float ) wanted/1875.0f, 1.0f/3.0f );
		float tl[12];
		for ( int i=0;i<3;++i ) tl[i]=v0[i];
		for ( int j=0;j<3;++j ) for ( int i=0;i<3;++i ) tl[3+3*j+i]=v0[i]+e[j][i]*s;
		f.tetraederlist2liquid ( tl,0,0 );
	}
	printf ( "surface  scene=jet particles=%lu steps=%d dt=0.004\n",
	         f.movingparticlecount(), steps );

	float * va = new float[16];
	float * na = new float[16];
	unsigned int * ti = new unsigned int[16];
	unsigned long vnlen=5, tilen=5, vncount=0;

	float  worstedge = 0.0f;
	double worstratio = 1e9;
	long   badindex = 0;
	int    peaktris = 0;

	for ( int s=0; s<steps; ++s )
	{
		f.progress ( 0.004f );
		const int tc = f.get_surfacegrid ( &va,&na,vnlen,vncount,&ti,tilen );
		if ( tc > peaktris ) peaktris = tc;
		if ( tc == 0 )
			continue;
		for ( int i=0;i<3*tc;++i )
			if ( ti[i] == 0 || ti[i] >= vncount )
				++badindex;
		if ( badindex )
			break;
		for ( int i=0;i<tc;++i )
			for ( int k=0;k<3;++k )
			{
				const unsigned int p=ti[3*i+k], q=ti[3*i+ ( k+1 ) %3];
				const float dx=va[3*p]-va[3*q], dy=va[3*p+1]-va[3*q+1], dz=va[3*p+2]-va[3*q+2];
				const float d=sqrtf ( dx*dx+dy*dy+dz*dz );
				if ( d>worstedge ) worstedge=d;
			}
		// Only judge the ratio once there is a real surface; a handful of
		// single cells is all boundary and says nothing.
		if ( tc > 500 )
		{
			const double r = ( double ) tc / ( double ) ( vncount-1 );
			if ( r < worstratio ) worstratio = r;
		}
	}

	printf ( "RESULT   peak triangles=%d  longest edge=%.2f (limit %.2f)  "
	         "worst triangles-per-vertex=%.2f  bad indices=%ld\n",
	         peaktris, worstedge, MAXEDGE, worstratio, badindex );

	int fail = 0;
	if ( badindex )                 { printf ( "FAIL     %ld indices point outside the vertex array\n",badindex ); fail=1; }
	if ( worstedge > MAXEDGE )      { printf ( "FAIL     a triangle is longer than a cube: %.2f > %.2f\n",worstedge,MAXEDGE ); fail=1; }
	if ( worstratio < 1.5 )         { printf ( "FAIL     only %.2f triangles per vertex; the mesh is full of holes\n",worstratio ); fail=1; }
	if ( peaktris < 1000 )          { printf ( "FAIL     the surface never grew past %d triangles\n",peaktris ); fail=1; }
	if ( !fail )
		printf ( "OK       closed, bounded, fully indexed\n" );
	return fail;
}
