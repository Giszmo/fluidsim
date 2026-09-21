#include "scenes.h"
#include <cmath>
#include <cstdio>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ---------------------------------------------------------------- helpers --

//One tetrahedron of liquid, given as a corner and three edges, scaled about
//that corner so that it holds about `wanted` particles instead of `natural`.
//tetraeder2particle() subdivides until each piece holds one particle, so the
//count is proportional to the volume and the edges scale with its cube root.
static void liquid_tetra ( Fluid & f, const float v0[3], const float e[3][3],
                           unsigned long wanted, unsigned long natural )
{
	float s = powf ( ( float ) wanted / ( float ) natural, 1.0f/3.0f );
	float t[12];
	int i,j;
	for ( i=0;i<3;++i )
		t[i] = v0[i];
	for ( j=0;j<3;++j )
		for ( i=0;i<3;++i )
			t[3+3*j+i] = v0[i] + e[j][i]*s;
	f.tetraederlist2liquid ( t,0,0 );
}

//The four control-particle kinds all take a triangle list, so a rectangle is
//two triangles wound so that the normal comes out where the caller wants it.
//n = (b-a) x (c-a), and a particle is affected when it is on the negative
//side of that normal.
static void rect ( float * t, int & c,
                   Vektor a, Vektor b, Vektor cc, Vektor d )
{
	Vektor v[6] = { a,b,cc, a,cc,d };
	for ( int i=0;i<6;++i )
	{
		t[c++]=v[i].x(); t[c++]=v[i].y(); t[c++]=v[i].z();
	}
}

//How far past a limit a value is, signed, and 0 between them. Both scenes with
//a rim build it out of this.
static float over ( float v, float lo, float hi )
{
	if ( v > hi ) return v-hi;
	if ( v < lo ) return v-lo;
	return 0;
}

// -------------------------------------------------------- scene: fountain --
// The scene this program has always started in: a parabolic bowl, a small
// plate at the bottom that throws anything crossing it straight up, and a
// tetrahedron of liquid dropped on it.

static float  fountain_h ( float a, float b ) { return .05f* ( a*a+b*b ); }
static Vektor fountain_n ( float a, float b ) { return Vektor ( -.1f*a,-.1f*b,1 ); }

static void fountain_build ( Fluid & f, unsigned long wanted )
{
	float t[18];
	int c = 0;
	rect ( t,c, Vektor ( -2,-2,0 ),Vektor ( -2,2,0 ),Vektor ( 2,2,0 ),Vektor ( 2,-2,0 ) );
	f.trianglelist2setspeed ( 0,0,70,t,0,1 );

	const float v0[3] = { -5,-5,55 };
	const float e[3][3] = { { 15, 0, 5 }, { 0, 15, 5 }, { 0, 0, 10 } };
	liquid_tetra ( f,v0,e,wanted,1875 );
}

// ---------------------------------------------------------- scene: funnel --
// The scene that was in main() in 2005 and has been sitting there commented
// out ever since. A closed loop: liquid falls through a funnel, is dropped out
// of its throat onto the floor, pushed east along it, lifted from a wall back
// up to z=100, and started downwards again by a plate up there.
//
// Two things about the loop are not the 2005 ones, because the 2005 ones do
// not close it:
//
// The return is a shift, not a teleport. A teleport sends every particle to
// one target point, and every control particle of the patch shares that
// target - so the whole arriving stream is squeezed into a fraction of a cell,
// whatever the patch's size. The wall here is 30 by 2 and the 2005 teleport
// mapped it onto about 3 by 0.2: the water came out as a blast that threw
// particles clear over the funnel's mouth, and they then slid away over a flat
// frictionless ground for good. Measured over 24 simulated seconds: 68% of the
// scene ended up outside any part of it that can return water. A shift moves
// the whole patch rigidly, so the stream arrives in the shape it left in.
//
// The ground is a basin rather than an infinite flat plane. Splashes leave the
// conveyor whatever the return path does, and on a flat frictionless ground
// anything that leaves is gone. The rim is part of the height function, so it
// is tested against every particle every step.

//Flat under the funnel, the conveyor and the wall; beyond that the ground turns
//up into a quadratic rim that rolls water back in. Same shape as the terrain
//scene's banks, per axis rather than radial, and flat where it meets them so
//there is no lip to bounce off. At 0.6 a particle leaving the conveyor at 30 -
//three times the conveyor's own speed - has climbed to a stop by r=31.
#define FUNNEL_FLAT 22.0f
#define FUNNEL_RIM  0.6f

static float funnel_h ( float x, float y )
{
	float ox = over ( x,-FUNNEL_FLAT,FUNNEL_FLAT ), oy = over ( y,-FUNNEL_FLAT,FUNNEL_FLAT );
	return FUNNEL_RIM* ( ox*ox + oy*oy );
}
static Vektor funnel_n ( float x, float y )
{
	float ox = over ( x,-FUNNEL_FLAT,FUNNEL_FLAT ), oy = over ( y,-FUNNEL_FLAT,FUNNEL_FLAT );
	return Vektor ( -2*FUNNEL_RIM*ox,-2*FUNNEL_RIM*oy,1 );
}

//The rings of the funnel, radius j*j at z=15+8j, j from 1 to 5. The viewer
//draws the same numbers into its TRICHTER display list.
void funnel_ring ( int j, int i, float * out )
{
	float a = ( float ) i/10*2*M_PI;
	out[0] = cosf ( a ) * ( float ) ( j*j );
	out[1] = sinf ( a ) * ( float ) ( j*j );
	out[2] = 15.0f + 8.0f* ( float ) j;
}

//The wall the water is lifted from, and where it is put down: x=15 back to
//x=0, z=0..2 up to z=98..100, which is the plate the 2005 scene dropped it
//from. y is untouched, so a stream arrives where it left.
#define FUNNEL_WALL_X   15.0f
#define FUNNEL_WALL_Y   15.0f
#define FUNNEL_TOP_Z    98.0f

static void funnel_build ( Fluid & f, unsigned long wanted )
{
	static float t[9*1700];
	int c = 0;
	int i,j;
	for ( j=2;j<6;j++ )
	{
		for ( i=0;i<10;i++ )
		{
			funnel_ring ( j,  i,   t+c ); c+=3;
			funnel_ring ( j,  i+1, t+c ); c+=3;
			funnel_ring ( j-1,i,   t+c ); c+=3;

			funnel_ring ( j,  i+1, t+c ); c+=3;
			funnel_ring ( j-1,i+1, t+c ); c+=3;
			funnel_ring ( j-1,i,   t+c ); c+=3;
		}
	}
	f.trianglelist2boundary ( t,0,c/9-1 );

	//the funnel's throat: drop whatever gets through it ten units
	c = 0;
	rect ( t,c, Vektor ( -3,-3,23 ),Vektor ( 3,-3,23 ),Vektor ( 3,3,23 ),Vektor ( -3,3,23 ) );
	f.trianglelist2shift ( 0,0,-10,t,0,1 );

	//the floor: push everything that lands on it east, all the way to the wall.
	//In 2005 this stopped at x=6 while the wall it feeds is at x=15, so water
	//that came down past the end of it was never pushed anywhere - it bounced
	//on the ground at 0.8 and stayed where it fell. 96% of all ground bounces
	//in the scene were beyond x=6.
	c = 0;
	rect ( t,c, Vektor ( -FUNNEL_WALL_Y,-FUNNEL_WALL_Y,1 ),Vektor ( FUNNEL_WALL_X,-FUNNEL_WALL_Y,1 ),
	       Vektor ( FUNNEL_WALL_X,FUNNEL_WALL_Y,1 ),Vektor ( -FUNNEL_WALL_Y,FUNNEL_WALL_Y,1 ) );
	f.trianglelist2setspeed ( 10,0,0,t,0,1 );

	//the wall it arrives at: back up to the top, rigidly
	c = 0;
	rect ( t,c, Vektor ( FUNNEL_WALL_X,FUNNEL_WALL_Y,2 ),Vektor ( FUNNEL_WALL_X,FUNNEL_WALL_Y,0 ),
	       Vektor ( FUNNEL_WALL_X,-FUNNEL_WALL_Y,0 ),Vektor ( FUNNEL_WALL_X,-FUNNEL_WALL_Y,2 ) );
	f.trianglelist2shift ( -FUNNEL_WALL_X,0,FUNNEL_TOP_Z,t,0,1 );

	//and start it falling again, into the funnel. The water arrives up here
	//still carrying the conveyor's 10 east, which over the 43 units it has to
	//fall to the funnel's mouth would carry it 30 sideways and clean past it;
	//this is what takes that off. It has to be wider than the 10x10 of 2005 for
	//the same reason: the arrival is a stream 30 long, not a point.
	c = 0;
	rect ( t,c, Vektor ( -8,-FUNNEL_WALL_Y-2,FUNNEL_TOP_Z-1 ),Vektor ( 8,-FUNNEL_WALL_Y-2,FUNNEL_TOP_Z-1 ),
	       Vektor ( 8,FUNNEL_WALL_Y+2,FUNNEL_TOP_Z-1 ),Vektor ( -8,FUNNEL_WALL_Y+2,FUNNEL_TOP_Z-1 ) );
	f.trianglelist2setspeed ( 0,0,-1,t,0,1 );

	const float v0[3] = { -5,-5,55 };
	const float e[3][3] = { { 10, 0, 5 }, { 0, 10, 5 }, { 0, 0, 10 } };
	liquid_tetra ( f,v0,e,wanted,833 );
}

// --------------------------------------------------------- scene: terrain --
// A channel tilted along x, rippled along and across it, with the low end
// wired back to the high end. Water poured in at the top runs down, fills the
// dips on the way as puddles, finds the gullies as rivers, and comes back
// round, so the amount of water in the channel is constant and the run never
// ends.
//
//     h(x,y) = -slope*x                      the tilt
//              - wave * sin(2*pi*RIPPLES*x/TERRAIN_LEN)   dips along it
//              + CROSS_AMP * cos(2*pi*GULLIES*y/TERRAIN_WID)  gullies across
//              + BANK * y*y                  banks, so nothing runs off sideways
//
// The sine along x has a whole number of periods over the channel, and both
// ends of the channel fall on one of its steepest downhill stretches, so the
// water drains freely into the gate at the low end instead of having to climb
// a ripple to reach it, and arrives back at the top already running. The
// periodicity is what matters most though:
//
//     h(x-TERRAIN_LEN,y) = h(x,y) + slope*TERRAIN_LEN
//
// exactly. That is what makes the loop seamless: the shift patch at the low
// end moves a particle back by (-TERRAIN_LEN, 0, +slope*TERRAIN_LEN), which
// leaves it the same height above the ground as it was, with the velocity it
// arrived with. A shift, unlike a teleport, keeps the particle's offset within
// the patch, so a river arrives back in its own gully.

#define TERRAIN_LEN 200.0f
#define TERRAIN_WID 50.0f
#define TERRAIN_X1  ( TERRAIN_LEN/2 )   //the low end
#define TERRAIN_X0  ( -TERRAIN_LEN/2 )  //the high end
#define RIPPLES     4
#define GULLIES     2
#define CROSS_AMP   2.0f
#define BANK        0.006f
//Where the channel ends and the ground starts climbing out of it, and how
//hard. The rim is quadratic and starts flat, so there is no kink at the lip
//and nothing to bounce off; by the time a particle is a few units up it it is
//climbing a 1:7 wall, which is more than the speed anything in the channel
//arrives with. This is what keeps the water in the scene: the terrain, which
//every particle is tested against every step, rather than more control
//particles, which only act within a cell of themselves.
#define RIM_X0      ( TERRAIN_X0-2.0f )
#define RIM_X1      ( TERRAIN_X1+2.0f )
#define RIM_Y       ( TERRAIN_WID*0.65f )
#define RIM         0.47f

static float terrain_slope = 0.12f;
static float terrain_wave  = 1.8f;

//Anything not positive means "leave it at the scene's own value", so that
//the defaults live here and not in two places.
void terrain_set ( float slope, float waviness )
{
	if ( slope > 0 )
		terrain_slope = slope;
	if ( waviness > 0 )
		terrain_wave = waviness;
}

static const float TERRAIN_KX = 2.0f*( float ) M_PI*RIPPLES/TERRAIN_LEN;
static const float TERRAIN_KY = 2.0f*( float ) M_PI*GULLIES/TERRAIN_WID;

static float terrain_h ( float x, float y )
{
	float ox = over ( x,RIM_X0,RIM_X1 ), oy = over ( y,-RIM_Y,RIM_Y );
	return -terrain_slope*x
	       - terrain_wave*sinf ( TERRAIN_KX*x )
	       + CROSS_AMP*cosf ( TERRAIN_KY*y )
	       + BANK*y*y
	       + RIM* ( ox*ox + oy*oy );
}

static Vektor terrain_n ( float x, float y )
{
	//(-dh/dx, -dh/dy, 1), the same convention the bowl's normal uses.
	float ox = over ( x,RIM_X0,RIM_X1 ), oy = over ( y,-RIM_Y,RIM_Y );
	return Vektor (
	           terrain_slope + terrain_wave*TERRAIN_KX*cosf ( TERRAIN_KX*x ) - 2*RIM*ox,
	           CROSS_AMP*TERRAIN_KY*sinf ( TERRAIN_KY*y ) - 2*BANK*y - 2*RIM*oy,
	           1 );
}

//The patch that closes the loop: it spans the whole channel, out to the rim,
//and follows the ground - one rectangle per strip across it, each flat
//bottomed so that its normal comes out exactly along -x. It does not have to
//be tall, because it is not a wall: it catches whatever crosses it, and the
//rim catches whatever gets past it and rolls it back within reach.
#define TERRAIN_STRIPS 24
#define TERRAIN_GATE_H 18.0f

static void terrain_gate ( Fluid & f, float x )
{
	static float t[9*2*TERRAIN_STRIPS];
	int c = 0;
	for ( int i=0;i<TERRAIN_STRIPS;++i )
	{
		float y0 = -RIM_Y + 2*RIM_Y* ( float ) i    /TERRAIN_STRIPS;
		float y1 = -RIM_Y + 2*RIM_Y* ( float ) ( i+1 ) /TERRAIN_STRIPS;
		float z0 = terrain_h ( x,y0 ) < terrain_h ( x,y1 ) ? terrain_h ( x,y0 ) : terrain_h ( x,y1 );
		z0 -= 2.0f;
		float z1 = z0 + TERRAIN_GATE_H;
		//wound so that (b-a)x(c-a) points along -x: a particle downstream of
		//the plane is on the negative side of it and is the one that is caught.
		rect ( t,c, Vektor ( x,y0,z0 ),Vektor ( x,y0,z1 ),Vektor ( x,y1,z1 ),Vektor ( x,y1,z0 ) );
	}
	f.trianglelist2shift ( -TERRAIN_LEN,0,terrain_slope*TERRAIN_LEN,t,0,c/9-1 );
}

//A box of liquid. Kuhn's decomposition: the six tetrahedra (0, e_p1,
//e_p1+e_p2, 1) over the six orderings of the three axes tile the unit cube
//exactly once, so the particles come out evenly and none is placed twice.
static void liquid_box ( Fluid & f, const float lo[3], const float hi[3] )
{
	static const int perm[6][3] = { {0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0} };
	for ( int p=0;p<6;++p )
	{
		float tet[12];
		float cur[3] = { lo[0],lo[1],lo[2] };
		for ( int k=0;k<3;++k ) tet[k] = cur[k];
		for ( int s=0;s<3;++s )
		{
			cur[perm[p][s]] = hi[perm[p][s]];
			for ( int k=0;k<3;++k ) tet[3* ( s+1 ) +k] = cur[k];
		}
		f.tetraederlist2liquid ( tet,0,0 );
	}
}

#define TERRAIN_SEGMENTS 24

static void terrain_build ( Fluid & f, unsigned long wanted )
{
	terrain_gate ( f,TERRAIN_X1 );

	//A sheet of water over the whole channel rather than a block at the top.
	//A block has to drain the length of the channel before anything
	//downstream of it is wet, and the first dip it meets holds tens of
	//thousands of particles, so most runs would spend all their water on one
	//puddle at the top. Spread out, every dip starts with its share and
	//whatever is over the brim runs.
	//
	//The sheet follows the ground, because the ground drops 43 over the length
	//of the channel and a flat slab would be a waterfall at one end. One box
	//per segment, each sitting just clear of the highest ground under it.
	const float wide = TERRAIN_WID*0.45f;              //half width
	const float x0 = TERRAIN_X0 + 4.0f;
	const float x1 = TERRAIN_X1 - 4.0f;
	float deep = 0.2f* ( float ) wanted / ( 2*wide* ( x1-x0 ) );
	if ( deep < 0.6f )
		deep = 0.6f;
	for ( int i=0;i<TERRAIN_SEGMENTS;++i )
	{
		float a = x0 + ( x1-x0 ) * ( float ) i    /TERRAIN_SEGMENTS;
		float b = x0 + ( x1-x0 ) * ( float ) ( i+1 ) /TERRAIN_SEGMENTS;
		float top = -1e30f;
		for ( int u=0;u<=8;++u )
			for ( int v=0;v<=8;++v )
			{
				float h = terrain_h ( a + ( b-a ) *u/8, -wide + 2*wide*v/8 );
				if ( h > top ) top = h;
			}
		const float lo[3] = { a,-wide,top+0.5f };
		const float hi[3] = { b, wide,top+0.5f+deep };
		liquid_box ( f,lo,hi );
	}
}

// ------------------------------------------------------------- the table --

static const Scene scenes[] =
{
	{
		"fountain",
		"a jet in a parabolic bowl throws a tetrahedron of liquid back up (the default)",
		fountain_h, fountain_n, fountain_build, 1875,
		false, -45,45,45, -45,45,45,
		false, 10.0f, 0.0f, 0.8f
	},
	{
		"funnel",
		"the 2005 scene: through a funnel, along the floor, lifted back to the top",
		funnel_h, funnel_n, funnel_build, 833,
		true, -30,30,30, -30,30,30,
		true, 110.0f, 30.0f, 0.8f
	},
	{
		"terrain",
		"a tilted rippled channel wired back on itself: rivers and puddles, forever",
		terrain_h, terrain_n, terrain_build, 25000,
		true, RIM_X0-6,RIM_X1+6,160, -RIM_Y-6,RIM_Y+6,48,
		false, 150.0f, 62.0f, 0.25f
	},
};
static const int scene_count = sizeof ( scenes ) /sizeof ( scenes[0] );

const Scene * scene_by_name ( const char * name )
{
	for ( int i=0;i<scene_count;++i )
		if ( !strcmp ( name,scenes[i].name ) )
			return &scenes[i];
	return 0;
}

const Scene * scene_default()
{
	return &scenes[0];
}

void scene_list()
{
	for ( int i=0;i<scene_count;++i )
		printf ( "  %-10s %s\n", scenes[i].name, scenes[i].summary );
}
