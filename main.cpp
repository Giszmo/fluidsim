#include "fluid.h"
#include "bicubic_bezier_surface.h"
#include "simthread.h"
#include "scenes.h"
#include <math.h>

#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;

#ifdef WIN32
#include <windows.h>
#include "glut.h"
#include "glati.h"
#include "wglati.h"
#include "ATIExtensions.h"
#else
#define DWORD unsigned long
#define WORD unsigned char
#define LONG long int
#include <GL/glut.h>
#include "shaders.h"
#endif
#define WIDTH  1024
#define HEIGHT 768
#define TITLE  "Waterworld :)"
#define M_PI 3.14159265358979323846

extern "C" void kbf ( unsigned char key,int x, int y );
extern "C" void kbf2 ( int key,int x, int y );
extern "C" void DisplayMain ( void );
extern "C" void zeitgeber ( int value );
extern "C" void mousemove ( int x, int y );
extern "C" void ReshapeMain ( int w, int h );

struct BMP_FILEHEADER
{
	WORD bfType; //Typ - muss BM=19778 sein
	DWORD bfSize; //Gr?e der Datei in Bytes
	WORD bfReserved1; //muss 0 sein
	WORD bfReserved2; //muss 0 sein
	DWORD bfOffBits; //Start der Bitmapdaten als offset von der BITMAPFILEHEADER Struktur
} mybmpfileheader;

struct BMP_INFOHEADER
{
	DWORD biSize; //gr?e der Bitmapinfo struct in bytes
	LONG biWidth; //breite des Bitmap in pixel
	LONG biHeight; //h?e des Bitmap in pixel
	WORD biPlanes; //muss 1 sein
	WORD biBitCount; //hier 24
	DWORD biCompression; //hier BI_RGB=0
	DWORD biSizeImage; //gr?e des Bildes in Bytes, kann wenn ohne Kompression auch 0 sein
	LONG biXPelsPerMeter; //pixel pro meter ->unwichtig
	LONG biYPelsPerMeter; //s.o.
	DWORD biClrUsed; //hier 0
	DWORD biClrImportant; //hier auch 0
} mybmpinfoheader;


void initMain ( void );
int DrawTextXY ( float x,float y,float scale,unsigned char *s );
void initCallLists ( void );
void makeMenu();
void Menufunc ( int value );

int Menu;
int winIdMain;
float distanz=10;
float xt;
float yt;
float speed=0.004f;
bool paused = true;
bool shownormals= false;

unsigned short framecount=0;
unsigned long trianglecount=0;
bool showcells = false;
bool showparticles = true;
//PN triangles on the surface. In 2005 this was TRUFORM, an ATI extension that
//no driver has any more; today it is the OpenGL 4 tessellation stage doing the
//same maths, see shaders.cpp. Off automatically where there is no OpenGL 4.
bool tessellation_on = true;
//How long one tessellated segment should come out on screen, in pixels, and how
//far a single triangle may be subdivided. Together they make the subdivision
//adaptive: close water is smooth, distant spray costs nothing.
float tess_pixels_per_segment = 12.0f;
float tess_max_level = 8.0f;
//The other way to draw the same water: no mesh, the particles splatted into a
//depth buffer which is then smoothed and shaded. See shaders.cpp. Off by
//default because the marching-cubes surface is what this program has always
//drawn; -S or the S key turns it on instead of the particle dots.
bool screenspace_on = false;
//Splat radius as a multiple of the particle radius. Particles sit about a
//diameter apart, so spheres of exactly one radius leave gaps between them and
//the smoothing pass has holes to fill; half again closes them without turning
//the spray into blobs.
float screenspace_radius = 1.5f;
//How many times the depth buffer is blurred. More is smoother and slower. Two
//passes of a filter that is itself as wide as a splat is plenty; it was four
//while the filter was capped far below that and had to be applied repeatedly to
//reach anywhere.
int screenspace_smoothing = 2;
bool showlight = false;
//The ground and the box the cells span. The ground was built into a display
//list in 2005 and then never drawn - the call was commented out - which is why
//the bowl has only ever been visible as the shape the water settles into. The
//box is the extent of the cell grid and nothing else: no particle is ever
//tested against it. It used to be drawn as a solid triangle strip through the
//eight corners, which comes out as three filled walls that the water appears
//to rest on; it is a wireframe now, which is what a reference box should be.
bool showground = true;
bool showbox = true;
//Which scene is running. Set from -X/--scene before the window exists.
const Scene * scene = 0;
float groundlevel ( float a, float b );
Vektor groundlevelnormal ( float a, float b );
void write_24bitbmp ( const char * filename, unsigned char * pixelinfo );
void read_24bitbmp ( const char * filename, unsigned char ** pixelinfo, unsigned long & width, unsigned long & height );

//What is on screen: one consistent view of the fluid, taken by the solver
//thread between two of its steps and handed over here. The viewer used to own
//these arrays and fill them itself from inside the display callback; now it
//only reads the newest one the solver has published. See simthread.h.
const FluidSnapshot * shown = 0;

//float g_normal_length = 1;

double rot[16]=
{
	1.0, 0.0, 0.0, 0.0,  // Rotationsmatrix fuer Ansicht
	0.0, 1.0, 0.0, 0.0,
	0.0, 0.0, 1.0, 0.0,
	0.0, 0.0, 0.0, 1.0
};
#define KUGELLIST 3
#define GROUND 5
#define GROUND_NORMALS 6
#define BOX 4
#define TRICHTER 7


#define MINX -200
#define MAXX 200
#define MINY -200
#define MAXY 200
#define MINZ -200
#define MAXZ 200

// The one array a scene lives in holds the water and the control particles that
// shape it both, and -n only asks for the water. Until 2026 the array was a
// constant 150000 that -n could be set equal to; the scene was then built off
// the end of it, which is a heap overrun and a segfault - -n 150000 laid down
// 150186 particles of water and then the funnel's control particles past them.
//
// Neither half of the count can be worked out in advance: the water is placed
// by subdividing tetrahedra until each piece holds one particle, so it lands
// near -n rather than on it, and how many control particles a scene needs
// depends on its geometry. So build the scene into a Fluid that only counts,
// and allocate the real one at what that says, plus a little for the solver's
// own use. It costs one extra layout pass, which is milliseconds.
//
// FLUIDSIM_MAX_PARTICLES overrides the size outright, for a run that wants room
// to spare. It is 120 bytes a particle.
static unsigned long room_for ( const Scene * sc, unsigned long wanted )
{
	const char * s = getenv ( "FLUIDSIM_MAX_PARTICLES" );
	unsigned long v = s ? strtoul ( s, 0, 10 ) : 0;
	if ( v )
		return v;
	Fluid counter ( 2, 0.5f, 0.0f, 0.0f, -9.81f, MINX, MINY, MINZ, MAXX, MAXY, MAXZ );
	counter.set_count_only ( true );
	sc->build ( counter, wanted );
	return counter.particlecount() + 64;
}

//Built in main() once the command line has been read, so that -n can size it.
static Fluid * fluid = 0;

using namespace std;
unsigned char ** rawdata;
bicubic_bezier_surface boden;

static void usage ( const char * argv0 )
{
	cout <<
	"usage: " << argv0 << " [options]\n"
	"\n"
	"  -X, --scene NAME    which scene to run (default fountain). --scene list\n"
	"                      names them all.\n"
	"  -n, --particles N   how much liquid, in particles. The default is what the\n"
	"                      scene is written for; this scales it.\n"
	"  -t, --threads N     OpenMP threads for the solver (default: as many as\n"
	"                      the machine has)\n"
	"  -r, --run           start running instead of paused\n"
	"  -s, --surface       start with the marching-cubes surface on\n"
	"  -F, --flat          do not tessellate the surface (PN triangles are on by\n"
	"                      default wherever the driver has OpenGL 4)\n"
	"  -S, --screenspace   draw the water with screen-space fluid rendering\n"
	"                      instead of particle dots - no mesh at all\n"
	"      --slope S       terrain scene: how steeply the channel is tilted\n"
	"                      (default 0.12)\n"
	"      --waves W       terrain scene: how deep the ripples along it are\n"
	"                      (default 1.8). There are dips for water to stand in\n"
	"                      only while W*0.126 > S; more W is deeper puddles and\n"
	"                      less running water.\n"
	"  -H, --sim-rate N    hold the solver to N steps per second (default: as\n"
	"                      many as it can take). The solver runs in a thread of\n"
	"                      its own and is not tied to the frame rate either way;\n"
	"                      at the default timestep, -H 250 is real time.\n"
	"  -h, --help          this\n"
	"\n"
	"The simulation is sized from -n, with room over for the scene's own control\n"
	"particles. FLUIDSIM_MAX_PARTICLES overrides that size outright, for a run\n"
	"that wants room to spare; it costs memory whether you fill it or not.\n"
	"Keys are listed in README.md. Examples:\n"
	"\n"
	"    " << argv0 << " -X funnel -n 150000 -r\n"
	"    " << argv0 << " -X terrain -n 40000 -r -s\n"
	"\n"
	"scenes:\n";
	scene_list();
}

int main ( int argc,char** argv )
{
	int i,j;

	unsigned long wanted_particles = 0;   //0: whatever the scene is written for
	int wanted_threads = 0;
	float wanted_sim_rate = 0.0f;   //0: as fast as the solver can go
	float wanted_slope = 0.0f, wanted_waves = 0.0f;   //0: the scene's own
	scene = scene_default();
	for ( i=1; i<argc; ++i )
	{
		string a = argv[i];
		if ( ( a=="-X" || a=="--scene" ) && i+1<argc )
		{
			string want = argv[++i];
			if ( want == "list" )
			{
				cout << "scenes:" << endl;
				scene_list();
				return 0;
			}
			scene = scene_by_name ( want.c_str() );
			if ( !scene )
			{
				cout << "no scene called " << want << ". There is:" << endl;
				scene_list();
				return 1;
			}
		}
		else if ( a=="--slope" && i+1<argc )
			wanted_slope = ( float ) atof ( argv[++i] );
		else if ( a=="--waves" && i+1<argc )
			wanted_waves = ( float ) atof ( argv[++i] );
		else if ( ( a=="-n" || a=="--particles" ) && i+1<argc )
			wanted_particles = strtoul ( argv[++i],0,10 );
		else if ( ( a=="-t" || a=="--threads" ) && i+1<argc )
			wanted_threads = atoi ( argv[++i] );
		else if ( a=="-r" || a=="--run" )
			paused = false;
		else if ( a=="-s" || a=="--surface" )
			showcells = true;
		else if ( a=="-F" || a=="--flat" )
			tessellation_on = false;
		else if ( a=="-S" || a=="--screenspace" )
			screenspace_on = true;
		else if ( ( a=="-H" || a=="--sim-rate" ) && i+1<argc )
			wanted_sim_rate = ( float ) atof ( argv[++i] );
		else if ( a=="-h" || a=="--help" )
		{
			usage ( argv[0] );
			return 0;
		}
	}
	if ( wanted_particles < 1 )
		wanted_particles = scene->natural_particles;

	rawdata=new unsigned char*;
	*rawdata=new unsigned char[10000];
	unsigned long w,h;
	for ( i=0; i<100; ++i )
	{
		for ( j=0; j<100; ++j )
		{
			//(*rawdata)[100*i+j]=(int)sink(15*(i-50),15*(j-50))%256;
			( *rawdata ) [100*i+j]= ( 5*j ) %256;
		}
	}
	write_24bitbmp ( "hight100x100.sink.bmp", *rawdata );
	read_24bitbmp ( "hight100x100.sink.bmp", rawdata, w, h );

	boden = bicubic_bezier_surface ( *rawdata );
	boden.callf ( 100,100 );
	boden.interpolate();

	//Everything about the run that is not the solver comes from the scene; see
	//scenes.cpp. Until 2026 there was one of these hard-coded here with the
	//others commented out in place, so picking one meant an edit and a rebuild.
	terrain_set ( wanted_slope, wanted_waves );
	//...and only now is the scene's geometry settled enough to be measured, so
	//this is where the array it lives in can be sized. See room_for().
	fluid = new Fluid ( room_for ( scene, wanted_particles ),
	                    0.5f,
	                    0.0f, 0.0f, -9.81f,
	                    MINX, MINY, MINZ,
	                    MAXX, MAXY, MAXZ );
	fluid->set_height_function ( scene->height );
	fluid->set_height_function_normal ( scene->height_normal );
	fluid->set_ground_restitution ( scene->ground_restitution );
	scene->build ( *fluid, wanted_particles );
	distanz = scene->camera_distance;
	showground = scene->ground_shown;
	{
		//The scene's own starting view, as the same rotation matrix the mouse
		//builds up, so moving the mouse carries on from here.
		//The rotation has to tip the scene towards the viewer, which means world z
		//comes out of the screen upwards: the column for z is (0,sa,ca), not
		//(0,-sa,ca). With the sign the other way round the view is not tipped but
		//turned over - the funnel's widest ring came out at the bottom of the
		//window, its throat above it and the ground at z=0 above that.
		double a = scene->camera_pitch*M_PI/180.0, ca = cos ( a ), sa = sin ( a );
		double r[16] = { 1,0,0,0,  0,ca,-sa,0,  0,sa,ca,0,  0,0,0,1 };
		for ( int m=0;m<16;++m ) rot[m] = r[m];
	}
	cout << "scene " << scene->name << ": " << scene->summary << endl;
	cout << fluid->movingparticlecount() << " liquid particles, "
	     << fluid->boundaryparticlecount() << " control particles." << endl;
	if ( fluid->refusedcount() )
		cout << "The scene did not fit: " << fluid->refusedcount() << " particles of it were "
		     "dropped, out of " << fluid->particlecount() +fluid->refusedcount() << ". Set "
		     "FLUIDSIM_MAX_PARTICLES above " << fluid->maxparticlecount() << "." << endl;

	//The solver gets a thread of its own, and everything the command line said
	//about how it should run. From here on it is the only thread that advances
	//the fluid; the viewer reads what it publishes, and otherwise only the
	//handful of things about it that never change once the scene is built.
	sim_set_dt ( speed );
	sim_set_paused ( paused );
	sim_set_rate_limit ( wanted_sim_rate );
	sim_start ( *fluid, wanted_threads );

	glutInit ( &argc, argv );
	initMain();
#ifdef WIN32
	SetupATIExtensions();
#endif
	glutTimerFunc ( 10,zeitgeber, 10 );
	glutMainLoop ();

	sim_stop();
	exit ( 554 );

	return 1;
}
void read_24bitbmp ( const char * filename, unsigned char ** pixelinfo, unsigned long & width, unsigned long & height )
{
	unsigned char buffer[3];
	ifstream stream ( filename, ios::in | ios::binary );
	if ( stream )
	{
		stream.read ( ( char * ) &mybmpfileheader.bfType, 2 );
		stream.read ( ( char * ) &mybmpfileheader.bfSize, 4 );
		stream.read ( ( char * ) &mybmpfileheader.bfReserved1, 2 );
		stream.read ( ( char * ) &mybmpfileheader.bfReserved2, 2 );
		stream.read ( ( char * ) &mybmpfileheader.bfOffBits, 4 );

		stream.read ( ( char * ) &mybmpinfoheader.biSize, 4 );
		stream.read ( ( char * ) &mybmpinfoheader.biWidth, 4 );
		width=abs ( mybmpinfoheader.biWidth );
		stream.read ( ( char * ) &mybmpinfoheader.biHeight, 4 );
		height=abs ( mybmpinfoheader.biHeight );
		stream.read ( ( char * ) &mybmpinfoheader.biPlanes, 2 );
		stream.read ( ( char * ) &mybmpinfoheader.biBitCount, 2 );
		stream.read ( ( char * ) &mybmpinfoheader.biCompression, 4 );
		stream.read ( ( char * ) &mybmpinfoheader.biSizeImage, 4 );
		stream.read ( ( char * ) &mybmpinfoheader.biXPelsPerMeter, 4 );
		stream.read ( ( char * ) &mybmpinfoheader.biYPelsPerMeter, 4 );
		stream.read ( ( char * ) &mybmpinfoheader.biClrUsed, 4 );
		stream.read ( ( char * ) &mybmpinfoheader.biClrImportant, 4 );
		if ( *pixelinfo != NULL )
			delete [] *pixelinfo;
		*pixelinfo = new unsigned char[height*width];
		for ( unsigned int i=0;i<height;i++ )
		{
			for ( unsigned int j=0;j<width;j++ )
			{
				stream.read ( ( char * ) buffer, 3 );
				( *pixelinfo ) [i*width+j]= ( unsigned char ) ( ( ( int ) buffer[0]+buffer[1]+buffer[2] ) /3 );
			}
		}
	}
	else
	{
		exit ( 76456 );
	}
	stream.close();
}
void write_24bitbmp ( const char * filename, unsigned char * pixelinfo )
{
	memset ( &mybmpfileheader,0,sizeof ( mybmpfileheader ) );
	mybmpfileheader.bfType=19778%255;//TODO: was 19778;
	mybmpfileheader.bfSize=sizeof ( BMP_FILEHEADER ) +sizeof ( BMP_INFOHEADER ) +100*100*24/8-2;
	mybmpfileheader.bfReserved1=mybmpfileheader.bfReserved2=0;
	mybmpfileheader.bfOffBits=sizeof ( BMP_FILEHEADER ) +sizeof ( BMP_INFOHEADER )-2;

	memset ( &mybmpinfoheader,0,sizeof ( mybmpinfoheader ) );
	mybmpinfoheader.biSize=sizeof ( BMP_INFOHEADER );
	mybmpinfoheader.biWidth=100;
	mybmpinfoheader.biHeight=100;
	mybmpinfoheader.biPlanes=1;
	mybmpinfoheader.biBitCount=24;
	mybmpinfoheader.biCompression=0;
	mybmpinfoheader.biSizeImage=0;
	mybmpinfoheader.biXPelsPerMeter=7874;
	mybmpinfoheader.biYPelsPerMeter=7874;
	mybmpinfoheader.biClrUsed=0;
	mybmpinfoheader.biClrImportant=0;

	ofstream stream2 ( filename, ofstream::binary );
	stream2.write ( ( char * ) &mybmpfileheader.bfType,2 );
	stream2.write ( ( char * ) &mybmpfileheader.bfSize,4 );
	stream2.write ( ( char * ) &mybmpfileheader.bfReserved1,2 );
	stream2.write ( ( char * ) &mybmpfileheader.bfReserved2,2 );
	stream2.write ( ( char * ) &mybmpfileheader.bfOffBits,4 );
	stream2.write ( ( char * ) &mybmpinfoheader.biSize,4 );
	stream2.write ( ( char * ) &mybmpinfoheader.biWidth,4 );
	stream2.write ( ( char * ) &mybmpinfoheader.biHeight,4 );
	stream2.write ( ( char * ) &mybmpinfoheader.biPlanes,2 );
	stream2.write ( ( char * ) &mybmpinfoheader.biBitCount,2 );
	stream2.write ( ( char * ) &mybmpinfoheader.biCompression,4 );
	stream2.write ( ( char * ) &mybmpinfoheader.biSizeImage,4 );
	stream2.write ( ( char * ) &mybmpinfoheader.biXPelsPerMeter,4 );
	stream2.write ( ( char * ) &mybmpinfoheader.biYPelsPerMeter,4 );
	stream2.write ( ( char * ) &mybmpinfoheader.biClrUsed,4 );
	stream2.write ( ( char * ) &mybmpinfoheader.biClrImportant,4 );

	//stream2.write((char *)&mybmpfileheader,sizeof(mybmpfileheader));;
	//stream2.write((char *)&mybmpinfoheader,sizeof(mybmpinfoheader));;
	for ( int i=0;i<100;i++ )
	{
		for ( int j=0;j<100;j++ )
		{
			stream2.write ( ( char * ) &pixelinfo[i*mybmpinfoheader.biWidth+j],1 );
			stream2.write ( ( char * ) &pixelinfo[i*mybmpinfoheader.biWidth+j],1 );
			stream2.write ( ( char * ) &pixelinfo[i*mybmpinfoheader.biWidth+j],1 );
		}
	}
	stream2.close();
}
//The height field read out of hight100x100.sink.bmp through a bicubic Bezier
//patch. No scene uses it - the ones that ship are analytic, in scenes.cpp -
//but the bitmap path is the only way this program can take a terrain from a
//file, so it stays.
float groundlevel ( float a, float b )
{
	a= ( a-MINX ) / ( MAXX-MINX );
	b= ( b-MINY ) / ( MAXY-MINY );
	//if((a>0)&&(a<1)&&(b>0)&&(b<1))
	return boden ( a,b ) /10;
	//return 128;
}
Vektor groundlevelnormal ( float a, float b )
{
	a= ( a-MINX ) / ( MAXX-MINX );
	b= ( b-MINY ) / ( MAXY-MINY );
	//if((a>0)&&(a<1)&&(b>0)&&(b<1))
	return boden.normal ( a,b );
	//return 128;
}

void initCallLists ( void )
{
	float k=M_PI/6;
	float i,j;
	float v[3];
	//The funnel, from the same ring numbers scenes.cpp hands the solver, so
	//what is drawn is what the water hits. In 2005 this block was commented
	//out while glCallList(TRICHTER) was not, so it called a list that had never
	//been defined and drew nothing.
	glNewList ( TRICHTER, GL_COMPILE );
	if ( scene->draw_funnel )
	{
		//As a wireframe, for the same reason the box is one: filled, it is an
		//opaque cone with the water inside it.
		float a[3],b[3];
		glDisable ( GL_LIGHTING );
		glBegin ( GL_LINES );
		glColor4f ( 0.55f,0.55f,0.62f,1 );
		for ( int jj=1;jj<6;jj++ )
		{
			for ( int ii=0;ii<10;ii++ )
			{
				funnel_ring ( jj,ii,  a );
				funnel_ring ( jj,ii+1,b );
				glVertex3fv ( a ); glVertex3fv ( b );
				if ( jj>1 )
				{
					funnel_ring ( jj-1,ii,b );
					glVertex3fv ( a ); glVertex3fv ( b );
				}
			}
		}
		glEnd();
	}
	glEndList();
	glNewList ( KUGELLIST, GL_COMPILE );
	glMaterialf ( GL_FRONT,GL_SHININESS,50 );
	for ( i=0;i<M_PI-k;i+=k )
	{
		glBegin ( GL_TRIANGLE_STRIP );
		for ( j=2*M_PI;j>=0;j-=k )
		{
			v[0]=sin ( i );
			v[1]=1-sin ( i );
			v[2]=0.5+sin ( j ) /2;
			glMaterialfv ( GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,v );
			glColor3fv ( v );
			v[0]=sin ( i ) *cos ( j );
			v[1]=cos ( i );
			v[2]=sin ( i ) *sin ( j );
			glNormal3fv ( v );
			glVertex3fv ( v );

			v[0]=sin ( i+k );
			v[1]=1-sin ( i+k );
			v[2]=0.5+sin ( j ) /2;
			glMaterialfv ( GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,v );
			glColor3fv ( v );
			glNormal3f ( sin ( i+k ) *cos ( j ),cos ( i+k ),sin ( i+k ) *sin ( j ) );
			glVertex3f ( sin ( i+k ) *cos ( j ),cos ( i+k ),sin ( i+k ) *sin ( j ) );
		}
		glEnd();
	}
	glEndList();
	//The ground, over whatever patch of it the scene says is worth looking at.
	//It was drawn over a fixed 20% of the domain before, and the call to draw
	//it was commented out, so the terrain has never actually been on screen.
	//The colour is baked in rather than lit, because lighting is off by
	//default and an unlit height field is a flat silhouette.
	Vektor tmp;
	glNewList ( GROUND, GL_COMPILE );
	{
		const float lx = 0.4360f, ly = 0.3270f, lz = 0.8384f;  //a fixed sun
		const float dx = ( scene->ground_x1 - scene->ground_x0 ) / scene->ground_nx;
		const float dy = ( scene->ground_y1 - scene->ground_y0 ) / scene->ground_ny;
		float zlo = 1e30f, zhi = -1e30f;
		for ( int a=0;a<=scene->ground_nx;++a )
			for ( int b=0;b<=scene->ground_ny;++b )
			{
				float z = fluid->height_function ( scene->ground_x0+a*dx, scene->ground_y0+b*dy );
				if ( z<zlo ) zlo=z;
				if ( z>zhi ) zhi=z;
			}
		if ( zhi-zlo < 1e-3f )
			zhi = zlo + 1;
		for ( int a=0;a<scene->ground_nx;++a )
		{
			glBegin ( GL_TRIANGLE_STRIP );
			for ( int b=0;b<=scene->ground_ny;++b )
			{
				for ( int e=0;e<2;++e )
				{
					float x = scene->ground_x0 + ( a+e ) *dx;
					float y = scene->ground_y0 + b*dy;
					float z = fluid->height_function ( x,y );
					Vektor n = fluid->height_function_normal ( x,y ).normed();
					float lit = 0.45f + 0.55f* ( n.x() *lx + n.y() *ly + n.z() *lz );
					if ( lit < 0 ) lit = 0;
					float t = ( z-zlo ) / ( zhi-zlo );
					glNormal3f ( n.x(),n.y(),n.z() );
					glColor4f ( ( 0.30f+0.45f*t ) *lit, ( 0.34f+0.30f*t ) *lit, ( 0.26f+0.28f*t ) *lit, 1 );
					glVertex3f ( x,y,z );
				}
			}
			glEnd();
		}
	}
	glEndList();

	glNewList ( GROUND_NORMALS, GL_COMPILE );
	{
		const float dx = ( scene->ground_x1 - scene->ground_x0 ) / scene->ground_nx;
		const float dy = ( scene->ground_y1 - scene->ground_y0 ) / scene->ground_ny;
		glBegin ( GL_LINES );
		glColor4f ( 0,1,0,1 );
		for ( int a=0;a<=scene->ground_nx;a+=4 )
			for ( int b=0;b<=scene->ground_ny;b+=4 )
			{
				float x = scene->ground_x0 + a*dx, y = scene->ground_y0 + b*dy;
				Vektor n = fluid->height_function_normal ( x,y ).norm ( 0.03f* ( scene->ground_x1-scene->ground_x0 ) );
				glVertex3f ( x,y,fluid->height_function ( x,y ) );
				glVertex3f ( x+n.x(),y+n.y(),fluid->height_function ( x,y ) +n.z() );
			}
		glEnd();
	}
	glEndList();

	//The extent of the cell grid, as twelve lines. It used to be a triangle
	//strip through the eight corners, which is three filled walls: opaque,
	//pastel, and close enough behind the water to look like a floor the
	//particles are resting on. Nothing in the solver ever tests a particle
	//against it - it is a ruler, not a wall.
	glNewList ( BOX, GL_COMPILE );
	{
		Vektor lo = fluid->minxyz(), hi = fluid->maxxyz();
		const float xs[2] = { lo.x(),hi.x() }, ys[2] = { lo.y(),hi.y() }, zs[2] = { lo.z(),hi.z() };
		glDisable ( GL_LIGHTING );
		glColor4f ( 0.35f,0.35f,0.40f,1 );
		glBegin ( GL_LINES );
		for ( int b=0;b<2;++b )
			for ( int c=0;c<2;++c )
			{
				glVertex3f ( xs[0],ys[b],zs[c] ); glVertex3f ( xs[1],ys[b],zs[c] );
				glVertex3f ( xs[b],ys[0],zs[c] ); glVertex3f ( xs[b],ys[1],zs[c] );
				glVertex3f ( xs[b],ys[c],zs[0] ); glVertex3f ( xs[b],ys[c],zs[1] );
			}
		glEnd();
	}
	glEndList();
}
void PrintVolumeOfClosedSurface()
{
	unsigned long i;
	const float * va;
	const unsigned int * ia;
	float volume,dv;
	if ( !showcells || !shown || !shown->has_surface )
	{
		cout << "Without Surfacevisualisation, surface is not calculated and thus makes no sence beeing asked its volume." << endl;
		return;
	}
	va = shown->vertices;
	ia = shown->indices;
	volume = 0;
	for ( i=0;i<shown->trianglecount;i++ )
	{
		dv =		va[3*ia[3*i+0]+0] * va[3*ia[3*i+1]+1] * va[3*ia[3*i+2]+2] +
		      va[3*ia[3*i+1]+0] * va[3*ia[3*i+2]+1] * va[3*ia[3*i+0]+2] +
		      va[3*ia[3*i+2]+0] * va[3*ia[3*i+0]+1] * va[3*ia[3*i+1]+2] -
		      va[3*ia[3*i+2]+0] * va[3*ia[3*i+1]+1] * va[3*ia[3*i+0]+2] -
		      va[3*ia[3*i+1]+0] * va[3*ia[3*i+0]+1] * va[3*ia[3*i+2]+2] -
		      va[3*ia[3*i+0]+0] * va[3*ia[3*i+2]+1] * va[3*ia[3*i+1]+2];
		volume -=dv;
	}
	cout << "Actual Surface encloses " << ( volume * 1000 ) << " liters" << endl;
}
void kbf ( unsigned char key,int x, int y )
{
	switch ( key )
	{
		case '+' :
			distanz*=0.8;
			break;
		case '-' :
			distanz/=0.8;
			break;
		case '*' :
			speed*=1.1;
			sim_set_dt ( speed );
			cout << "Speed = " << speed << endl;
			break;
		case '/' :
			speed/=1.1;
			sim_set_dt ( speed );
			cout << "Speed = " << speed << endl;
			break;
		case 9 ://tab-key
			glPolygonMode ( GL_FRONT,GL_NONE );
			glPolygonMode ( GL_BACK,GL_FILL );
			break;
		case 'n' :
			shownormals = !shownormals;
			break;
		case 'v' :
			showparticles = !showparticles;
			break;
		case 'y' :
			PrintVolumeOfClosedSurface();
			break;
		case 'c' :
			showcells = !showcells;
			if ( showcells && screenspace_on && screenspace_available() )
				cout << "the marching-cubes surface is on, but screen-space is drawing the water; S to go back to it" << endl;
			break;
		case 'b' :
			showground = !showground;
			break;
		case 'B' :
			showbox = !showbox;
			break;
		case 'e' :
			glPolygonMode ( GL_FRONT,GL_FILL );
			glPolygonMode ( GL_BACK,GL_LINE );
			break;
		case 'f' :
			glShadeModel ( GL_FLAT );
			break;
		case 'g' :
			glShadeModel ( GL_SMOOTH );
			break;
		case 'l' :
			showlight = !showlight;
			break;
		case 't' :
			tessellation_on = !tessellation_on;
			if ( tessellation_on && !surface_tessellation_available() )
				cout << "no tessellation here: " << surface_shaders_error() << endl;
			else
				cout << "PN triangle tessellation " << ( tessellation_on ? "on" : "off" ) << endl;
			break;
		case 'S' :
			screenspace_on = !screenspace_on;
			if ( screenspace_on && !screenspace_available() )
			{
				cout << "no screen-space fluid here: " << screenspace_error() << endl;
				screenspace_on = false;
			}
			else
				cout << "screen-space fluid " << ( screenspace_on ? "on" : "off" )
				     << ( screenspace_on && showcells ? " (the marching-cubes surface is off while it is)" : "" ) << endl;
			break;
		case 'T' :
			tess_max_level *= 2.0f;
			if ( tess_max_level > 16.0f )
				tess_max_level = 2.0f;
			cout << "tessellation level up to " << tess_max_level << endl;
			break;
		case 'p' :
			paused = !paused;
			sim_set_paused ( paused );
			break;
			/*    case 'n' :
			        fluid->deb_colliderecording() ? fluid->deb_hidecollide() : fluid->deb_showcollide();
					ReshapeMain(glutGet(GLUT_WINDOW_WIDTH),glutGet(GLUT_WINDOW_HEIGHT));
			        break;*/
		case 'q' :
			glPolygonMode ( GL_FRONT_AND_BACK,GL_POINT );
			break;
//#ifdef debug
//    case 't' :
//        timechecker(10);
//        break;
//#endif
		case 'w' :
			glPolygonMode ( GL_FRONT_AND_BACK,GL_LINE );
			break;
		case 'o' :
		{
			if ( !shown || !shown->has_surface )
				break;
			const float * va = shown->vertices;
			const float * na = shown->normals;
			cout << shown->vn_alloc << endl;
			for ( unsigned long i=0;i<shown->vn_count;++i )
			{
				float l;
				l = sqrt ( na[3*i+1]*na[3*i+1]+na[3*i]*na[3*i]+na[3*i+2]*na[3*i+2] );
				printf ( "V %f %f %f %f %f %f 0.5 0.5\n", va[3*i+1], va[3*i], va[3*i+2], na[3*i+1]/l, na[3*i]/l, na[3*i+2]/l );
			}
			for ( unsigned long i=0;i<shown->trianglecount*3;++i )
			{
				printf ( "I %u\n", shown->indices[i] );
			}
			break;
		}
		case 27  :
//#ifdef debug
//     timechecker(10);
//#endif
			exit ( 4 );
			break;
	}
}

void kbf2 ( int key,int x, int y )
{
	switch ( key )
	{
			/*    case GLUT_KEY_F1    :
			        fluid->deb_print_relevant_grid();
			        break;*/
		case GLUT_KEY_F2    :
			initMain();
			printf ( "WinMode\n" );
			break;
		case GLUT_KEY_F3    :
		case GLUT_KEY_F4    :
		case GLUT_KEY_F5    :
		case GLUT_KEY_F6    :
		case GLUT_KEY_F7    :
		case GLUT_KEY_F8    :
		case GLUT_KEY_F9    :
		case GLUT_KEY_F11   :
		case GLUT_KEY_F12   :
			break;
//	#ifdef debug
//    case GLUT_KEY_LEFT  :
//        fluid->deb_showcollideprev();
////        cout << "Showing collisions with particle " << fluid->deb_get_bumpingparticle() << endl;
//        break;
//    case GLUT_KEY_RIGHT :
//        fluid->deb_showcollidenext();
////        cout << "Showing collisions with particle " << fluid->deb_get_bumpingparticle() << endl;
//        break;
//	#else
//    case GLUT_KEY_LEFT  :
//        g_normal_length/=1.05f;
//        cout << "g_normal_length is set to " << g_normal_length << endl;
////        cout << "Showing collisions with particle " << fluid->deb_get_bumpingparticle() << endl;
//        break;
//    case GLUT_KEY_RIGHT :
//        g_normal_length*=1.05;
//        cout << "g_normal_length is set to " << g_normal_length << endl;
////        cout << "Showing collisions with particle " << fluid->deb_get_bumpingparticle() << endl;
//        break;
			//#endif
			/*    case GLUT_KEY_UP    :
			        if (deb_checkbit < 512)
			            deb_checkbit <<= 1;
			        cout << "checkbit=" << deb_checkbit << endl;
			        break;
			    case GLUT_KEY_DOWN  :
			        if (deb_checkbit > 1)
			            deb_checkbit >>= 1;
			        cout << "checkbit=" << deb_checkbit << endl;
			        break;
			#endif*/
		default  :
			break;
	}
}


void kbuf ( unsigned char key,int x, int y )
{
	switch ( key )
	{
		case 'n' :
			break;
	}

}

void kbuf2 ( int key,int x, int y )
{
	switch ( key )
	{
		case GLUT_KEY_DOWN  :
			break;
		default  :
			break;
	}
}
void DisplayMain ( void )
{
	//The solver runs on its own thread (simthread.cpp) and this takes the newest
	//state it has published, then asks it for the next one. Nothing here steps
	//the simulation any more, so the vertical retrace this function ends on no
	//longer decides how fast the fluid moves.
	//Screen-space fluid and the marching-cubes surface are two ways of drawing
	//the same water, not two things to draw, so the one that is on wins and the
	//other is not built at all. Asking the solver for no surface is what makes
	//-S cheap: it is the mesh rebuild, once per frame, that screen-space exists
	//to avoid.
	const bool screenspace = screenspace_on && screenspace_available();
	const bool surface = showcells && !screenspace;
	shown = sim_acquire_snapshot ( surface );

	if ( showlight )
		glEnable ( GL_LIGHTING );
	else
		glDisable ( GL_LIGHTING );


	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glLoadIdentity ();//	gluPerspective (90,1.4,0.01,1000);
	glColor3f ( 1.0,1.0,1.0 );

	glTranslatef ( 0,0,-distanz );
	glMultMatrixd ( rot );


//					float texcoord[6]={0,0,  1,0,  1,1}; //u-v-Koordinate fr 1., 2. & 3.Vertice

//					float normal[9]={0,0,1,  0,0,1,  0,0,1}; //x-y-z-Koordinate fr 1.-3.verticesnormale

	/*
						for(i=0;i<vertexcount;indices[i]=i++);*/

//					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
//					glEnableClientState(GL_NORMAL_ARRAY);
	//definieren von VA-Format
//					glNormalPointer(3, GL_FLOAT, 0, normal);
//					glTexCoordPointer(2, GL_FLOAT, 0,texcoord);
	//					glDrawElements(GL_TRIANGLES,vertexcount / 3, GL_UNSIGNED_INT, indices);
	if ( surface && shown && shown->has_surface )
	{
		trianglecount = shown->trianglecount;
		glEnableClientState ( GL_VERTEX_ARRAY );
		glEnableClientState ( GL_NORMAL_ARRAY );
		glVertexPointer ( 3, GL_FLOAT, 0, shown->vertices );
		glNormalPointer ( GL_FLOAT, 0, shown->normals );
//#ifdef debug
//		glEnableClientState(GL_COLOR_ARRAY );
//		glColorPointer(3, GL_FLOAT, 3*sizeof(GLfloat), *normal_array_pointer);
//#else
		glDisableClientState ( GL_COLOR_ARRAY );
//#endif

		glColor4f ( 0.66f,0.66f,1.0f,0.4f );

		//The same triangles either way. With tessellation they go in as patches
		//and come back out of the evaluation shader as a curved surface; without
		//it they are drawn flat, exactly as before.
		const bool tess = tessellation_on && surface_tessellation_available();
		if ( tess )
		{
			surface_tessellation_bind ( tess_pixels_per_segment, tess_max_level );
			glDrawElements ( GL_PATCHES, 3*trianglecount, GL_UNSIGNED_INT, shown->indices );
			surface_tessellation_unbind();
		}
		else
		{
			glDrawElements ( GL_TRIANGLES, 3*trianglecount, GL_UNSIGNED_INT, shown->indices );
		}
		//glDrawArrays(GL_POINTS,0,vertexcount);
		/*glBegin(GL_LINES);
		for (i=0; i < vertexcount-3; i+=3)
			glVertex3f (va[i],va[i+1],va[i+2]);
		glEnd();*/
		if ( shownormals )
		{
			glDisable ( GL_LIGHTING );
			const float *tmp_norm=shown->normals;
			const float *tmp_vert=shown->vertices;
			glBegin ( GL_LINES );
			glColor4f ( 1.0f,0.86f,0.86f,1.0f );
			for ( unsigned long i=0; i<shown->vn_count; i++ )
			{
				glVertex3f ( tmp_vert[3*i+0],tmp_vert[3*i+1],tmp_vert[3*i+2] );
				glVertex3f ( tmp_vert[3*i+0]+tmp_norm[3*i+0]*1.01f,tmp_vert[3*i+1]+tmp_norm[3*i+1]*1.01f,tmp_vert[3*i+2]+tmp_norm[3*i+2]*1.01f );
			}
			glEnd();
			if ( showlight )
				glEnable ( GL_LIGHTING );
			else
				glDisable ( GL_LIGHTING );
		}
	}
	//Screen-space fluid takes the place of the particle dots: it is the same
	//particles, drawn as a surface instead of as points. Drawn further down,
	//after the rest of the scene, because it is transparent.
	if ( !screenspace && showparticles && shown )
	{
		/* set up the array data */
		glVertexPointer ( 3, GL_FLOAT, 3*sizeof ( GLfloat ), shown->particles );
		//glColorPointer(3, GL_FLOAT, 3*sizeof(GLfloat), *particle_colors_pointer);

		/* enable vertex arrays */
		glDisableClientState ( GL_NORMAL_ARRAY );
		glEnableClientState ( GL_VERTEX_ARRAY );
		//glEnableClientState( GL_COLOR_ARRAY );

		/* draw a polygon using the arrays sequentially */
		//Water-coloured rather than white: on a terrain scene a thin sheet of
		//white dots on pale ground is invisible, and the thin sheet is the
		//interesting part. The control particles - the barriers, teleports,
		//shifts and set-speed patches the scene is built out of - come after
		//the moving ones in the same array and are drawn dim, so that a scene
		//with two thousand of them shows its plumbing without burying the water
		//in it.
		glColor3f ( 0.55f,0.80f,1.0f );
		glDrawArrays ( GL_POINTS,0,shown->movingparticlecount );
		if ( shown->particlecount > shown->movingparticlecount )
		{
			glColor3f ( 0.22f,0.22f,0.26f );
			glDrawArrays ( GL_POINTS,shown->movingparticlecount,
			               shown->particlecount - shown->movingparticlecount );
		}
		glColor3f ( 1.0f,1.0f,1.0f );
	}
	if ( showground )
		glCallList ( GROUND );
	glCallList ( TRICHTER );
	if ( shownormals )
		glCallList ( GROUND_NORMALS );
	glDisable ( GL_LIGHTING );
	if ( showbox )
		glCallList ( BOX );
	glBegin ( GL_LINES );
	glColor4f ( 1,0,0,1 );
	glVertex3f ( 0,0,0 );
	glVertex3f ( 4*fluid->particleradius(),0,0 );
	glColor4f ( 0,1,0,1 );
	glVertex3f ( 0,0,0 );
	glVertex3f ( 0,4*fluid->particleradius(),0 );
	glColor4f ( 0,0,1,1 );
	glVertex3f ( 0,0,0 );
	glVertex3f ( 0,0,4*fluid->particleradius() );
	glEnd();

	//Last, because the water is transparent and has to be composited over
	//whatever is behind it. Drawn any earlier it blends against the cleared
	//background instead of against the scene, and comes out dark.
	if ( screenspace && shown )
	{
		screenspace_render ( shown->particles, shown->movingparticlecount,
		                     fluid->particleradius() *screenspace_radius, screenspace_smoothing,
		                     ( float ) distanz, showlight );
	}

	glutSwapBuffers ();
	framecount++;
}

void mousemove ( int x, int y )
{
	static int xalt;
	static int yalt;
	int dx,dy;
	xt=x;
	yt=y;
	dx=xalt-x;
	dy=yalt-y;

	glLoadIdentity ();
	glRotated ( -dx, 0.0 , 0.01, 0.0 );
	glRotated ( -dy, 0.01, 0.0 , 0.0 );
	glMultMatrixd ( rot );
	glGetDoublev ( GL_MODELVIEW_MATRIX,rot );
	xalt=x;
	yalt=y;
}

void zeitgeber ( int value )
{
	static int talt = 0;
	int t = glutGet ( GLUT_ELAPSED_TIME );
	talt = t;
	static int sekunden = 0;
	if ( ( sekunden+1 ) *1000 <= t )
	{

		cout << sekunden << ". Sekunde: " << framecount << "frames, "
		     << sim_steps_since_last_call() << " steps, "
		     << fluid->particlecount() << " particles, " << fluid->movingparticlecount()
		     << " moving particles, " << trianglecount << "triangles." << endl;
		trianglecount=0;
		sekunden ++;
		framecount = 0;
	}
	glutTimerFunc ( 200,zeitgeber, 1000 );
}

void ReshapeMain ( GLint width, GLint height )
{
	/*	if(fluid->deb_colliderecording()) {
			GLsizei tmp1 = glutGet(GLUT_WINDOW_WIDTH);
			GLsizei tmp2 = glutGet(GLUT_WINDOW_HEIGHT);
			glutSetWindow(winIdMain);
			glViewport(0, 0, tmp1, tmp2);
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
		    gluPerspective(90.0, (float)tmp1 / tmp2, distanz - fluid->particleradius(),distanz + fluid->particleradius());
			glMatrixMode(GL_MODELVIEW);
		} else {*/
	glViewport ( 0, 0, width, height );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective ( 90.0, ( float ) width / height, 0.01, 1000 );
	glMatrixMode ( GL_MODELVIEW );
//	}
}
int DrawTextXY ( float x,float y,float scale,char *s )
{
	unsigned int i;
	int l=0;
	glPushMatrix();
	glTranslatef ( x,y,0.3 );
	glScalef ( scale,scale,scale );
	for ( i=0;i<strlen ( s );i++ )
	{
		l+=glutStrokeWidth ( GLUT_STROKE_ROMAN, s[i] );
		glutStrokeCharacter ( GLUT_STROKE_ROMAN,s[i] );
	}
	glPopMatrix();
	return l;
}
void initMain()
{
	if ( glutGameModeGet ( GLUT_GAME_MODE_ACTIVE ) )
		glutLeaveGameMode();
	glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
	glutInitWindowPosition ( 0, 0 );
	glutInitWindowSize ( WIDTH, HEIGHT );
	glPolygonMode ( GL_FRONT,GL_FILL );
	glPolygonMode ( GL_BACK,GL_LINE );
	glClearDepth ( 0.0 );
	if ( glutGetWindow() )
		glutDestroyWindow ( winIdMain );
	winIdMain = glutCreateWindow ( ( char* ) TITLE );
	//Needs the context, so not before the window exists.
	surface_shaders_init();

	float mat_specular[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
	float mat_shininess = 50.0;
	float light0_position[4] = { -40.0, 20.0, 155.0, 0.0 };
	float light0_ambient[4] = { 1.0, 0.1, 0.0, 1.0 }; // Define some ambient light to add to the scene
	float light1_position[4] = { 30.0, -80.0, -55.0, 0.0 };
	float light1_ambient[4] = { 0.2, 0.3, 1.0, 1.0 }; // Define some ambient light to add to the scene


	glClearColor ( 0.0, 0.0, 0.0, 1.0 ); // Black Background
	glShadeModel ( GL_SMOOTH );       // Use Smooth shading ( This is the Default so we dont actually have to set it)

	glMaterialfv ( GL_FRONT, GL_SPECULAR, mat_specular );
	glMaterialfv ( GL_FRONT, GL_SHININESS, &mat_shininess );

	//glLightfv(GL_LIGHT0,GL_AMBIENT, light0_ambient);
	glLightfv ( GL_LIGHT0,GL_DIFFUSE, light0_ambient );
	glLightfv ( GL_LIGHT0, GL_POSITION, light0_position );
	// Light 1
	//glLightfv(GL_LIGHT1,GL_AMBIENT, light1_ambient);
	glLightfv ( GL_LIGHT1,GL_DIFFUSE, light1_ambient );
	glLightfv ( GL_LIGHT1, GL_POSITION, light1_position );

	glEnable ( GL_LIGHTING );
	glEnable ( GL_LIGHT0 );
	glEnable ( GL_LIGHT1 );
	glEnable ( GL_COLOR_MATERIAL );
	glEnable ( GL_NORMALIZE );
	glLightModelf ( GL_LIGHT_MODEL_LOCAL_VIEWER,1 );
	glLightModelf ( GL_LIGHT_MODEL_TWO_SIDE,0 );

	glDepthFunc ( GL_LEQUAL );
	glEnable ( GL_DEPTH_TEST );

	glEnable ( GL_CULL_FACE );
	glDepthFunc ( GL_LESS );
	glEnable ( GL_DEPTH_TEST );

	glEnable ( GL_BLEND );
	//	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glutDisplayFunc ( DisplayMain );
	glutReshapeFunc ( ReshapeMain );
	glutIdleFunc ( DisplayMain );

	glutIgnoreKeyRepeat ( 1 );
	glutPassiveMotionFunc ( mousemove );
	glutKeyboardFunc ( kbf );
	glutKeyboardUpFunc ( kbuf );
	glutSpecialFunc ( kbf2 );
	glutSpecialUpFunc ( kbuf2 );
	initCallLists();
	glPointSize ( 2.0 );
	//    glEnable(GL_POLYGON_SMOOTH);
	glCullFace ( GL_BACK );

	makeMenu();
}

void makeMenu()
{
	Menu = glutCreateMenu ( Menufunc );
	glutSetMenu ( Menu );
	glutAttachMenu ( GLUT_RIGHT_BUTTON );
	glutAddMenuEntry ( "dots",1 );
	glutAddMenuEntry ( "wire",2 );
	glutAddMenuEntry ( "filled",3 );
	glutAddMenuEntry ( "flat/smooth",4 );

	glutAddMenuEntry ( "exit",99 );
}

void Menufunc ( int value )
{
	int z[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	switch ( value )
	{
		case 1 :
			glPolygonMode ( GL_FRONT_AND_BACK,GL_POINT );
			break;
		case 2 :
			glPolygonMode ( GL_FRONT_AND_BACK,GL_LINE );
			break;
		case 3 :
			glPolygonMode ( GL_FRONT_AND_BACK,GL_FILL );
			break;
		case 4 :
			glGetIntegerv ( GL_SHADE_MODEL, z );
			z[0]==7425 ? glShadeModel ( GL_FLAT ) : glShadeModel ( GL_SMOOTH );
			break;
		case 99:
			exit ( 5 );
			break;
		default:
			break;
	}
}
