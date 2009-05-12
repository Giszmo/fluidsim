#include "fluid.h"
#include "bicubic_bezier_surface.h"
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
bool TRUFORM_ON = false;
bool showlight = false;
float groundlevel ( float a, float b );
Vektor groundlevelnormal ( float a, float b );
float sink ( float a, float b );
Vektor sinknormal ( float a, float b );
float flat ( float a, float b );
Vektor flatnormal ( float a, float b );
void write_24bitbmp ( const char * filename, unsigned char * pixelinfo );
void read_24bitbmp ( const char * filename, unsigned char ** pixelinfo, unsigned long & width, unsigned long & height );

float ** vertex_array_pointer, ** normal_array_pointer;
unsigned long vn_arraylength, vn_count;

unsigned long ** index_array_pointer;
unsigned long index_arraylength;

float ** particle_coords_pointer, ** particle_colors_pointer;
unsigned long cc_arraylength;

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

Fluid f = Fluid ( 150000,
                  0.5f,
                  0.0f, 0.0f, -9.81f,
                  MINX, MINY, MINZ,
                  MAXX, MAXY, MAXZ );

using namespace std;
unsigned char ** rawdata;
bicubic_bezier_surface boden;

int main ( int argc,char** argv )
{
	int i,j;
	for ( i=0;i<200;i++ )
	{
		cout << ( float ) min ( powf ( ( float ) i/40.0f,0.3333f ),0.99f ) << "f, ";
	}

	vertex_array_pointer=new float*;
	normal_array_pointer=new float*;
	*vertex_array_pointer = new float[100];
	*normal_array_pointer = new float[100];
	vn_arraylength=10;

	index_array_pointer=new unsigned long*;
	*index_array_pointer = new unsigned long[100];
	index_arraylength=10;

	particle_coords_pointer=new float*;
	particle_colors_pointer=new float*;
	*particle_coords_pointer = new float[100];
	*particle_colors_pointer = new float[100];
	cc_arraylength=10;

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

	//f.set_height_function(groundlevel);	f.set_height_function_normal(groundlevelnormal);
	f.set_height_function(sink);	f.set_height_function_normal(sinknormal);
	//f.set_height_function ( flat );	f.set_height_function_normal ( flatnormal );
	float testlist[9*1700];
	int cnt=0;
/*	for ( j=2;j<6;j++ )
	{
		for ( i=0;i<10;i++ )
		{
			testlist[cnt++]=cos ( ( float ) i/10*2*M_PI ) *j*j;		testlist[cnt++]=sin ( ( float ) i/10*2*M_PI ) *j*j;		testlist[cnt++]=15+8*j;
			testlist[cnt++]=cos ( ( float ) ( i+1 ) /10*2*M_PI ) *j*j;	testlist[cnt++]=sin ( ( float ) ( i+1 ) /10*2*M_PI ) *j*j;	testlist[cnt++]=15+8*j;
			testlist[cnt++]=cos ( ( float ) i/10*2*M_PI ) * ( j-1 ) * ( j-1 );	testlist[cnt++]=sin ( ( float ) i/10*2*M_PI ) * ( j-1 ) * ( j-1 );	testlist[cnt++]=15+8* ( j-1 );

			testlist[cnt++]=cos ( ( float ) ( i+1 ) /10*2*M_PI ) *j*j;	testlist[cnt++]=sin ( ( float ) ( i+1 ) /10*2*M_PI ) *j*j;	testlist[cnt++]=15+8*j;
			testlist[cnt++]=cos ( ( float ) ( i+1 ) /10*2*M_PI ) * ( j-1 ) * ( j-1 );		testlist[cnt++]=sin ( ( float ) ( i+1 ) /10*2*M_PI ) * ( j-1 ) * ( j-1 );		testlist[cnt++]=15+8* ( j-1 );
			testlist[cnt++]=cos ( ( float ) i/10*2*M_PI ) * ( j-1 ) * ( j-1 );	testlist[cnt++]=sin ( ( float ) i/10*2*M_PI ) * ( j-1 ) * ( j-1 );	testlist[cnt++]=15+8* ( j-1 );
		}
	}
	f.trianglelist2boundary ( testlist,0,cnt/9-1 );
	
	cnt=0;
	testlist[cnt++]=-3;	testlist[cnt++]=-3;	testlist[cnt++]=23;
	testlist[cnt++]=3;	testlist[cnt++]=-3;	testlist[cnt++]=23;
	testlist[cnt++]=-3;	testlist[cnt++]=3;	testlist[cnt++]=23;
	testlist[cnt++]=-3;	testlist[cnt++]=3;	testlist[cnt++]=23;
	testlist[cnt++]=3;	testlist[cnt++]=-3;	testlist[cnt++]=23;
	testlist[cnt++]=3;	testlist[cnt++]=3;	testlist[cnt++]=23;
	f.trianglelist2shift ( 0,0,-10,testlist,0,1 );

	cnt=0;
	testlist[cnt++]=15;	testlist[cnt++]=15;	testlist[cnt++]=2;
	testlist[cnt++]=15;	testlist[cnt++]=15;	testlist[cnt++]=0;
	testlist[cnt++]=15;	testlist[cnt++]=-15;	testlist[cnt++]=2;
	testlist[cnt++]=15;	testlist[cnt++]=15;	testlist[cnt++]=0;
	testlist[cnt++]=15;	testlist[cnt++]=-15;	testlist[cnt++]=0;
	testlist[cnt++]=15;	testlist[cnt++]=-15;	testlist[cnt++]=2;
	f.trianglelist2teleport ( 0,0,100,testlist,0,1 );

	cnt=0;
	testlist[cnt++]=-5;	testlist[cnt++]=-5;	testlist[cnt++]=100;
	testlist[cnt++]=5;	testlist[cnt++]=-5;	testlist[cnt++]=100;
	testlist[cnt++]=-5;	testlist[cnt++]=5;	testlist[cnt++]=100;
	testlist[cnt++]=5;	testlist[cnt++]=-5;	testlist[cnt++]=100;
	testlist[cnt++]=-5;	testlist[cnt++]=5;	testlist[cnt++]=100;
	testlist[cnt++]=5;	testlist[cnt++]=5;	testlist[cnt++]=100;
	f.trianglelist2setspeed ( 0,0,-1,testlist,0,1 );
*/

	cnt=0;
	testlist[cnt++]=-2;	testlist[cnt++]=-2;	testlist[cnt++]=0;
	testlist[cnt++]=-2;	testlist[cnt++]=2;	testlist[cnt++]=0;
	testlist[cnt++]=2;	testlist[cnt++]=-2;	testlist[cnt++]=0;
	testlist[cnt++]=2;	testlist[cnt++]=-2;	testlist[cnt++]=0;
	testlist[cnt++]=-2;	testlist[cnt++]=2;	testlist[cnt++]=0;
	testlist[cnt++]=2;	testlist[cnt++]=2;	testlist[cnt++]=0;
	f.trianglelist2setspeed ( 0,0,70,testlist,0,1 );

	float tetralist[12*3] = {	-5,-5,55,		10,-5,60,		-5,10,60,		-5,-5,65,
								6,-6,160,		25,-6,160,		6,-25,160,		6,-6,164
	                        };
	f.tetraederlist2liquid ( tetralist,0,0 );

	glutInit ( &argc, argv );
	initMain();
#ifdef WIN32
	SetupATIExtensions();
#endif
	glutTimerFunc ( 10,zeitgeber, 10 );
	glutMainLoop ();

	delete [] *vertex_array_pointer;
	delete [] *normal_array_pointer;
//	delete [] *color_array_pointer;
	delete [] *index_array_pointer;
	delete [] *particle_coords_pointer;
	delete [] *particle_colors_pointer;

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
float sink ( float a, float b )
{
	return .05* ( a*a+b*b );
}
Vektor sinknormal ( float a, float b )
{
	return Vektor ( -.1*a,-.1*b,1 );
}
float flat ( float a, float b )
{
	return 0.0f;
}
Vektor flatnormal ( float a, float b )
{
	return Vektor ( 0.0f,0.0f,1.0f );
}
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
	/*
	glNewList ( TRICHTER, GL_COMPILE );
	for ( j=2;j<6;j++ )
	{
		glBegin ( GL_TRIANGLES );
		glColor4f ( 1,1,1,1 );
		for ( i=0;i<10;i++ )
		{
			glVertex3f ( cos ( ( float ) i/10*2*M_PI ) *j*j,sin ( ( float ) i/10*2*M_PI ) *j*j,15+8*j );
			glVertex3f ( cos ( ( float ) ( i+1 ) /10*2*M_PI ) *j*j,sin ( ( float ) ( i+1 ) /10*2*M_PI ) *j*j,15+8*j );
			glVertex3f ( cos ( ( float ) i/10*2*M_PI ) * ( j-1 ) * ( j-1 ),sin ( ( float ) i/10*2*M_PI ) * ( j-1 ) * ( j-1 ),15+8* ( j-1 ) );

			glVertex3f ( cos ( ( float ) ( i+1 ) /10*2*M_PI ) *j*j,sin ( ( float ) ( i+1 ) /10*2*M_PI ) *j*j,15+8*j );
			glVertex3f ( cos ( ( float ) ( i+1 ) /10*2*M_PI ) * ( j-1 ) * ( j-1 ),sin ( ( float ) ( i+1 ) /10*2*M_PI ) * ( j-1 ) * ( j-1 ),15+8* ( j-1 ) );
			glVertex3f ( cos ( ( float ) i/10*2*M_PI ) * ( j-1 ) * ( j-1 ),sin ( ( float ) i/10*2*M_PI ) * ( j-1 ) * ( j-1 ),15+8* ( j-1 ) );
		}
		glEnd();
	}
	glEndList();*/
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
	Vektor tmp;
	glNewList ( GROUND, GL_COMPILE );
	for ( float i = MINX+0.40* ( MAXX-MINX );i<MINX+0.60* ( MAXX-MINX );i+=0.05* ( MAXX-MINX ) )
	{
		glBegin ( GL_TRIANGLE_STRIP );
		glColor4f ( 1,0,0,1 );
		for ( float j = MINY+0.40* ( MAXY-MINY );j<MINY+0.60* ( MAXY-MINY );j+=0.05* ( MAXY-MINY ) )
		{
			tmp=f.height_function_normal ( i,j ).normed();
			glNormal3f ( tmp.x(),tmp.y(),tmp.z() );
			glVertex3f ( i,j,f.height_function ( i,j ) );
			tmp=f.height_function_normal ( i+0.05* ( MAXX-MINX ),j ).normed();
			glNormal3f ( tmp.x(),tmp.y(),tmp.z() );
			glVertex3f ( i+0.05* ( MAXX-MINX ),j,f.height_function ( i+0.05* ( MAXX-MINX ),j ) );
		}
		glEnd();
	}
	glEndList();

	glNewList ( GROUND_NORMALS, GL_COMPILE );
	for ( float i = MINX+0.40* ( MAXX-MINX );i<MINX+0.60* ( MAXX-MINX );i+=0.05* ( MAXX-MINX ) )
	{
		glBegin ( GL_LINES );
		glColor4f ( 0,1,0,1 );
		for ( float j = MINY+0.40* ( MAXY-MINY );j<MINY+0.60* ( MAXY-MINY );j+=0.05* ( MAXY-MINY ) )
		{
			Vektor tmp=f.height_function_normal ( i,j ).norm ( 0.5f );
			glVertex3f ( i,j,f.height_function ( i,j ) );
			glVertex3f ( i+tmp.x(),j+tmp.y(),f.height_function ( i,j ) +tmp.z() );
		}
		glEnd();
	}
	glEndList();

	glNewList ( BOX, GL_COMPILE );
	glMaterialf ( GL_FRONT,GL_SHININESS,50 );
	glColor3f ( 1,1,1 );
	glBegin ( GL_TRIANGLE_STRIP );
	glVertex3f ( f.minxyz().x(),f.minxyz().y(),f.minxyz().z() );
	glColor3f ( 1,1,0.8 );
	glVertex3f ( f.minxyz().x(),f.minxyz().y(),f.maxxyz().z() );
	glColor3f ( 1,0.8,1 );
	glVertex3f ( f.maxxyz().x(),f.minxyz().y(),f.minxyz().z() );
	glColor3f ( 0.8,1,1 );
	glVertex3f ( f.maxxyz().x(),f.minxyz().y(),f.maxxyz().z() );
	glColor3f ( 1,0.8,0.8 );
	glVertex3f ( f.maxxyz().x(),f.maxxyz().y(),f.minxyz().z() );
	glColor3f ( 0.8,1,0.8 );
	glVertex3f ( f.maxxyz().x(),f.maxxyz().y(),f.maxxyz().z() );
	glColor3f ( 1,0.8,0.8 );
	glVertex3f ( f.minxyz().x(),f.maxxyz().y(),f.minxyz().z() );
	glColor3f ( 1,1,1 );
	glVertex3f ( f.minxyz().x(),f.maxxyz().y(),f.maxxyz().z() );
	glEnd();
	glEndList();
}
void PrintVolumeOfClosedSurface()
{
	unsigned long i;
	float * va;
	unsigned long * ia;
	float volume,dv;
	va = ( *vertex_array_pointer );
	ia = ( *index_array_pointer );
	volume = 0;
	if ( showcells )
	{
		for ( i=0;i<trianglecount;i++ )
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
	else
	{
		cout << "Without Surfacevisualisation, surface is not calculated and thus makes no sence beeing asked its volume." << endl;
	}
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
			cout << "Speed = " << speed << endl;
			break;
		case '/' :
			speed/=1.1;
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
			TRUFORM_ON = !TRUFORM_ON;
			break;
		case 'p' :
			paused = !paused;
			break;
			/*    case 'n' :
			        f.deb_colliderecording() ? f.deb_hidecollide() : f.deb_showcollide();
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
			unsigned int i;
			cout << vn_arraylength << endl;
			for ( i=0;i<vn_count;++i )
			{
				float l;
				l = sqrt ( ( *normal_array_pointer ) [3*i+1]* ( *normal_array_pointer ) [3*i+1]+ ( *normal_array_pointer ) [3*i]* ( *normal_array_pointer ) [3*i]+ ( *normal_array_pointer ) [3*i+2]* ( *normal_array_pointer ) [3*i+2] );
				printf ( "V %f %f %f %f %f %f 0.5 0.5\n", ( *vertex_array_pointer ) [3*i+1], ( *vertex_array_pointer ) [3*i], ( *vertex_array_pointer ) [3*i+2], ( *normal_array_pointer ) [3*i+1]/l, ( *normal_array_pointer ) [3*i]/l, ( *normal_array_pointer ) [3*i+2]/l );
			}
			for ( i=0;i<trianglecount*3;++i )
			{
				printf ( "I %lu\n", ( *index_array_pointer ) [i] );
			}
			break;
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
			        f.deb_print_relevant_grid();
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
//        f.deb_showcollideprev();
////        cout << "Showing collisions with particle " << f.deb_get_bumpingparticle() << endl;
//        break;
//    case GLUT_KEY_RIGHT :
//        f.deb_showcollidenext();
////        cout << "Showing collisions with particle " << f.deb_get_bumpingparticle() << endl;
//        break;
//	#else
//    case GLUT_KEY_LEFT  :
//        g_normal_length/=1.05f;
//        cout << "g_normal_length is set to " << g_normal_length << endl;
////        cout << "Showing collisions with particle " << f.deb_get_bumpingparticle() << endl;
//        break;
//    case GLUT_KEY_RIGHT :
//        g_normal_length*=1.05;
//        cout << "g_normal_length is set to " << g_normal_length << endl;
////        cout << "Showing collisions with particle " << f.deb_get_bumpingparticle() << endl;
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
	static int talt = glutGet ( GLUT_ELAPSED_TIME );
	//float dt;
	int t = glutGet ( GLUT_ELAPSED_TIME );
	talt = t;
	if ( !paused )
		f.progress ( speed );

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
	if ( showcells )
	{
#ifdef WIN32
		if ( TRUFORM_ON )
		{
			glEnable ( GL_PN_TRIANGLES_ATI );
			glPNTrianglesiATI ( GL_PN_TRIANGLES_POINT_MODE_ATI, GL_PN_TRIANGLES_POINT_MODE_CUBIC_ATI );
			glPNTrianglesiATI ( GL_PN_TRIANGLES_NORMAL_MODE_ATI, GL_PN_TRIANGLES_NORMAL_MODE_QUADRATIC_ATI );
			glPNTrianglesiATI ( GL_PN_TRIANGLES_TESSELATION_LEVEL_ATI, 7 );
		}
#endif
		//glCallList(BOX);
		//glCallList(GROUND);

		trianglecount = f.get_surfacegrid (
		                    vertex_array_pointer,
		                    normal_array_pointer,
		                    vn_arraylength,
		                    vn_count,
		                    index_array_pointer,
		                    index_arraylength );
		glEnableClientState ( GL_VERTEX_ARRAY );
		glEnableClientState ( GL_NORMAL_ARRAY );
		glVertexPointer ( 3, GL_FLOAT, 0, *vertex_array_pointer );
		glNormalPointer ( GL_FLOAT, 0, *normal_array_pointer );
//#ifdef debug
//		glEnableClientState(GL_COLOR_ARRAY );
//		glColorPointer(3, GL_FLOAT, 3*sizeof(GLfloat), *normal_array_pointer);
//#else
		glDisableClientState ( GL_COLOR_ARRAY );
//#endif

		glColor4f ( 0.66f,0.66f,1.0f,0.4f );

		glDrawElements ( GL_TRIANGLES, 3*trianglecount, GL_UNSIGNED_INT, *index_array_pointer );
		//glDrawArrays(GL_POINTS,0,vertexcount);
		/*glBegin(GL_LINES);
		for (i=0; i < vertexcount-3; i+=3)
			glVertex3f (va[i],va[i+1],va[i+2]);
		glEnd();*/
#ifdef WIN32
		if ( TRUFORM_ON )
		{
			glDisable ( GL_PN_TRIANGLES_ATI );
		}
#endif
		if ( shownormals )
		{
			glDisable ( GL_LIGHTING );
			float *tmp_norm=*normal_array_pointer;
			float *tmp_vert=*vertex_array_pointer;
			glBegin ( GL_LINES );
			glColor4f ( 1.0f,0.86f,0.86f,1.0f );
			for ( unsigned long i=0; i<vn_count; i++ )
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
	if ( showparticles )
	{
		f.get_particlearray ( particle_coords_pointer/*,particle_colors_pointer*/, cc_arraylength );

		/* set up the array data */
		glVertexPointer ( 3, GL_FLOAT, 3*sizeof ( GLfloat ), *particle_coords_pointer );
		//glColorPointer(3, GL_FLOAT, 3*sizeof(GLfloat), *particle_colors_pointer);

		/* enable vertex arrays */
		glDisableClientState ( GL_NORMAL_ARRAY );
		glEnableClientState ( GL_VERTEX_ARRAY );
		//glEnableClientState( GL_COLOR_ARRAY );

		/* draw a polygon using the arrays sequentially */
		glDrawArrays ( GL_POINTS,0,f.particlecount() );
	}
	glDisable ( GL_LIGHTING );
	glCallList ( BOX );
	glCallList ( TRICHTER );
	glBegin ( GL_LINES );
	glColor4f ( 1,0,0,1 );
	glVertex3f ( 0,0,0 );
	glVertex3f ( 4*f.particleradius(),0,0 );
	glColor4f ( 0,1,0,1 );
	glVertex3f ( 0,0,0 );
	glVertex3f ( 0,4*f.particleradius(),0 );
	glColor4f ( 0,0,1,1 );
	glVertex3f ( 0,0,0 );
	glVertex3f ( 0,0,4*f.particleradius() );
	glEnd();

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

		cout << sekunden << ". Sekunde: " << framecount << "frames, " << f.particlecount() << " particles, " << f.movingparticlecount() << " moving particles, " << trianglecount << "triangles." << endl;
		trianglecount=0;
		sekunden ++;
		framecount = 0;
	}
	glutTimerFunc ( 200,zeitgeber, 1000 );
}

void ReshapeMain ( GLint width, GLint height )
{
	/*	if(f.deb_colliderecording()) {
			GLsizei tmp1 = glutGet(GLUT_WINDOW_WIDTH);
			GLsizei tmp2 = glutGet(GLUT_WINDOW_HEIGHT);
			glutSetWindow(winIdMain);
			glViewport(0, 0, tmp1, tmp2);
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
		    gluPerspective(90.0, (float)tmp1 / tmp2, distanz - f.particleradius(),distanz + f.particleradius());
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
