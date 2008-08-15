#ifndef FLUID_H
#define FLUID_H

//#include "little_tools.h"
#include "particle.h"
#include "vektor.h"
#include <math.h>
#include <float.h>
#include <cstring>
#include <cstdlib>

#include <iostream>
using namespace std;

#ifndef WIN32
#define __int64 long long
#endif

#define BOUNDING_PARTICLES_PER_CELL 2.5f

//8 corners (m | p)(x | y | z)^3
#define m000 0
#define m100 1
#define m110 2
#define m010 3
#define m001 4
#define m101 5
#define m111 6
#define m011 7

//8 corners (m | p)(x | y | z)^3
#define m000bit 1
#define m100bit 2
#define m110bit 4
#define m010bit 8
#define m001bit 16
#define m101bit 32
#define m111bit 64
#define m011bit 128

//12 edges (m | p)(x | y | z)^2
#define mymz 0
#define pxmz 1
#define pymz 2
#define mxmz 3
#define mypz 4
#define pxpz 5
#define pypz 6
#define mxpz 7
#define mxmy 8
#define pxmy 9
#define pxpy 10
#define mxpy 11

struct CUBE //size: 8++1+32+48
{
	unsigned __int64 cell_id;
	Vektor coord_alt;
	unsigned char ecken;//1 bit per corner. 1==inside fluid. 0==outside fluid
	unsigned long nodelist_id[8];//indices of the according corner in the node_list
	unsigned long vn_id[12];//indices of the according vertices/normals
};
struct NODE_DATA
{
	unsigned __int64 cell_id;
	unsigned char count;
	Vektor coord;
};
inline void make_zero_cube_list ( CUBE * cubelist, unsigned long & count )
{
	memset ( cubelist,0x00,count*sizeof ( CUBE ) );
	count=0;
}
inline void resize_node_list ( NODE_DATA ** nodelist, unsigned long count, unsigned long & max_count, unsigned long newsize )
{
	NODE_DATA * tmp = new NODE_DATA[newsize];
	memcpy ( tmp,*nodelist,count*sizeof ( NODE_DATA ) );
	memset ( &tmp[count],0x00, ( newsize-count ) *sizeof ( NODE_DATA ) );
	delete [] *nodelist;
	*nodelist=tmp;
	max_count=newsize;
}

inline void resize_cube_list ( CUBE ** cubelist, unsigned long count, unsigned long & max_count, unsigned long newsize )
{
	CUBE * tmp;
	tmp = new CUBE[newsize];
	memcpy ( tmp,*cubelist,count*sizeof ( CUBE ) );
	memset ( &tmp[count],0x00, ( newsize-count ) *sizeof ( CUBE ) );
	delete [] *cubelist;
	*cubelist=tmp;
	max_count=newsize;
}
inline void resize_unsigned_long_list ( unsigned long ** ullist, unsigned long count, unsigned long & max_count, unsigned long newsize )
{
	unsigned long * tmp;
	tmp = new unsigned long[newsize];
	memcpy ( tmp,*ullist,count*sizeof ( unsigned long ) );
	delete [] *ullist;
	*ullist=tmp;
	max_count=newsize;
}
inline void swap_cubes_without_vn_id ( CUBE & A, CUBE & B )
{
	static CUBE tmp;
	memcpy ( &tmp,   &A, sizeof ( CUBE )-12*sizeof ( long ) );
	memcpy ( &A,   &B, sizeof ( CUBE )-12*sizeof ( long ) );
	memcpy ( &B, &tmp, sizeof ( CUBE )-12*sizeof ( long ) );
}
inline void resize_float ( float ** list, unsigned long count, unsigned long & max_count, unsigned long newsize )
{
	float * tmp;
	tmp = new float[newsize];
	memcpy ( tmp,*list,count*sizeof ( float ) );
	memset ( &tmp[count],0x00, ( newsize-count ) *sizeof ( float ) );
	delete [] *list;
	*list=tmp;
	max_count=newsize;
}
inline void resize_2float ( float ** list1, float ** list2, unsigned long count, unsigned long & max_count, unsigned long newsize )
{
	float * tmp;

	tmp = new float[newsize];
	memcpy ( tmp,*list1,count*sizeof ( float ) );
	memset ( &tmp[count],0x00, ( newsize-count ) *sizeof ( float ) );
	delete [] *list1;
	*list1=tmp;

	tmp = new float[newsize];
	memcpy ( tmp,*list2,count*sizeof ( float ) );
	memset ( &tmp[count],0x00, ( newsize-count ) *sizeof ( float ) );
	delete [] *list2;
	*list2=tmp;
	max_count=newsize;
}

inline void make_zero_node_list ( NODE_DATA * nodelist, unsigned long & count )
{
	memset ( nodelist,0x00,count*sizeof ( NODE_DATA ) );
	count=0;
}
inline float _min ( float a, float b )
{
	return ( a>b ) ?b:a;
}
inline float _max ( float a, float b )
{
	return ( a<b ) ?b:a;
}

int compare ( const void * a, const void * b )
{
	return ( ( unsigned __int64 * ) a > ( unsigned __int64 * ) b ) ?1:-1;
}

class Fluid
{

	private:
//	static int edgeTable[256];
		static int triTable[256][19];
		//static int triTable2[256][13];

		float _t; //time
		//particlebox _particles;
		CUBE ** _tmp_cube_list_pointer;
		unsigned long _tmp_cube_count,_tmp_max_cube_count;
		//Vektor _vistapoint;
		//float _triangle_view_size; //if negative, no refinement will occure

		const unsigned long _maxparticlecount;//limit, das bei Instanzierung festgelegt wird
		const float _particlesize; //Partikeldurchmesser in m (Abstand < Durchmesser => Interaktion)
		const float _minx;//minimales x des Simulationsgebiets in m
		const float _miny;//minimales y des Simulationsgebiets in m
		const float _minz;//minimales z des Simulationsgebiets in m
		const float _cellsize;//4*_particlesize
		const unsigned char _cellsxplus2bits;//resultierend aus Simulationsgebietsgr��e und Partikeldurchmesser werden bits zur Bestimmung des x-Zellindex ben�tigt. +2, da in x-Richtung nicht staggered
		const unsigned char _cellsybits;//resultierend aus Simulationsgebietsgr��e und Partikeldurchmesser werden bits zur Bestimmung des y-Zellindex ben�tigt.
		const unsigned char _cellszbits;//resultierend aus Simulationsgebietsgr��e und Partikeldurchmesser werden bits zur Bestimmung des z-Zellindex ben�tigt.
		const unsigned __int64 _cellsx;//2^(_cellsxplus2bits-2)
		const unsigned __int64 _cellsy;//2^_cellsybits
		const unsigned __int64 _cellsz;//2^_cellszbits
		const float _lenx;//_cellsx*_cellsize
		const float _leny;//_cellsy*_cellsize
		const float _lenz;//_cellsz*_cellsize
		const float _maxx;//_minx+_lenx
		const float _maxy;//_miny+_leny
		const float _maxz;//_minz+_lenz
		const unsigned __int64 _cellsyz;//_cellsy*_cellsz
		const unsigned __int64 _cellsxyz;//_cellsx*_cellsy*_cellsz
		const unsigned char _cellsyzbits;//_cellsybits+_cellszbits
		const unsigned char _cellsxplus2ybits;//_cellsxplus2bits+_cellsybits
		const unsigned char _cellsxplus2yzbits;//_cellsxplus2bits+_cellsybits+_cellszbits
		const unsigned char _maxparticlecountbits;//log2(_maxparticlecount)+1
		const unsigned char _cellsxplus2particlebits;//_cellsxplus2bits+_maxparticlecountbits
		const unsigned char _cellsxplus2yparticlebits;//_cellsxplus2bits+_cellsybits+_maxparticlecountbits
		const unsigned char _cellsxplus2yzparticlebits; //constants throughout one simulation;
		const unsigned __int64 _dx;//Schrittweite zur n�chsten Zelle in x-Richtung
		const unsigned __int64 _dy;//Schrittweite zur n�chsten Zelle in y-Richtung
		const unsigned __int64 _dz;//Schrittweite zur n�chsten Zelle in z-Richtung
		const unsigned __int64 _dxy;//_dx+_dy
		const unsigned __int64 _dxz;//_dx+_dz
		const unsigned __int64 _dyz;//_dy+_dz
		const unsigned __int64 _dxyz;//_dx+_dy+_dz
		const float _particlesize_q;//_particlesize*_particlesize
		const float _m_poly6;//um es nur einmal berechnen zu m�ssen
		const float _m_spiky;//um es nur einmal berechnen zu m�ssen
		const float _m_viscosity;//um es nur einmal berechnen zu m�ssen

		Particle * _particle;//Liste der Partikel
		unsigned __int64 * _sortlist;	//4 Permutationen der Partikelliste (4 staggered grids)
		//[000000000][z-bits][y-bits][x-bits][2 more x-bits][Ballindex-bits]
		//Beispiel 1:
		//39bit cell-id (2^39 cells of edge-length 4*r)
		//23bit particle-id (2^23 = 8M particles)
		//assumeing particles of granularity 1cm we now have approximately
		//4M cm^3 = 4.000 liters of liquid
		//in 2^39*64cm^3
		//Beispiel 2(Computerspiel: Sanit?e Einrichtungen im Level haben "interaktives" Wasser):
		//45bit cell-id (2^15 cells in jeder Dim = 32k^3 cells a 4r kantenlaenge
		//17bit particle-id (2^17 = 125k Baelle)
		//= 60 l (1/2 cm^3 pro 1cm-particle)
		//in 1.5 km^3
		unsigned long ** _tmp_pending_list_pointer;
		unsigned long _tmp_pending_count;
		unsigned long _tmp_max_pending_count;

		NODE_DATA ** _tmp_node_list_pointer;
		unsigned long _tmp_node_list_count,_tmp_max_node_list_count;

		bool _listssorted;
		unsigned long _particlecount;
		unsigned long _particlecount_physik;
		unsigned long _particlecount_boundary;
		unsigned __int64 _BoxIdBitmask;
		unsigned __int64 _ParticleIdBitmask;
		Vektor _a;
		float ( *_height_function ) ( float x, float y );
		Vektor ( *_height_function_normal ) ( float x, float y );
	public:
		//#ifdef debug
		//bool deb_colliderecording;
		//#endif
		Fluid ( const long maxparticlecount,
		        const float particlesize,
		        const float ax, const float ay, const float az,
		        const float minx, const float miny, const float minz,
		        const float maxx, const float maxy, const float maxz ) :
				_maxparticlecount ( maxparticlecount ),
				_particlesize ( particlesize ),
				_minx ( minx ), _miny ( miny ), _minz ( minz ),
				_cellsize ( _particlesize* ( float ) 4 ),
				_cellsxplus2bits ( ( unsigned char ) ( log ( ( double ) ( ( maxx-minx ) / _cellsize ) ) /log ( 2.0 ) +1+2 ) ),
				_cellsybits ( ( unsigned char ) ( log ( ( double ) ( ( maxy-miny ) / _cellsize ) ) /log ( 2.0 ) +1 ) ),
				_cellszbits ( ( unsigned char ) ( log ( ( double ) ( ( maxz-minz ) / _cellsize ) ) /log ( 2.0 ) +1 ) ),
				_cellsx ( ( 1 << _cellsxplus2bits ) >> 2 ),
				_cellsy ( 1 << _cellsybits ),
				_cellsz ( 1 << _cellszbits ),
				_lenx ( ( float ) _cellsx*_cellsize ),
				_leny ( ( float ) _cellsy*_cellsize ),
				_lenz ( ( float ) _cellsz*_cellsize ),
				_maxx ( _minx+_lenx ),
				_maxy ( _miny+_leny ),
				_maxz ( _minz+_lenz ),
				_cellsyz ( _cellsy * _cellsz ),
				_cellsxyz ( _cellsx * _cellsyz ),

				_cellsyzbits ( _cellsybits + _cellszbits ),
				_cellsxplus2ybits ( _cellsxplus2bits  + _cellsybits ),
				_cellsxplus2yzbits ( _cellsxplus2ybits + _cellszbits ),
				_maxparticlecountbits ( ( unsigned char ) ( logf ( ( float ) maxparticlecount ) /logf ( 2.0 ) + 1 ) ),
				_cellsxplus2particlebits ( _cellsxplus2bits + _maxparticlecountbits ),
				_cellsxplus2yparticlebits ( _cellsxplus2particlebits  + _cellsybits ),
				_cellsxplus2yzparticlebits ( _cellsxplus2yparticlebits + _cellszbits ),

				_dx ( ( unsigned __int64 ) 1<< ( _maxparticlecountbits + 2 ) ),//////////////+1????????
				_dy ( ( unsigned __int64 ) 1<< ( _cellsxplus2particlebits ) ),
				_dz ( ( unsigned __int64 ) 1<< ( _cellsxplus2yparticlebits ) ),
				_dxy ( _dx | _dy ),
				_dxz ( _dx | _dz ),
				_dyz ( _dy | _dz ),
				_dxyz ( _dx | _dy | _dz ),
				_particlesize_q ( _particlesize*_particlesize ),
				_m_poly6 ( 315.0/64.0/M_PI/pow ( _particlesize,9 ) ),
				_m_spiky ( 15.0/M_PI/pow ( _particlesize,6 ) ),
				_m_viscosity ( ( float ) 15/2.0/M_PI/pow ( _particlesize,3 ) ),
				_a ( Vektor ( ax,ay,az ) )
		{
			_particle = new Particle[maxparticlecount];
			_listssorted = false;

			_particlecount=_particlecount_physik=_particlecount_boundary=0;
			if ( _cellsxplus2yzparticlebits > 64 )
			{
				cout << "2 + CellsXBits(" << ( int ) _cellsxplus2bits-2 <<
				") + CellsYBits(" << ( int ) _cellsybits <<
				") + CellsZBits(" << ( int ) _cellszbits <<
				") + MaxParticleCountBits(" << ( int ) _maxparticlecountbits <<
				") = " << ( ( int ) _cellsxplus2yzparticlebits ) <<
				" > 64. Please either choose bigger particles or smaller space." <<
				" Dimensions are rounded up to be power-of-two multiples of (4*BallSize)" << endl;
				exit ( 1 );
			}

			_sortlist = new unsigned __int64[ _maxparticlecount * 4 ];
			_ParticleIdBitmask = ( 0xffffffffffffffffULL ) >> ( 64-_maxparticlecountbits );
			_BoxIdBitmask =
			    ( ( ( ( ( unsigned __int64 ) 1 <<  _cellszbits )-1 ) ) << _cellsxplus2yparticlebits ) |
			    ( ( ( ( ( unsigned __int64 ) 1 <<  _cellsybits )-1 ) ) << _cellsxplus2particlebits ) |
			    ( ( ( ( ( unsigned __int64 ) 1 << ( _cellsxplus2bits-2 ) )-1 ) ) << ( _maxparticlecountbits + 2 ) );

			_tmp_pending_count=_tmp_max_pending_count=100;
			_tmp_pending_list_pointer = new unsigned long*;
			*_tmp_pending_list_pointer = new unsigned long[_tmp_max_pending_count];

			_tmp_node_list_count=_tmp_max_node_list_count=100;
			_tmp_node_list_pointer = new NODE_DATA*;
			*_tmp_node_list_pointer = new NODE_DATA[_tmp_max_node_list_count];
			make_zero_node_list ( *_tmp_node_list_pointer,_tmp_node_list_count );
			//_vistapoint=Vektor(maxx/2,maxy/2,maxz/2);
			//_triangle_view_size=3;
			_tmp_cube_count = _tmp_max_cube_count = 100;
			_tmp_cube_list_pointer = new CUBE*;
			*_tmp_cube_list_pointer = new CUBE[_tmp_max_cube_count];
			make_zero_cube_list ( ( *_tmp_cube_list_pointer ),_tmp_cube_count );//also sets _tmp_cube_count to 0.

			if ( ( maxx <= minx ) || ( maxy <= miny ) || ( maxz <= minz ) || ( particlesize <= 0 ) )
			{
				cout << "strange setup (min>max or particlesize<=0)!" << endl;
				exit ( 1 );
			}

		}
		~Fluid()
		{
//		fluid2file("test.txt");
//		file2fluid("test.txt");
			delete [] _particle;
			delete [] _sortlist;
			delete [] *_tmp_pending_list_pointer;
			delete [] *_tmp_node_list_pointer;
			delete [] *_tmp_cube_list_pointer;
		}
		float height_function ( float x, float y )
		{
			return _height_function ( x,y );
		}
		Vektor height_function_normal ( float x, float y )
		{
			return _height_function_normal ( x,y );
		}
		void set_height_function ( float ( *height_function ) ( float x, float y ) )
		{
			_height_function=height_function;
		}
		void set_height_function_normal ( Vektor ( *height_function_normal ) ( float x, float y ) )
		{
			_height_function_normal=height_function_normal;
		}
		//void fluid2file(char * filename) {
		//	_particles.fluid2file(filename);
		//}
		int progress ( const float dt )
		{
			if ( dt>0 )
			{
				unsigned long i;
				static float t=0;
				t+=dt;
				putparticles2cells();
				sortparticlelists();
				for ( i=0; i<_particlecount_physik; i++ )
				{
					_particle[i].resetrho();
					_particle[i].acc ( _a );
				}
				for ( i=0; i<_particlecount_physik; i++ )
				{
					collidewithneighboursforrho ( i );
				}
				for ( i=0; i<_particlecount_physik; i++ )
				{
					collidewithneighbours ( i );
				}

				for ( i=0; i<_particlecount_physik; i++ )
				{
					_particle[i].move ( dt );
					if ( _particle[i].x().z() < _height_function ( _particle[i].x().x(),_particle[i].x().y() ) )
					{
						_particle[i].x().setz ( _height_function ( _particle[i].x().x(),_particle[i].x().y() ) +0.002f );
						Vektor tmp_norm = _height_function_normal ( _particle[i].x().x(),_particle[i].x().y() );
						_particle[i].setv ( _particle[i].v() +	tmp_norm.normed() * ( ( _particle[i].v() *tmp_norm.norm ( -1.8f ) ) ) );
					}
				}
				return 1;
			}
			return 0;
		}
		unsigned long particlecount ( void )
		{
			return _particlecount;
		}
		unsigned long boundaryparticlecount ( void )
		{
			return _particlecount_boundary;
		}
		unsigned long movingparticlecount ( void )
		{
			return _particlecount_physik;
		}
		float particleradius ( void )
		{
			return _particlesize;
		}
		Vektor minxyz ( void )
		{
			return Vektor ( _minx,_miny,_minz );
		}
		Vektor maxxyz ( void )
		{
			return Vektor ( _maxx,_maxy,_maxz );
		}
		void get_particlearray ( float ** VertexArray, unsigned long & arraylength )
		{
			unsigned long i;
			float * va;
			if ( arraylength < 3*particlecount() )
			{
				resize_float ( VertexArray,0,arraylength,3*particlecount() +3 );
			}
			va = *VertexArray;
			for ( i = 0; i<particlecount();i++ )
			{
				va[i*3]=get_particlecoord ( i ).x();
				va[i*3+1]=get_particlecoord ( i ).y();
				va[i*3+2]=get_particlecoord ( i ).z();
			}
		}
		int get_surfacegrid ( float ** VertexArrayPointer, float ** NormalArrayPointer, unsigned long & vn_arraylength, unsigned long & vn_count, unsigned long ** triangleindicesPointer, unsigned long & ti_arraylength )   // returns trianglecount
		{
			unsigned long j=0,i,tii=0;
			unsigned long * tmp_ti;
			float * tmp_na;
			CUBE * tmp_cl;
			tmp_ti = *triangleindicesPointer;

			//vn_count = _particles.get_surface_cubes(_tmp_cube_list_pointer,_tmp_cube_count,_tmp_max_cube_count,VertexArrayPointer,NormalArrayPointer,vn_arraylength);
			fill_node_list ( /*,maxparticlesincell*/ );// gets densities out of sorted list of particles. particles in same cell contribute to cells density.
			create_cubes();//creates all needed cubes with information on wich corners are inside/outside fluid and per-corner-index of according cell
			//remove_cubes_inside_liquid(cube_list,cube_count);
			quicksortcubes ( ( *_tmp_cube_list_pointer ),0,_tmp_cube_count-1 );//as cubes got inserted unsorted, the list needs to be sorted
			vn_count = make_normals_and_vertices_for_cubes ( _tmp_cube_list_pointer, _tmp_cube_count,VertexArrayPointer,NormalArrayPointer,vn_arraylength );

			tmp_na=*NormalArrayPointer;
			for ( i=0; i<vn_count; ++i )
			{
				float l=sqrtf ( tmp_na[3*i]*tmp_na[3*i]+tmp_na[3*i+1]*tmp_na[3*i+1]+tmp_na[3*i+2]*tmp_na[3*i+2] );
				tmp_na[3*i]/=l;
				tmp_na[3*i+1]/=l;
				tmp_na[3*i+2]/=l;
			}

			tmp_cl = *_tmp_cube_list_pointer;
			for ( j=1;j<_tmp_cube_count;++j )
			{
				if ( tii + 72 >= ti_arraylength )
				{
					resize_unsigned_long_list ( triangleindicesPointer,tii,ti_arraylength, ( unsigned long ) ( ti_arraylength*1.2 + 72 ) );
					tmp_ti = *triangleindicesPointer;
				}
				for ( i=0;triTable[tmp_cl[j].ecken][i]!=-1;i+=3 )
				{
					//#ifdef debug
					//if((tii+3)>=ti_arraylength)
					//myexit(5453);
					//#endif
					if ( tmp_cl[j].vn_id[triTable[tmp_cl[j].ecken][i  ]]*tmp_cl[j].vn_id[triTable[tmp_cl[j].ecken][i+1]]*tmp_cl[j].vn_id[triTable[tmp_cl[j].ecken][i+2]] == 0 )
					{
						cout << "j: " << j << "\ti: " << i << "\ttmp_cl[j].ecken: " << ( int ) tmp_cl[j].ecken << "\teck1: " << tmp_cl[j].vn_id[triTable[tmp_cl[j].ecken][i  ]] << "\teck2: " << tmp_cl[j].vn_id[triTable[tmp_cl[j].ecken][i+1]] << "\teck3: " << tmp_cl[j].vn_id[triTable[tmp_cl[j].ecken][i+2]] << endl;
					}
					else
					{
						tmp_ti[tii++] = tmp_cl[j].vn_id[triTable[tmp_cl[j].ecken][i  ]];
						tmp_ti[tii++] = tmp_cl[j].vn_id[triTable[tmp_cl[j].ecken][i+1]];
						tmp_ti[tii++] = tmp_cl[j].vn_id[triTable[tmp_cl[j].ecken][i+2]];
					}
				}
			}
			return tii/3;
		}
		void trianglelist2boundary ( float * trianglelist, unsigned int firsttriangleid, unsigned int lasttriangleid )
		{
			Vektor v[3];
			for ( unsigned int i = firsttriangleid; i<=lasttriangleid;++i )
			{
				v[0]=Vektor ( trianglelist[9*i+0],trianglelist[9*i+1],trianglelist[9*i+2] );
				v[1]=Vektor ( trianglelist[9*i+3],trianglelist[9*i+4],trianglelist[9*i+5] );
				v[2]=Vektor ( trianglelist[9*i+6],trianglelist[9*i+7],trianglelist[9*i+8] );
				triangle2boundary ( v );
			}
		}
		void trianglelist2teleport ( float target_x, float target_y, float target_z, float * trianglelist, unsigned int firsttriangleid, unsigned int lasttriangleid )
		{
			Vektor v[3];
			for ( unsigned int i = firsttriangleid; i<=lasttriangleid;++i )
			{
				v[0]=Vektor ( trianglelist[9*i+0],trianglelist[9*i+1],trianglelist[9*i+2] );
				v[1]=Vektor ( trianglelist[9*i+3],trianglelist[9*i+4],trianglelist[9*i+5] );
				v[2]=Vektor ( trianglelist[9*i+6],trianglelist[9*i+7],trianglelist[9*i+8] );
				triangle2teleport ( v,Vektor ( target_x,target_y,target_z ) );
			}
		}
		void trianglelist2shift ( float shift_x, float shift_y, float shift_z, float * trianglelist, unsigned int firsttriangleid, unsigned int lasttriangleid )
		{
			Vektor v[3];
			for ( unsigned int i = firsttriangleid; i<=lasttriangleid;++i )
			{
				v[0]=Vektor ( trianglelist[9*i+0],trianglelist[9*i+1],trianglelist[9*i+2] );
				v[1]=Vektor ( trianglelist[9*i+3],trianglelist[9*i+4],trianglelist[9*i+5] );
				v[2]=Vektor ( trianglelist[9*i+6],trianglelist[9*i+7],trianglelist[9*i+8] );
				triangle2shift ( v,Vektor ( shift_x,shift_y,shift_z ) );
			}
		}
		void trianglelist2setspeed ( float speed_x, float speed_y, float speed_z, float * trianglelist, unsigned int firsttriangleid, unsigned int lasttriangleid )
		{
			Vektor v[3];
			for ( unsigned int i = firsttriangleid; i<=lasttriangleid;++i )
			{
				v[0]=Vektor ( trianglelist[9*i+0],trianglelist[9*i+1],trianglelist[9*i+2] );
				v[1]=Vektor ( trianglelist[9*i+3],trianglelist[9*i+4],trianglelist[9*i+5] );
				v[2]=Vektor ( trianglelist[9*i+6],trianglelist[9*i+7],trianglelist[9*i+8] );
				triangle2setspeed ( v,Vektor ( speed_x,speed_y,speed_z ) );
			}
		}
		void tetraederlist2liquid ( float * tetraederlist, unsigned int firsttetraid, unsigned int lasttetraid )
		{
			Vektor v[4];
			for ( unsigned int i = firsttetraid; i<=lasttetraid;++i )
			{
				v[0]=Vektor ( tetraederlist[12*i+ 0],tetraederlist[12*i+ 1],tetraederlist[12*i+ 2] );
				v[1]=Vektor ( tetraederlist[12*i+ 3],tetraederlist[12*i+ 4],tetraederlist[12*i+ 5] );
				v[2]=Vektor ( tetraederlist[12*i+ 6],tetraederlist[12*i+ 7],tetraederlist[12*i+ 8] );
				v[3]=Vektor ( tetraederlist[12*i+ 9],tetraederlist[12*i+10],tetraederlist[12*i+11] );
				tetraeder2particle ( v );
			}
		}
	void fluid2file(char * filename) {
//		//TODO: implement
//          fstream binary_file(filename,ios::out|ios::binary);
//		  binary_file << particlecount();
///*          binary_file.write(reinterpret_cast<const char *>(&_cellsize),sizeof(float));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsx),sizeof(__int64));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsxplus2bits),sizeof(char));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsxplus2particlebits),sizeof(char));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsxplus2ybits),sizeof(char));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsxplus2yparticlebits),sizeof(char));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsxplus2yzbits),sizeof(char));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsxplus2yzparticlebits),sizeof(char));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsxyz),sizeof(__int64));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsy),sizeof(__int64));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsybits),sizeof(char));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsyz),sizeof(__int64));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsyzbits),sizeof(char));
//          binary_file.write(reinterpret_cast<const char *>(&_cellsz),sizeof(__int64));
//          binary_file.write(reinterpret_cast<const char *>(&_cellszbits),sizeof(char));*/
//          binary_file.close();
	}
	void file2fluid(char * filename) {
//		//TODO: implement
	}
//
	private:
		Vektor get_particlecoord ( const unsigned long n )
		{
			return _particle[n].x();
		}
		void triangle2boundary ( Vektor * v )
		{
			Vektor n;
			n= ( v[1]-v[0] ) % ( v[2]-v[0] );

			float A = n.abs() /2;
			float A0 = _cellsize * _cellsize;
			float ratio = BOUNDING_PARTICLES_PER_CELL*A/A0;
			if ( ratio > 1.0002f )
			{
				Vektor v2[3];
				memcpy ( v2,v,sizeof ( Vektor ) *3 );
				float teilverhaeltnis;
				teilverhaeltnis = floorf ( ( 1.0f+ratio ) /2.0f ) /ratio;
				if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[0] ).absabs() )
				{
					if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[1] ).absabs() )  //01
					{
						v2[0]=v[1]=v[0]* ( 1.0f-teilverhaeltnis ) +v[1]*teilverhaeltnis;
					}
					else  //12
					{
						v2[1]=v[2]=v[1]* ( 1.0f-teilverhaeltnis ) +v[2]*teilverhaeltnis;
					}
				}
				else
				{
					if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[0] ).absabs() )  //01
					{
						v2[0]=v[1]=v[0]* ( 1.0f-teilverhaeltnis ) +v[1]*teilverhaeltnis;
					}
					else  //02
					{
						v2[0]=v[2]=v[0]* ( 1.0f-teilverhaeltnis ) +v[2]*teilverhaeltnis;
					}
				}
				triangle2boundary ( v2 );
				triangle2boundary ( v );
			}
			else
			{
				Vektor coord = ( v[0]+v[1]+v[2] ) *0.33333333333f;
				n=n.normed();
				newcontrollparticle ( coord.x(),coord.y(),coord.z(),n.x(),n.y(),n.z(),n.x(),n.y(),n.z(),boundary );
			}
		}
		void triangle2teleport ( Vektor * v, Vektor target )
		{
			Vektor n;
			n= ( v[1]-v[0] ) % ( v[2]-v[0] );

			float A = n.abs() /2;
			float A0 = _cellsize * _cellsize;
			float ratio = BOUNDING_PARTICLES_PER_CELL*A/A0;
			if ( ratio > 1.0002f )
			{
				Vektor v2[3];
				memcpy ( v2,v,sizeof ( Vektor ) *3 );
				float teilverhaeltnis;
				teilverhaeltnis = floorf ( ( 1.0f+ratio ) /2.0f ) /ratio;
				if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[0] ).absabs() )
				{
					if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[1] ).absabs() )  //01
					{
						v2[0]=v[1]=v[0]* ( 1.0f-teilverhaeltnis ) +v[1]*teilverhaeltnis;
					}
					else  //12
					{
						v2[1]=v[2]=v[1]* ( 1.0f-teilverhaeltnis ) +v[2]*teilverhaeltnis;
					}
				}
				else
				{
					if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[0] ).absabs() )  //01
					{
						v2[0]=v[1]=v[0]* ( 1.0f-teilverhaeltnis ) +v[1]*teilverhaeltnis;
					}
					else  //02
					{
						v2[0]=v[2]=v[0]* ( 1.0f-teilverhaeltnis ) +v[2]*teilverhaeltnis;
					}
				}
				triangle2teleport ( v2,target );
				triangle2teleport ( v,target );
			}
			else
			{
				Vektor coord = ( v[0]+v[1]+v[2] ) *0.33333333333f;
				newcontrollparticle ( coord.x(),coord.y(),coord.z(),n.x(),n.y(),n.z(),target.x(),target.y(),target.z(),teleport );
			}
		}
		void triangle2shift ( Vektor * v, Vektor vshift )
		{
			Vektor n;
			n= ( v[1]-v[0] ) % ( v[2]-v[0] );

			float A = n.abs() /2;
			float A0 = _cellsize * _cellsize;
			float ratio = BOUNDING_PARTICLES_PER_CELL*A/A0;
			if ( ratio > 1.0002f )
			{
				Vektor v2[3];
				memcpy ( v2,v,sizeof ( Vektor ) *3 );
				float teilverhaeltnis;
				teilverhaeltnis = floorf ( ( 1.0f+ratio ) /2.0f ) /ratio;
				if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[0] ).absabs() )
				{
					if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[1] ).absabs() )  //01
					{
						v2[0]=v[1]=v[0]* ( 1.0f-teilverhaeltnis ) +v[1]*teilverhaeltnis;
					}
					else  //12
					{
						v2[1]=v[2]=v[1]* ( 1.0f-teilverhaeltnis ) +v[2]*teilverhaeltnis;
					}
				}
				else
				{
					if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[0] ).absabs() )  //01
					{
						v2[0]=v[1]=v[0]* ( 1.0f-teilverhaeltnis ) +v[1]*teilverhaeltnis;
					}
					else  //02
					{
						v2[0]=v[2]=v[0]* ( 1.0f-teilverhaeltnis ) +v[2]*teilverhaeltnis;
					}
				}
				triangle2shift ( v2,vshift );
				triangle2shift ( v,vshift );
			}
			else
			{
				Vektor coord = ( v[0]+v[1]+v[2] ) *0.33333333333f;
				newcontrollparticle ( coord.x(),coord.y(),coord.z(),n.x(),n.y(),n.z(),vshift.x(),vshift.y(),vshift.z(),shift );
			}
		}
		void triangle2setspeed ( Vektor * v, Vektor speed )
		{
			Vektor n;
			n= ( v[1]-v[0] ) % ( v[2]-v[0] );

			float A = n.abs() /2;
			float A0 = _cellsize * _cellsize;
			float ratio = BOUNDING_PARTICLES_PER_CELL*A/A0;
			if ( ratio > 1.0002f )
			{
				Vektor v2[3];
				memcpy ( v2,v,sizeof ( Vektor ) *3 );
				float teilverhaeltnis;
				teilverhaeltnis = floorf ( ( 1.0f+ratio ) /2.0f ) /ratio;
				if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[0] ).absabs() )
				{
					if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[1] ).absabs() )  //01
					{
						v2[0]=v[1]=v[0]* ( 1.0f-teilverhaeltnis ) +v[1]*teilverhaeltnis;
					}
					else  //12
					{
						v2[1]=v[2]=v[1]* ( 1.0f-teilverhaeltnis ) +v[2]*teilverhaeltnis;
					}
				}
				else
				{
					if ( ( v[1]-v[0] ).absabs() > ( v[2]-v[0] ).absabs() )  //01
					{
						v2[0]=v[1]=v[0]* ( 1.0f-teilverhaeltnis ) +v[1]*teilverhaeltnis;
					}
					else  //02
					{
						v2[0]=v[2]=v[0]* ( 1.0f-teilverhaeltnis ) +v[2]*teilverhaeltnis;
					}
				}
				triangle2setspeed ( v2,speed );
				triangle2setspeed ( v,speed );
			}
			else
			{
				Vektor coord = ( v[0]+v[1]+v[2] ) *0.33333333333f;
				newcontrollparticle ( coord.x(),coord.y(),coord.z(),n.x(),n.y(),n.z(),speed.x(),speed.y(),speed.z(),setspeed );
			}
		}

		void fill_node_list ( /*, unsigned short & maxparticlesincell*/ )
		{
			unsigned __int64 actualboxid;
			NODE_DATA * tmp_nl = *_tmp_node_list_pointer;
			static unsigned int maxparticlesincell=0;
			unsigned int cnt;
			unsigned long i=0;
			Vektor coord;
			make_zero_node_list ( *_tmp_node_list_pointer,_tmp_node_list_count );
			while ( i<_particlecount )
			{
				if ( sortlist2particleid ( _sortlist[i] ) <_particlecount_physik )
				{
					actualboxid = _sortlist[i] & _BoxIdBitmask;
					cnt=1;
					coord=_particle[sortlist2particleid ( _sortlist[i] ) ].x();
					++i;
					while ( actualboxid== ( _sortlist[i] & _BoxIdBitmask ) )
					{
						if ( sortlist2particleid ( _sortlist[i] ) <_particlecount_physik )
						{
							++cnt;
							coord += _particle[sortlist2particleid ( _sortlist[i] ) ].x();
							++i;
						}
						else
						{
							++i;
						}
					}
					if ( _tmp_node_list_count + 12 >= _tmp_max_node_list_count )
					{
						resize_node_list ( _tmp_node_list_pointer,_tmp_node_list_count,_tmp_max_node_list_count, ( unsigned long ) ( _tmp_max_node_list_count*1.2+12 ) );
						tmp_nl = *_tmp_node_list_pointer;
					}

					tmp_nl[++_tmp_node_list_count].cell_id = actualboxid;
					tmp_nl[  _tmp_node_list_count].count  = cnt;
					tmp_nl[  _tmp_node_list_count].coord  = coord / ( float ) cnt;
					if ( maxparticlesincell<cnt )
					{
						maxparticlesincell=cnt;
						cout << "Maxparticles in Cell: " << ( int ) cnt << endl;
					}
				}
				else
				{
					++i;
				}
			}
			++_tmp_node_list_count;//node Null is not a valid node.
		}
		void create_cubes()
		{
			CUBE * tmp_cl;
			NODE_DATA * tmp_nl;
			unsigned long * tmp_pl;
			tmp_cl = *_tmp_cube_list_pointer;
			tmp_nl = *_tmp_node_list_pointer;
			tmp_pl = *_tmp_pending_list_pointer;
			unsigned long i,j;
			bool found_neighbourcube[8];
			unsigned __int64 actualboxid;
			make_zero_cube_list ( tmp_cl,_tmp_cube_count );// -> cube_count=0
			_tmp_pending_count=0;
			//TODO:?make_zero_cube_list is not needed, as cubes have already been created!!?
			for ( i=1; i<_tmp_node_list_count; ++i )
			{
				if ( ( _tmp_cube_count + 20 ) >= _tmp_max_cube_count )
				{
					resize_cube_list ( _tmp_cube_list_pointer,_tmp_cube_count,_tmp_max_cube_count, ( unsigned long ) ( _tmp_max_cube_count*1.2 + 20 ) );
					tmp_cl= ( *_tmp_cube_list_pointer );
				}
				actualboxid=tmp_nl[i].cell_id;
				memset ( found_neighbourcube, 0x00,8*sizeof ( bool ) );
				tmp_cl[_tmp_cube_count].cell_id=actualboxid;
				tmp_cl[_tmp_cube_count].nodelist_id[m000]=i;
				tmp_cl[_tmp_cube_count].ecken=m000bit;

				for ( j=0;j<_tmp_pending_count;j++ )
				{
					unsigned __int64 tmp=actualboxid - tmp_cl[tmp_pl[j]].cell_id;
					if ( tmp & ( 0xffffffffffffffffULL ^ _dxyz ) )
					{
						if ( ( tmp_cl[tmp_pl[j]].cell_id + _dxyz ) < actualboxid )
						{
							tmp_pl[j--]=tmp_pl[--_tmp_pending_count];//bei geschachtelten if: _tmp_pending[j--]=_tmp_pending[_tmp_pending_count--];
						}
					}
					else
					{
						if ( tmp == _dx )
						{
							found_neighbourcube[m100]=true;
							tmp_cl[tmp_pl[j]].ecken |= m100bit;
							tmp_cl[tmp_pl[j]].nodelist_id[m100] = i;
						}
						else if ( tmp == _dy )
						{
							found_neighbourcube[m010]=true;
							tmp_cl[tmp_pl[j]].ecken |= m010bit;
							tmp_cl[tmp_pl[j]].nodelist_id[m010] = i;
						}
						else if ( tmp == _dz )
						{
							found_neighbourcube[m001]=true;
							tmp_cl[tmp_pl[j]].ecken |= m001bit;
							tmp_cl[tmp_pl[j]].nodelist_id[m001] = i;
						}
						else if ( tmp == _dxy )
						{
							found_neighbourcube[m110]=true;
							tmp_cl[tmp_pl[j]].ecken |= m110bit;
							tmp_cl[tmp_pl[j]].nodelist_id[m110] = i;
						}
						else if ( tmp == _dxz )
						{
							found_neighbourcube[m101]=true;
							tmp_cl[tmp_pl[j]].ecken |= m101bit;
							tmp_cl[tmp_pl[j]].nodelist_id[m101] = i;
						}
						else if ( tmp == _dyz )
						{
							found_neighbourcube[m011]=true;
							tmp_cl[tmp_pl[j]].ecken |= m011bit;
							tmp_cl[tmp_pl[j]].nodelist_id[m011] = i;
						}
						else if ( tmp == _dxyz )
						{
							found_neighbourcube[m111]=true;
							tmp_cl[tmp_pl[j]].ecken |= m111bit;
							tmp_cl[tmp_pl[j]].nodelist_id[m111] = i;
						}
					}
				}
				tmp_pl[_tmp_pending_count++]=_tmp_cube_count;//new cube to _tmp_pending-list
				++_tmp_cube_count;
				//if neighbour not found: create cube, set corner-data accordingly, make it _tmp_pending.
				if ( _tmp_pending_count+10 >= _tmp_max_pending_count )
				{
					resize_unsigned_long_list ( _tmp_pending_list_pointer,_tmp_pending_count,_tmp_max_pending_count, ( unsigned long ) ( _tmp_max_pending_count*1.2+10 ) );
					tmp_pl=*_tmp_pending_list_pointer;
				}
				if ( !found_neighbourcube[m001] )
				{
					tmp_cl[_tmp_cube_count].cell_id=actualboxid-_dz;
					tmp_cl[_tmp_cube_count].nodelist_id[m001]=i;
					tmp_cl[_tmp_cube_count].ecken=m001bit;
					tmp_pl[_tmp_pending_count++]=_tmp_cube_count++;
				}
				if ( !found_neighbourcube[m010] )
				{
					tmp_cl[_tmp_cube_count].cell_id=actualboxid-_dy;
					tmp_cl[_tmp_cube_count].nodelist_id[m010]=i;
					tmp_cl[_tmp_cube_count].ecken=m010bit;
					tmp_pl[_tmp_pending_count++]=_tmp_cube_count++;
				}
				if ( !found_neighbourcube[m011] )
				{
					tmp_cl[_tmp_cube_count].cell_id=actualboxid-_dyz;
					tmp_cl[_tmp_cube_count].nodelist_id[m011]=i;
					tmp_cl[_tmp_cube_count].ecken=m011bit;
					tmp_pl[_tmp_pending_count++]=_tmp_cube_count++;
				}
				if ( !found_neighbourcube[m100] )
				{
					tmp_cl[_tmp_cube_count].cell_id=actualboxid-_dx;
					tmp_cl[_tmp_cube_count].nodelist_id[m100]=i;
					tmp_cl[_tmp_cube_count].ecken=m100bit;
					tmp_pl[_tmp_pending_count++]=_tmp_cube_count++;
				}
				if ( !found_neighbourcube[m101] )
				{
					tmp_cl[_tmp_cube_count].cell_id=actualboxid-_dxz;
					tmp_cl[_tmp_cube_count].nodelist_id[m101]=i;
					tmp_cl[_tmp_cube_count].ecken=m101bit;
					tmp_pl[_tmp_pending_count++]=_tmp_cube_count++;
				}
				if ( !found_neighbourcube[m110] )
				{
					tmp_cl[_tmp_cube_count].cell_id=actualboxid-_dxy;
					tmp_cl[_tmp_cube_count].nodelist_id[m110]=i;
					tmp_cl[_tmp_cube_count].ecken=m110bit;
					tmp_pl[_tmp_pending_count++]=_tmp_cube_count++;
				}
				if ( !found_neighbourcube[m111] )
				{
					tmp_cl[_tmp_cube_count].cell_id=actualboxid-_dxyz;
					tmp_cl[_tmp_cube_count].nodelist_id[m111]=i;
					tmp_cl[_tmp_cube_count].ecken=m111bit;
					tmp_pl[_tmp_pending_count++]=_tmp_cube_count++;
				}
			}
		}
		void quicksortcubes ( CUBE * cube_list, const unsigned long lo, const unsigned long hi )
		{
			unsigned long i=lo, j=hi;
			unsigned __int64 x=cube_list[ ( lo+hi ) /2].cell_id;
			while ( i + 1 <= j + 1 )
			{
				while ( cube_list[i].cell_id < x )
					i++;
				while ( cube_list[j].cell_id > x )
					j--;
				if ( i + 1 <= j + 1 )
				{
					swap_cubes_without_vn_id ( cube_list[i++], cube_list[j--] );
				}
			}
			if ( lo + 1 < j  + 1 )
				quicksortcubes ( cube_list,lo, j );
			if ( i  + 1 < hi + 1 )
				quicksortcubes ( cube_list,i, hi );
		}
		void remove_cubes_inside_liquid_scrambleing_order ( CUBE ** cube_list_pointer, unsigned long & cube_count )
		{
			CUBE * tmp_cl;
			tmp_cl= ( *cube_list_pointer );
			unsigned long i=0;
			unsigned long j=cube_count - 1;
			while ( i<j )
			{
				while ( ( tmp_cl[i].ecken!=0xFF ) && ( i<j ) )
				{
					++i;
				}
				while ( ( tmp_cl[j].ecken==0xFF ) && ( i<j ) )
				{
					--j;
				}
				if ( i<j )
				{
					tmp_cl[i++]=tmp_cl[j--];
				}
			}
			cube_count = j+1;
		}
		void remove_cubes_inside_liquid_keeping_order ( CUBE ** cube_list_pointer, unsigned long & cube_count )
		{
			cout << "remove_cubes_inside_liquid_keeping_order not tested so far" << endl;
			CUBE * tmp_cl;
			tmp_cl= ( *cube_list_pointer );
			unsigned long i=0;
			unsigned long j=1;
			while ( j<cube_count )
			{
				while ( ( tmp_cl[i].ecken!=0xFF ) && ( i<cube_count ) )
				{
					++i;
				}
				if ( i>=j )
				{
					j=i+1;
				}
				while ( ( tmp_cl[j].ecken==0xFF ) && ( j<cube_count ) )
				{
					++j;
				}
				if ( j<cube_count )
				{
					tmp_cl[i++]=tmp_cl[j++];
				}
			}
			cube_count = j+1;
		}
		unsigned long make_normals_and_vertices_for_cubes (	CUBE ** cube_list_pointer,
		        const unsigned long cube_count,
		        float ** vertex_list_pointer,
		        float ** normal_list_pointer,
		        unsigned long & vn_arraylength )
		{
			CUBE * tmp_cl;
			unsigned long * tmp_pl;
			NODE_DATA * tmp_nl;
			float * tmp_va;
			float * tmp_na;
			static float dichte2size[200] = {0.0f,0.292438f,0.36844f,0.421753f,0.464195f,0.500035f,0.531363f,0.559377f,0.584835f,0.60825f,0.62999f,0.650324f,0.66946f,0.68756f,0.704755f,0.721148f,0.736829f,0.751869f,0.76633f,0.780265f,0.793719f,0.806732f,0.819338f,0.831567f,0.843447f,0.855001f,0.866252f,0.877217f,0.887915f,0.898361f,0.908569f,0.918553f,0.928325f,0.937895f,0.947273f,0.95647f,0.965493f,0.97435f,0.983049f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f,0.99f
			                                };

			tmp_cl= ( *cube_list_pointer );
			tmp_pl = *_tmp_pending_list_pointer;
			tmp_nl = *_tmp_node_list_pointer;
			tmp_va = *vertex_list_pointer;
			tmp_na = *normal_list_pointer;

			unsigned long vn_count=1;
			_tmp_pending_count=0;
			unsigned long i,j;
			unsigned __int64 actualboxid;
			for ( i=0;i<cube_count;++i )
			{
				actualboxid=tmp_cl[i].cell_id;
				if ( ( 3*vn_count + 360 ) >= vn_arraylength )
				{
					resize_2float ( vertex_list_pointer,normal_list_pointer,vn_count,vn_arraylength, ( unsigned long ) ( vn_arraylength*1.2+360 ) );
					tmp_va=*vertex_list_pointer;
					tmp_na=*normal_list_pointer;
				}

#include "autogenerated.h"

				/*			if(((tmp_cl[i].ecken >> m111)&1) && (!((tmp_cl[i].ecken >> m110)&1))) {//kante m111(1)-m110(0)
								tmp_cl[i].vn_id[pxpy]=vn_count;
								Vektor relevante_ecke = tmp_nl[tmp_cl[i].nodelist_id[m111]].coord;
								tmp_va[3 * vn_count    ] = relevante_ecke.x(); //x
								tmp_va[3 * vn_count + 1] = relevante_ecke.y(); //y
								tmp_va[3 * vn_count + 2] = relevante_ecke.z() - _cellsize * dichte2size[tmp_nl[tmp_cl[i].nodelist_id[m111]].count]; //z

								tmp_na[3 * vn_count    ] =
									0.20000f*( ((tmp_cl[i].ecken & m010bit ) >> m010)+((tmp_cl[i].ecken & m011bit ) >> m011))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m001bit ) >> m001)); //x
								tmp_na[3 * vn_count + 1] =
										 0.20000f*( ((tmp_cl[i].ecken & m100bit ) >> m100)+((tmp_cl[i].ecken & m101bit ) >> m101))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m001bit ) >> m001)); //y
								tmp_na[3 * vn_count + 2] =
										 1.00000f*( (float)((tmp_cl[i].ecken & m110bit ) >> m110)-(float)((tmp_cl[i].ecken & m111bit ) >> m111))
										+0.20000f*( (float)((tmp_cl[i].ecken & m100bit ) >> m100)-(float)((tmp_cl[i].ecken & m101bit ) >> m101)+
													(float)((tmp_cl[i].ecken & m010bit ) >> m010)-(float)((tmp_cl[i].ecken & m011bit ) >> m011))
										+0.11111f*( (float)((tmp_cl[i].ecken & m000bit ) >> m000)-(float)((tmp_cl[i].ecken & m001bit ) >> m001)); //z

								vn_count++;
							} else if((!((tmp_cl[i].ecken >> m111)&1)) && ((tmp_cl[i].ecken >> m110)&1)) {//kante m110(1)-m111(0)
								tmp_cl[i].vn_id[pxpy]=vn_count;
								Vektor relevante_ecke = tmp_nl[tmp_cl[i].nodelist_id[m110]].coord;
								tmp_va[3 * vn_count    ] = relevante_ecke.x(); //x
								tmp_va[3 * vn_count + 1] = relevante_ecke.y(); //y
								tmp_va[3 * vn_count + 2] = relevante_ecke.z() + _cellsize * dichte2size[tmp_nl[tmp_cl[i].nodelist_id[m110]].count]; //z

								tmp_na[3 * vn_count    ] =
										 0.20000f*( ((tmp_cl[i].ecken & m010bit ) >> m010)+((tmp_cl[i].ecken & m011bit ) >> m011))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m001bit ) >> m001)); //x
								tmp_na[3 * vn_count + 1] =
										 0.20000f*( ((tmp_cl[i].ecken & m100bit ) >> m100)+((tmp_cl[i].ecken & m101bit ) >> m101))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m001bit ) >> m001)); //y
								tmp_na[3 * vn_count + 2] =
										 1.00000f*( (float)((tmp_cl[i].ecken & m110bit ) >> m110)-(float)((tmp_cl[i].ecken & m111bit ) >> m111))
										+0.20000f*( (float)((tmp_cl[i].ecken & m100bit ) >> m100)-(float)((tmp_cl[i].ecken & m101bit ) >> m101)+
													(float)((tmp_cl[i].ecken & m010bit ) >> m010)-(float)((tmp_cl[i].ecken & m011bit ) >> m011))
										+0.11111f*( (float)((tmp_cl[i].ecken & m000bit ) >> m000)-(float)((tmp_cl[i].ecken & m001bit ) >> m001)); //z

								vn_count++;
							}
							if(((tmp_cl[i].ecken >> m111)&1) && (!((tmp_cl[i].ecken >> m101)&1))) {//kante m111(1)-m101(0)
								tmp_cl[i].vn_id[pxpz]=vn_count;
								Vektor relevante_ecke = tmp_nl[tmp_cl[i].nodelist_id[m111]].coord;
								tmp_va[3 * vn_count    ] = relevante_ecke.x(); //x
								tmp_va[3 * vn_count + 1] = relevante_ecke.y() - _cellsize * dichte2size[tmp_nl[tmp_cl[i].nodelist_id[m111]].count]; //y
								tmp_va[3 * vn_count + 2] = relevante_ecke.z(); //z

								tmp_na[3 * vn_count    ] =
										 0.20000f*( ((tmp_cl[i].ecken & m001bit ) >> m001)+((tmp_cl[i].ecken & m011bit ) >> m011))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m010bit ) >> m010)); //x
								tmp_na[3 * vn_count + 1] =
										 1.00000f*( (float)((tmp_cl[i].ecken & m101bit ) >> m101)-(float)((tmp_cl[i].ecken & m111bit ) >> m111))
										+0.20000f*( (float)((tmp_cl[i].ecken & m100bit ) >> m100)-(float)((tmp_cl[i].ecken & m110bit ) >> m110)+
													(float)((tmp_cl[i].ecken & m001bit ) >> m001)-(float)((tmp_cl[i].ecken & m011bit ) >> m011))
										+0.11111f*( (float)((tmp_cl[i].ecken & m000bit ) >> m000)-(float)((tmp_cl[i].ecken & m010bit ) >> m010)); //y
								tmp_na[3 * vn_count + 2] =
										 0.20000f*( ((tmp_cl[i].ecken & m100bit ) >> m100)+((tmp_cl[i].ecken & m110bit ) >> m110))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m010bit ) >> m010)); //z

								vn_count++;
							} else if((!((tmp_cl[i].ecken >> m111)&1)) && ((tmp_cl[i].ecken >> m101)&1)) {//kante m111(0)-m101(1)
								tmp_cl[i].vn_id[pxpz]=vn_count;
								Vektor relevante_ecke = tmp_nl[tmp_cl[i].nodelist_id[m101]].coord;
								tmp_va[3 * vn_count    ] = relevante_ecke.x(); //x
								tmp_va[3 * vn_count + 1] = relevante_ecke.y() + _cellsize * dichte2size[tmp_nl[tmp_cl[i].nodelist_id[m101]].count]; //y
								tmp_va[3 * vn_count + 2] = relevante_ecke.z(); //z

								tmp_na[3 * vn_count    ] =
										 0.20000f*( ((tmp_cl[i].ecken & m001bit ) >> m001)+((tmp_cl[i].ecken & m011bit ) >> m011))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m010bit ) >> m010)); //x
								tmp_na[3 * vn_count + 1] =
										 1.00000f*( (float)((tmp_cl[i].ecken & m101bit ) >> m101)-(float)((tmp_cl[i].ecken & m111bit ) >> m111))
										+0.20000f*( (float)((tmp_cl[i].ecken & m100bit ) >> m100)-(float)((tmp_cl[i].ecken & m110bit ) >> m110)+
													(float)((tmp_cl[i].ecken & m001bit ) >> m001)-(float)((tmp_cl[i].ecken & m011bit ) >> m011))
										+0.11111f*( (float)((tmp_cl[i].ecken & m000bit ) >> m000)-(float)((tmp_cl[i].ecken & m010bit ) >> m010)); //y
								tmp_na[3 * vn_count + 2] =
										 0.20000f*( ((tmp_cl[i].ecken & m100bit ) >> m100)+((tmp_cl[i].ecken & m110bit ) >> m110))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m010bit ) >> m010)); //z

								vn_count++;
							}
							if(((tmp_cl[i].ecken >> m111)&1) && (!((tmp_cl[i].ecken >> m011)&1))) {//kante m111(1)-m011(0)
								tmp_cl[i].vn_id[pypz]=vn_count;
								Vektor relevante_ecke = tmp_nl[tmp_cl[i].nodelist_id[m111]].coord;
								tmp_va[3 * vn_count    ] = relevante_ecke.x() - _cellsize * dichte2size[tmp_nl[tmp_cl[i].nodelist_id[m111]].count]; //x
								tmp_va[3 * vn_count + 1] = relevante_ecke.y(); //y
								tmp_va[3 * vn_count + 2] = relevante_ecke.z(); //z

								tmp_na[3 * vn_count    ] =
										 1.00000f*( (float)((tmp_cl[i].ecken & m011bit ) >> m011)-(float)((tmp_cl[i].ecken & m111bit ) >> m111))
										+0.20000f*( (float)((tmp_cl[i].ecken & m001bit ) >> m001)-(float)((tmp_cl[i].ecken & m101bit ) >> m101)+
													(float)((tmp_cl[i].ecken & m010bit ) >> m010)-(float)((tmp_cl[i].ecken & m110bit ) >> m110))
										+0.11111f*( (float)((tmp_cl[i].ecken & m000bit ) >> m000)-(float)((tmp_cl[i].ecken & m100bit ) >> m100)); //x
								tmp_na[3 * vn_count + 1] =
										 0.20000f*( ((tmp_cl[i].ecken & m001bit ) >> m001)+((tmp_cl[i].ecken & m101bit ) >> m101))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m100bit ) >> m100)); //y
								tmp_na[3 * vn_count + 2] =
										 0.20000f*( ((tmp_cl[i].ecken & m010bit ) >> m010)+((tmp_cl[i].ecken & m110bit ) >> m110))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m100bit ) >> m100)); //z

								vn_count++;
							} else if((!((tmp_cl[i].ecken >> m111)&1)) && ((tmp_cl[i].ecken >> m011)&1)) {//kante m111(0)-m011(1)
								tmp_cl[i].vn_id[pypz]=vn_count;
								Vektor relevante_ecke = tmp_nl[tmp_cl[i].nodelist_id[m011]].coord;
								tmp_va[3 * vn_count    ] = relevante_ecke.x() + _cellsize * dichte2size[tmp_nl[tmp_cl[i].nodelist_id[m011]].count]; //x
								tmp_va[3 * vn_count + 1] = relevante_ecke.y(); //y
								tmp_va[3 * vn_count + 2] = relevante_ecke.z(); //z

								tmp_na[3 * vn_count    ] =
										 1.00000f*( (float)((tmp_cl[i].ecken & m011bit ) >> m011)-(float)((tmp_cl[i].ecken & m111bit ) >> m111))
										+0.20000f*( (float)((tmp_cl[i].ecken & m001bit ) >> m001)-(float)((tmp_cl[i].ecken & m101bit ) >> m101)+
													(float)((tmp_cl[i].ecken & m010bit ) >> m010)-(float)((tmp_cl[i].ecken & m110bit ) >> m110))
										+0.11111f*( (float)((tmp_cl[i].ecken & m000bit ) >> m000)-(float)((tmp_cl[i].ecken & m100bit ) >> m100)); //x
								tmp_na[3 * vn_count + 1] =
										 0.20000f*( ((tmp_cl[i].ecken & m001bit ) >> m001)+((tmp_cl[i].ecken & m101bit ) >> m101))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m100bit ) >> m100)); //y
								tmp_na[3 * vn_count + 2] =
										 0.20000f*( ((tmp_cl[i].ecken & m010bit ) >> m010)+((tmp_cl[i].ecken & m110bit ) >> m110))
										+0.11111f*( ((tmp_cl[i].ecken & m000bit ) >> m000)+((tmp_cl[i].ecken & m100bit ) >> m100)); //z

								vn_count++;
							}*/

				for ( j=0;j<_tmp_pending_count;j++ )
				{
					unsigned __int64 tmp=actualboxid - tmp_cl[tmp_pl[j]].cell_id;
					if ( tmp & ( 0xffffffffffffffffULL ^ _dxyz ) )
					{
						if ( ( tmp_cl[tmp_pl[j]].cell_id + _dxyz ) < actualboxid )
						{
							//fr die nachbarn ins innere die normalenbeitr?e einbringen!
							tmp_pl[j--]=tmp_pl[--_tmp_pending_count];//bei geschachtelten if: _tmp_pending[j--]=_tmp_pending[_tmp_pending_count--];
						}
					}
					else
					{
						if ( tmp == _dx )
						{
							if ( tmp_cl[tmp_pl[j]].vn_id[pxmy] )
							{
								tmp_cl[i].vn_id[mxmy]=tmp_cl[tmp_pl[j]].vn_id[pxmy];
								tmp_na[3 * tmp_cl[i].vn_id[mxmy]    ] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mxmy] + 1] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[mxmy] + 2] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //z
							}
							if ( tmp_cl[tmp_pl[j]].vn_id[pxmz] )
							{
								tmp_cl[i].vn_id[mxmz]=tmp_cl[tmp_pl[j]].vn_id[pxmz];
								tmp_na[3 * tmp_cl[i].vn_id[mxmz]    ] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mxmz] + 1] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[mxmz] + 2] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //z
							}
							if ( tmp_cl[tmp_pl[j]].vn_id[pxpy] )
							{
								tmp_cl[i].vn_id[mxpy]=tmp_cl[tmp_pl[j]].vn_id[pxpy];
								tmp_na[3 * tmp_cl[i].vn_id[mxpy]    ] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mxpy] + 1] +=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m000bit ) >> m000 ) + ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[mxpy] + 2] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) ); //z
							}
							if ( tmp_cl[tmp_pl[j]].vn_id[pxpz] )
							{
								tmp_cl[i].vn_id[mxpz]=tmp_cl[tmp_pl[j]].vn_id[pxpz];
								tmp_na[3 * tmp_cl[i].vn_id[mxpz]    ] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mxpz] + 1] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[mxpz] + 2] +=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m000bit ) >> m000 ) + ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) ); //z
							}
						}
						else if ( tmp == _dy )
						{
							if ( tmp_cl[tmp_pl[j]].vn_id[mxpy] )
							{
								tmp_cl[i].vn_id[mxmy]=tmp_cl[tmp_pl[j]].vn_id[mxpy];
								tmp_na[3 * tmp_cl[i].vn_id[mxmy]    ] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mxmy] + 1] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[mxmy] + 2] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //z
							}
							if ( tmp_cl[tmp_pl[j]].vn_id[pxpy] )
							{
								tmp_cl[i].vn_id[pxmy]=tmp_cl[tmp_pl[j]].vn_id[pxpy];
								tmp_na[3 * tmp_cl[i].vn_id[pxmy]    ] +=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m000bit ) >> m000 ) + ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[pxmy] + 1] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[pxmy] + 2] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) ); //z
							}
							if ( tmp_cl[tmp_pl[j]].vn_id[pymz] )
							{
								tmp_cl[i].vn_id[mymz]=tmp_cl[tmp_pl[j]].vn_id[pymz];
								tmp_na[3 * tmp_cl[i].vn_id[mymz]    ] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[mymz] + 1] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mymz] + 2] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
							}
							if ( tmp_cl[tmp_pl[j]].vn_id[pypz] )
							{
								tmp_cl[i].vn_id[mypz]=tmp_cl[tmp_pl[j]].vn_id[pypz];
								tmp_na[3 * tmp_cl[i].vn_id[mypz]    ] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) ); //z
								tmp_na[3 * tmp_cl[i].vn_id[mypz] + 1] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mypz] + 2] +=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m000bit ) >> m000 ) + ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) ); //y
							}
						}
						else if ( tmp == _dz )
						{
							if ( tmp_cl[tmp_pl[j]].vn_id[mxpz] )
							{
								tmp_cl[i].vn_id[mxmz]=tmp_cl[tmp_pl[j]].vn_id[mxpz];
								tmp_na[3 * tmp_cl[i].vn_id[mxmz]    ] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mxmz] + 1] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[mxmz] + 2] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //z
							}
							if ( tmp_cl[tmp_pl[j]].vn_id[pxpz] )
							{
								tmp_cl[i].vn_id[pxmz]=tmp_cl[tmp_pl[j]].vn_id[pxpz];
								tmp_na[3 * tmp_cl[i].vn_id[pxmz]    ] +=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m000bit ) >> m000 ) + ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[pxmz] + 1] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) ); //z
								tmp_na[3 * tmp_cl[i].vn_id[pxmz] + 2] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) ); //y
							}
							if ( tmp_cl[tmp_pl[j]].vn_id[mypz] )
							{
								tmp_cl[i].vn_id[mymz]=tmp_cl[tmp_pl[j]].vn_id[mypz];
								tmp_na[3 * tmp_cl[i].vn_id[mymz]    ] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //z
								tmp_na[3 * tmp_cl[i].vn_id[mymz] + 1] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mymz] + 2] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
							}
							if ( tmp_cl[tmp_pl[j]].vn_id[pypz] )
							{
								tmp_cl[i].vn_id[pymz]=tmp_cl[tmp_pl[j]].vn_id[pypz];
								tmp_na[3 * tmp_cl[i].vn_id[pymz]    ] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) ); //z
								tmp_na[3 * tmp_cl[i].vn_id[pymz] + 1] +=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m000bit ) >> m000 ) + ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[pymz] + 2] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) ); //y
							}
						}
						else if ( tmp == _dxy )
						{
							if ( tmp_cl[tmp_pl[j]].vn_id[pxpy] )
							{
								tmp_cl[i].vn_id[mxmy]=tmp_cl[tmp_pl[j]].vn_id[pxpy];
								tmp_na[3 * tmp_cl[i].vn_id[mxmy]    ] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mxmy] + 1] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[mxmy] + 2] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //z
							}
						}
						else if ( tmp == _dxz )
						{
							if ( tmp_cl[tmp_pl[j]].vn_id[pxpz] )
							{
								tmp_cl[i].vn_id[mxmz]=tmp_cl[tmp_pl[j]].vn_id[pxpz];
								tmp_na[3 * tmp_cl[i].vn_id[mxmz]    ] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mxmz] + 1] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
								tmp_na[3 * tmp_cl[i].vn_id[mxmz] + 2] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //z
							}
						}
						else if ( tmp == _dyz )
						{
							if ( tmp_cl[tmp_pl[j]].vn_id[pypz] )
							{
								tmp_cl[i].vn_id[mymz]=tmp_cl[tmp_pl[j]].vn_id[pypz];
								tmp_na[3 * tmp_cl[i].vn_id[mymz]    ] +=
								    1.00000f* ( ( float ) ( ( tmp_cl[i].ecken & m000bit ) >> m000 )- ( float ) ( ( tmp_cl[i].ecken & m100bit ) >> m100 ) )
								    +0.20000f* ( ( float ) ( ( tmp_cl[i].ecken & m010bit ) >> m010 )- ( float ) ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) +
								                 ( float ) ( ( tmp_cl[i].ecken & m001bit ) >> m001 )- ( float ) ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.11111f* ( ( float ) ( ( tmp_cl[i].ecken & m011bit ) >> m011 )- ( float ) ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //z
								tmp_na[3 * tmp_cl[i].vn_id[mymz] + 1] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m010bit ) >> m010 ) + ( ( tmp_cl[i].ecken & m110bit ) >> m110 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //x
								tmp_na[3 * tmp_cl[i].vn_id[mymz] + 2] -=
								    0.20000f* ( ( ( tmp_cl[i].ecken & m001bit ) >> m001 ) + ( ( tmp_cl[i].ecken & m101bit ) >> m101 ) )
								    +0.11111f* ( ( ( tmp_cl[i].ecken & m011bit ) >> m011 ) + ( ( tmp_cl[i].ecken & m111bit ) >> m111 ) ); //y
							}
						}
					}
				}
				if ( _tmp_pending_count+10 >= _tmp_max_pending_count )
				{
					resize_unsigned_long_list ( _tmp_pending_list_pointer,_tmp_pending_count,_tmp_max_pending_count, ( unsigned long ) ( _tmp_max_pending_count*1.2+10 ) );
					tmp_pl=*_tmp_pending_list_pointer;
				}
				tmp_pl[_tmp_pending_count++]=i;//new cube to _tmp_pending-list
			}
			float r;
			for ( i=0;i<vn_count;i++ )
			{
				r=/*g_normal_length/*/sqrtf ( tmp_na[3*i]*tmp_na[3*i]+tmp_na[3*i+1]*tmp_na[3*i+1]+tmp_na[3*i+2]*tmp_na[3*i+2] );
				tmp_na[3*i]*=r;
				tmp_na[3*i+1]*=r;
				tmp_na[3*i+2]*=r;

			}
			return vn_count;
		}
		/*	TODO:int destroyparticle(const unsigned long i) {
		        if (i>_particlecount - 1) {
		            return 0;
		        } else {
		            _particlecount--;

					_particle[i]=_particle[_particlecount];

					_particle[_particle[_particlecount].get_sortlistindex(0)].set_sortlistindex(0,get_sortlistindex(_particlecount,0));
		            _particle[_particle[_particlecount].get_sortlistindex(1)].set_sortlistindex(1,get_sortlistindex(_particlecount,1));
		            _particle[_particle[_particlecount].get_sortlistindex(2)].set_sortlistindex(2,get_sortlistindex(_particlecount,2));
		            _particle[_particle[_particlecount].get_sortlistindex(3)].set_sortlistindex(3,get_sortlistindex(_particlecount,3));

					_particle[sortlist2particleid(_sortlist[_particle[i].get_sortlistindex(0)])].set_sortlistindex(0,_particle[i].get_sortlistindex(0));
					_particle[sortlist2particleid(_sortlist[_particle[i].get_sortlistindex(1)])].set_sortlistindex(1,_particle[i].get_sortlistindex(1));
					_particle[sortlist2particleid(_sortlist[_particle[i].get_sortlistindex(2)])].set_sortlistindex(2,_particle[i].get_sortlistindex(2));
					_particle[sortlist2particleid(_sortlist[_particle[i].get_sortlistindex(3)])].set_sortlistindex(3,_particle[i].get_sortlistindex(3));

		            return 1;
		        }
		    }*/
		unsigned __int64 get_sortlistvalue ( const unsigned short listid, const unsigned long id )
		{
			return _sortlist[listid*_maxparticlecount + id];
		}
		unsigned __int64 cell6bitplus ( const float x, const float y, const float z )  	//3bit higher resolution to get the 4(8) staggered grids at once
		{
			return	( ( ( ( unsigned __int64 ) ( _cellsz* ( ( z-_minz ) /_lenz ) *4 ) ) << _cellsxplus2ybits ) << 4 ) +
			       ( ( ( ( unsigned __int64 ) ( _cellsy* ( ( y-_miny ) /_leny ) *4 ) ) << _cellsxplus2bits ) ) +
			       ( unsigned __int64 ) ( _cellsx* ( ( x-_minx ) /_lenx ) *4 );
		}
		unsigned __int64 cell6bitplus ( Vektor x )
		{
			return cell6bitplus ( x.x(), x.y(), x.z() );
		}
		unsigned __int64 cellx2bitplus ( Vektor x )
		{
			return	( unsigned __int64 ) ( _cellsx* ( ( x.x()-_minx ) /_lenx ) *4 );
		}
		unsigned __int64 celly2bitplus ( Vektor x )
		{
			return	( unsigned __int64 ) ( _cellsy* ( ( x.y()-_miny ) /_leny ) *4 );
		}
		unsigned __int64 cellz2bitplus ( Vektor x )
		{
			return	( unsigned __int64 ) ( _cellsz* ( ( x.z()-_minz ) /_lenz ) *4 );
		}
		unsigned char get_relevantsortlist ( const unsigned long i )
		{
			return _particle[i].relevantsortlist();
		}
		unsigned long get_sortlistindex ( const unsigned long particle_id, const unsigned char sortlist_id )
		{
			return _particle[particle_id].get_sortlistindex ( sortlist_id );
		}
		void collidewithneighbours ( const unsigned long id )
		{

			unsigned char activesortlist = get_relevantsortlist ( id );
			unsigned long startsortlistid = get_sortlistindex ( id,activesortlist );
			unsigned __int64 startcellid = sortlist2cellid ( _sortlist[startsortlistid] );

			unsigned __int64 sortlistbase = startsortlistid - ( startsortlistid % _maxparticlecount );
			unsigned __int64 sortlistindex = startsortlistid - sortlistbase + 1;

			while ( ( ( sortlist2cellid ( _sortlist[sortlistbase + sortlistindex] ) - startcellid ) <= 1 ) && ( sortlistindex < _particlecount ) )
			{
				if ( condcollide ( id,sortlist2particleid ( _sortlist[sortlistbase + sortlistindex] ) ) == 1 )
					return;
				++sortlistindex;
			}
			sortlistindex=startsortlistid - sortlistbase - 1;
			while ( ( ( startcellid - sortlist2cellid ( _sortlist[sortlistbase + sortlistindex] ) ) <= 1 ) && ( sortlistindex < _particlecount ) )
			{
				if ( condcollide ( id,sortlist2particleid ( _sortlist[sortlistbase + sortlistindex] ) ) == 1 )
					return;
				--sortlistindex;
			}
		}
		int condcollide ( const unsigned long first, const unsigned long second )
		{
			if ( ( first < second ) & ( second < _particlecount ) )   //TODO: workaround?? exit outer loop
			{
				return collide ( first,second );
			}
			return 10;
		}
		int collide ( const unsigned long first, const unsigned long second )
		{
			Vektor f; //force
			Vektor r;
			float offset;
			r = _particle[first].x()-_particle[second].x(); //connecting Vektor
			offset = r*_particle[second].a();
			if ( second >=_particlecount_physik )
			{
				switch ( _particle[second].kind() )
				{
					case ( moving ) :
									cout << "MOVING particle with wrong id " << second << "!" << endl;
						exit ( 9856 );
						return 1;
						break;
					case ( boundary ) ://flip velocity vertical to the boundary
									if ( offset < 0 )
							{
								_particle[first].setx ( _particle[first].x() - _particle[second].v().norm ( offset ) );
								_particle[first].setv ( _particle[first].v() + _particle[second].v() * ( ( _particle[first].v() *_particle[second].v().norm ( -1.01f ) ) ) );
							}
						return 10;
						break;
					case ( teleport ) :
									if ( offset < 0 )
							{
								_particle[first].setx ( _particle[second].v() + ( _particle[first].x()-_particle[second].x() ) /10.0f );
								return 1;
							}
						break;
					case ( shift ) :
									if ( offset < 0 )
							{
								_particle[first].setx ( _particle[first].x() +_particle[second].v() );
								return 1;
							}
						break;
					case ( setspeed ) :
									if ( offset < 0 )
							{
								_particle[first].setv ( _particle[second].v() );
								return 10;
							}
						break;
						//case DELETE:
						//	destroyparticle(first);
						//	break;
					default:
						break;
				}
			}
			else
			{
				float d = r.abs(); //distance
				if ( d < _particlesize )  //collide
				{
#define _gaskonstante 1.05
#define _rho0 0.100
#define _m 1
					f=r.norm ( _gaskonstante* ( _particle[first].rho() +_particle[second].rho()-2*_rho0 ) /2/_particle[second].rho() *w_poly6_grad ( d ) /*-5*/ ) + ( _particle[first].v()-_particle[second].v() ) * ( 14.72f );
					_particle[first].push ( f* ( -1 ) );
					_particle[second].push ( f );
				}
			}
			return 10;
		}
		void collidewithneighboursforrho ( const unsigned long id )
		{
			unsigned char activesortlist = get_relevantsortlist ( id );
			unsigned long startsortlistid = get_sortlistindex ( id,activesortlist );
			unsigned __int64 startcellid = sortlist2cellid ( _sortlist[startsortlistid] );

			unsigned __int64 sortlistbase = startsortlistid - ( startsortlistid % _maxparticlecount );
			unsigned __int64 sortlistindex = startsortlistid - sortlistbase + 1;

			while ( ( ( sortlist2cellid ( _sortlist[sortlistbase + sortlistindex] ) - startcellid ) <= 1 ) && ( sortlistindex < _particlecount ) )
			{
				condcollideforrho ( id,sortlist2particleid ( _sortlist[sortlistbase + sortlistindex] ) );
				++sortlistindex;
			}
			sortlistindex=startsortlistid - sortlistbase - 1;
			while ( ( ( startcellid - sortlist2cellid ( _sortlist[sortlistbase + sortlistindex] ) ) <= 1 ) && ( sortlistindex < _particlecount ) )
			{
				condcollideforrho ( id,sortlist2particleid ( _sortlist[sortlistbase + sortlistindex] ) );
				--sortlistindex;
			}
		}
		void condcollideforrho ( const unsigned long first, const unsigned long second )
		{
			if ( ( first < second ) & ( second < _particlecount_physik ) )   //TODO: workaround?? exit outer loop
			{
				collideforrho ( first,second );
			}
		}
		void collideforrho ( const unsigned long first, const unsigned long second )
		{
			Vektor r = _particle[first].x()-_particle[second].x(); //connecting Vektor
			float dq = r.absabs(); //distance square
			if ( dq < _particlesize_q )  //collide
			{
				float rho=_m*w_poly6 ( dq );
				_particle[first].addrho ( rho );
				_particle[second].addrho ( rho );
			}
		}
		unsigned __int64 sortlist2cellid ( unsigned __int64 in )
		{
			return ( in >> _maxparticlecountbits );
		}
		unsigned __int64 sortlist2cellxid ( unsigned __int64 in )
		{
			return ( ( in >> _maxparticlecountbits ) & ( ( 0xffffffffffffffffULL ) >> ( 64-_cellsxplus2bits ) ) );
		}
		unsigned __int64 sortlist2cellyid ( unsigned __int64 in )
		{
			return ( ( in >> _cellsxplus2particlebits ) & ( ( 0xffffffffffffffffULL ) >> ( 64-_cellsybits ) ) );
		}
		unsigned __int64 sortlist2cellzid ( unsigned __int64 in )
		{
			return ( in >> _cellsxplus2yparticlebits );
		}
		void putparticles2cells()   //only workaround. to be removed. putparticle2cell brings in invalid particle-indices
		{
			unsigned long i;
			unsigned __int64 tmpx, tmpy, tmpz;
			for ( i=0; i<_particlecount;i++ )
			{
				//TODO:!!!!!!!		for(i=0; i<_particlecount_physik;i++) {
				tmpx = cellx2bitplus ( _particle[i].x() );
				tmpy = celly2bitplus ( _particle[i].x() );
				tmpz = cellz2bitplus ( _particle[i].x() );



				_sortlist[                    i] = ( ( ( ( tmpz    >> 2 ) << _cellsxplus2ybits ) | ( ( tmpy      >> 2 ) << _cellsxplus2bits ) | tmpx ) <<_maxparticlecountbits ) | i;
				_sortlist[  _maxparticlecount+i] = ( ( ( ( ( tmpz+2 ) >> 2 ) << _cellsxplus2ybits ) | ( ( tmpy      >> 2 ) << _cellsxplus2bits ) | tmpx ) <<_maxparticlecountbits ) | i;
				_sortlist[2*_maxparticlecount+i] = ( ( ( ( tmpz    >> 2 ) << _cellsxplus2ybits ) | ( ( ( tmpy + 2 ) >> 2 ) << _cellsxplus2bits ) | tmpx ) <<_maxparticlecountbits ) | i;
				_sortlist[3*_maxparticlecount+i] = ( ( ( ( ( tmpz+2 ) >> 2 ) << _cellsxplus2ybits ) | ( ( ( tmpy + 2 ) >> 2 ) << _cellsxplus2bits ) | tmpx ) <<_maxparticlecountbits ) | i;

				_particle[i].set_relevantsortlist ( ( ( 1- ( unsigned char ) ( ( ( tmpy ) & 1 ) ^ ( ( ( tmpy ) & 2 ) >> 1 ) ) ) << 1 ) | ( 1- ( unsigned char ) ( ( tmpz & 1 ) ^ ( ( tmpz & 2 ) >> 1 ) ) ) );//yz?zy
			}
			set_unsorted();
		}
		void putparticle2cell ( const unsigned long particleindex )
		{
			unsigned __int64 tmpx = cellx2bitplus ( _particle[particleindex].x() );
			unsigned __int64 tmpy = celly2bitplus ( _particle[particleindex].x() );
			unsigned __int64 tmpz = cellz2bitplus ( _particle[particleindex].x() );

			_sortlist[_particle[particleindex].get_sortlistindex ( 0 ) ] = ( ( ( ( tmpz     >> 2 ) << _cellsxplus2ybits ) | ( ( tmpy      >> 2 ) << _cellsxplus2bits ) | tmpx ) <<_maxparticlecountbits ) | particleindex;
			_sortlist[_particle[particleindex].get_sortlistindex ( 1 ) ] = ( ( ( ( ( tmpz+2 ) >> 2 ) << _cellsxplus2ybits ) | ( ( tmpy      >> 2 ) << _cellsxplus2bits ) | tmpx ) <<_maxparticlecountbits ) | particleindex;
			_sortlist[_particle[particleindex].get_sortlistindex ( 2 ) ] = ( ( ( ( tmpz     >> 2 ) << _cellsxplus2ybits ) | ( ( ( tmpy + 2 ) >> 2 ) << _cellsxplus2bits ) | tmpx ) <<_maxparticlecountbits ) | particleindex;
			_sortlist[_particle[particleindex].get_sortlistindex ( 3 ) ] = ( ( ( ( ( tmpz+2 ) >> 2 ) << _cellsxplus2ybits ) | ( ( ( tmpy + 2 ) >> 2 ) << _cellsxplus2bits ) | tmpx ) <<_maxparticlecountbits ) | particleindex;
			_particle[particleindex].set_relevantsortlist ( ( ( 1- ( unsigned char ) ( ( ( tmpy ) & 1 ) ^ ( ( ( tmpy ) & 2 ) >> 1 ) ) ) << 1 ) | ( 1- ( unsigned char ) ( ( tmpz & 1 ) ^ ( ( tmpz & 2 ) >> 1 ) ) ) );//yz?zy
		}
		void quickbubblesort ( const unsigned long lo, const unsigned long hi )
		{
			unsigned long i = lo;
			unsigned long j = hi;
			unsigned __int64 x = _sortlist[ ( lo+hi ) /2];

			while ( i + 2000 <= j + 2000 )
			{
				while ( _sortlist[i] < x )
					i++;
				while ( _sortlist[j] > x )
					j--;
				if ( i + 2000 <= j + 2000 )
				{
					exchange ( i++, j-- );
				}
			}

			if ( j - lo + 2000 > 2000 )
			{
				if ( j - lo > 10 )
				{
					quickbubblesort ( lo, j );
				}
				else
				{
					bubblesort ( lo, j );
				}
			}
			if ( hi - i  + 2000 > 2000 )
			{
				if ( hi - i > 10 )
				{
					quickbubblesort ( i, hi );
				}
				else
				{
					bubblesort ( i, hi );
				}
			}
		}
		void bubblesort ( unsigned long lo, unsigned long hi )  //wird von quickbubblesort gebraucht
		{
			unsigned long i = lo, tmp;
			while ( lo < hi )
			{
				tmp = lo;
				while ( i<hi )
				{
					if ( _sortlist[i] > _sortlist[i+1] )
					{
						exchange ( i, i+1 );
						tmp = i;
					}
					i++;
				}
				hi=tmp;
				while ( i>lo )
				{
					if ( _sortlist[i] < _sortlist[i-1] )
					{
						exchange ( i, i-1 );
						tmp = i;
					}
					i--;
				}
				lo=tmp;
			}
		}
		void exchange ( const unsigned long i, const unsigned long j )
		{
			unsigned __int64 t = _sortlist[i];
			_sortlist[i] = _sortlist[j];
			_sortlist[j] = t;
		}
		void sortparticlelists()
		{
			quickbubblesort ( 0                ,                     _particlecount - 1 );
			quickbubblesort ( _maxparticlecount,     _maxparticlecount + _particlecount - 1 );
			quickbubblesort ( 2 * _maxparticlecount, 2 * _maxparticlecount + _particlecount - 1 );
			quickbubblesort ( 3 * _maxparticlecount, 3 * _maxparticlecount + _particlecount - 1 );
			fixsortlistindices();
			_listssorted = true;
		}
		void fixsortlistindices()
		{
			unsigned long i;
			for ( i=0; i<_particlecount; i++ )
			{
				_particle[sortlist2particleid ( _sortlist[i                    ] ) ].setsortlistindex ( 0,i );
				_particle[sortlist2particleid ( _sortlist[i+  _maxparticlecount] ) ].setsortlistindex ( 1,i+  _maxparticlecount );
				_particle[sortlist2particleid ( _sortlist[i+2*_maxparticlecount] ) ].setsortlistindex ( 2,i+2*_maxparticlecount );
				_particle[sortlist2particleid ( _sortlist[i+3*_maxparticlecount] ) ].setsortlistindex ( 3,i+3*_maxparticlecount );
			}
		}
		unsigned long sortlist2particleid ( unsigned __int64 in )
		{
			return ( unsigned long ) ( in & _ParticleIdBitmask );
		}
		void set_unsorted()
		{
			_listssorted = false;
		}
		void mergeparticles ( const unsigned long first, const unsigned long second )
		{
//TODO: implement
//		float m1 = _particle[first].m();
//		float m2 = _particle[second].m();
//		float m3 = m1 + m2;
//		_particle[first].setm(m3);
//		_particle[first].setx((_particle[first].x()*m1 + _particle[second].x()*m2)/m3);
//		_particle[first].setv((_particle[first].v()*m1 + _particle[second].v()*m2)/m3);
////		_particle[first].sett((_particle[first].t()*m1 + _particle[second].t()*m2)/m3);
//		putparticle2cell(first);
//		destroyparticle(second);
		}
		Vektor sortlistindex2coord ( unsigned __int64 sortlistindex )
		{
			return Vektor (
			           ( float ) ( ( sortlistindex >> ( _maxparticlecountbits+2 ) ) & ( 0xffffffffffffffffULL >> ( 66 - _cellsxplus2bits ) ) ) *_cellsize + _minx,
			           ( float ) ( ( sortlistindex >> _cellsxplus2particlebits ) & ( 0xffffffffffffffffULL >> ( 64 - _cellsybits ) ) ) *_cellsize + _miny,
			           ( float ) ( sortlistindex >> _cellsxplus2yparticlebits ) *_cellsize + _minz );
		}
		void setsortlistvalue ( const unsigned char listid, const unsigned long id, const unsigned __int64 value )
		{
			_sortlist[listid*_maxparticlecount + id] = value;
		}
		void tetraeder2particle ( Vektor * v )
		{
			float Volumen;
			Volumen=fabs ( ( v[1]-v[0] ) % ( v[2]-v[0] ) * ( v[3]-v[0] ) /6.0f );
			int max_lq_i,max_lq_j,i,j;
			float max_lq=0;
			for ( i=0;i<4;++i )
			{
				for ( j=0;j<4;++j )
				{
					if ( ( v[i]-v[j] ).absabs() > max_lq )
					{
						max_lq_i = i;
						max_lq_j = j;
						max_lq = ( v[i]-v[j] ).absabs();
					}
				}
			}
			float Volumen0 = 0.2;//TODO: Dringend hier .... ja was?
			float ratio = Volumen/Volumen0;
			if ( ratio > 1.0002f )
			{
				Vektor v2[4];
				memcpy ( v2,v,sizeof ( Vektor ) *4 );
				float teilverhaeltnis;
				teilverhaeltnis = floorf ( ( 1.0f+ratio ) /2.0f ) /ratio;

				v[max_lq_i]=v2[max_lq_j]=v[max_lq_i]* ( 1.0f-teilverhaeltnis ) +v[max_lq_j]*teilverhaeltnis;
				tetraeder2particle ( v2 );
				tetraeder2particle ( v );
			}
			else
			{
				Vektor coord = ( v[0]+v[1]+v[2]+v[3] ) *0.25f;
				newparticle ( coord.x(),coord.y(),coord.z(),0.0001f,0.0001f,0.0001f,1,1,moving );
			}
		}
		void newparticle ( const float x, const float y, const float z, const float vx, const float vy, const float vz, const float m, const float t, const PARTICLE_TYPE kind )
		{
			if ( _particlecount_boundary > 0 )
			{
				newcontrollparticle (
				    _particle[_particlecount_physik].x().x(),
				    _particle[_particlecount_physik].x().y(),
				    _particle[_particlecount_physik].x().z(),
				    _particle[_particlecount_physik].a().x(),
				    _particle[_particlecount_physik].a().y(),
				    _particle[_particlecount_physik].a().z(),
				    _particle[_particlecount_physik].v().x(),
				    _particle[_particlecount_physik].v().y(),
				    _particle[_particlecount_physik].v().z(),
				    _particle[_particlecount_physik].kind() );
			}
			_particle[_particlecount_physik].set ( x,y,z,vx,vy,vz,_particlesize,t,kind );

			_sortlist[                    _particlecount_physik]=_particlecount_physik;
			_sortlist[_maxparticlecount  +_particlecount_physik]=_particlecount_physik;
			_sortlist[_maxparticlecount*2+_particlecount_physik]=_particlecount_physik;
			_sortlist[_maxparticlecount*3+_particlecount_physik]=_particlecount_physik;

			_particle[_particlecount_physik].set_sortlistindex ( 0,                    _particlecount_physik );
			_particle[_particlecount_physik].set_sortlistindex ( 1,_maxparticlecount  +_particlecount_physik );
			_particle[_particlecount_physik].set_sortlistindex ( 2,_maxparticlecount*2+_particlecount_physik );
			_particle[_particlecount_physik].set_sortlistindex ( 3,_maxparticlecount*3+_particlecount_physik );

			putparticle2cell ( _particlecount_physik++ );
			_particlecount++;
		}
		void newcontrollparticle ( const float x, const float y, const float z, const float bx, const float by, const float bz, const float dx, const float dy, const float dz, PARTICLE_TYPE kind )
		{
			_particle[_particlecount].set ( x,y,z,bx,by,bz,dx,dy,dz,kind );

			_sortlist[                    _particlecount]=_particlecount;
			_sortlist[_maxparticlecount  +_particlecount]=_particlecount;
			_sortlist[_maxparticlecount*2+_particlecount]=_particlecount;
			_sortlist[_maxparticlecount*3+_particlecount]=_particlecount;

			_particle[_particlecount].set_sortlistindex ( 0,                    _particlecount );
			_particle[_particlecount].set_sortlistindex ( 1,_maxparticlecount  +_particlecount );
			_particle[_particlecount].set_sortlistindex ( 2,_maxparticlecount*2+_particlecount );
			_particle[_particlecount].set_sortlistindex ( 3,_maxparticlecount*3+_particlecount );

			putparticle2cell ( _particlecount++ );
			_particlecount_boundary++;
		}
		float w_poly6 ( float & rq )
		{
			return _m_poly6*powf ( ( _particlesize_q-rq ),3 );
		}
		float w_poly6_grad ( float & r )
		{
			return -6 * _m_poly6 * powf ( ( _particlesize_q-r*r ),2 ) * r;
		}
		float w_poly6_laplace ( float & rq )
		{
			return 1;
		}

		float w_spiky ( float & r )
		{
			return _m_spiky*powf ( ( _particlesize-r ),3 );
		}
		float w_spiky_grad ( float & r )
		{
			return -3 * _m_spiky * powf ( ( _particlesize-r ),2 );
		}
		float w_spiky_laplace ( float & r )
		{
			return 1;
		}

		float w_viscosity ( float & r )
		{
			return _m_viscosity* (
			           -powf ( r/_particlesize,3 ) /2
			           +r*r/_particlesize_q
			           +_particlesize/r/2
			           -1 );
		}
		float w_viscosity_grad ( float & r )
		{
			return _m_viscosity* (
			           -1.5f*r*r/powf ( _particlesize,3 )
			           +2.0f*r/_particlesize_q
			           -_particlesize/r/r/2.0f );
		}
		float w_viscosity_laplace ( float & r )
		{
			return 1;
		}

};
//
//int fluid::edgeTable[256] = {
//	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
//	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
//	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
//	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
//	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
//	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
//	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
//	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
//	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
//	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
//	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
//	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
//	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
//	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
//	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
//	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
//	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
//	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
//	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
//	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
//	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
//	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
//	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
//	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
//	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
//	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
//	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
//	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
//	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
//	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
//	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
//	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
//};
//
//int fluid::triTable2[256][13] = {//[256][19]
//    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 3, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 0, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 3, 1, 8, 1, 9,-1,-1,-1,-1,-1,-1,-1},
//    {10, 1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 3, 0, 1, 2,10,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 0, 2, 9, 2,10,-1,-1,-1,-1,-1,-1,-1},
//    { 3, 2, 8, 2,10, 8, 8,10, 9,-1,-1,-1,-1},
//    {11, 2, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    {11, 2, 0,11, 0, 8,-1,-1,-1,-1,-1,-1,-1},
//    {11, 2, 3, 0, 1, 9,-1,-1,-1,-1,-1,-1,-1},
//    { 2, 1,11, 1, 9,11,11, 9, 8,-1,-1,-1,-1},
//    {10, 1, 3,10, 3,11,-1,-1,-1,-1,-1,-1,-1},
//    { 1, 0,10, 0, 8,10,10, 8,11,-1,-1,-1,-1},
//    { 0, 3, 9, 3,11, 9, 9,11,10,-1,-1,-1,-1},
//    { 8,10, 9, 8,11,10,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 3, 0, 4, 3, 4, 7,-1,-1,-1,-1,-1,-1,-1},
//    { 1, 9, 0, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 4, 1, 4, 7, 1, 1, 7, 3,-1,-1,-1,-1},
//    {10, 1, 2, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1},
//    { 2,10, 1, 0, 4, 7, 0, 7, 3,-1,-1,-1,-1},
//    { 4, 7, 8, 0, 2,10, 0,10, 9,-1,-1,-1,-1},
//    { 2, 7, 3, 2, 9, 7, 7, 9, 4, 2,10, 9,-1},
//    { 2, 3,11, 7, 8, 4,-1,-1,-1,-1,-1,-1,-1},
//    { 7,11, 4,11, 2, 4, 4, 2, 0,-1,-1,-1,-1},
//    { 3,11, 2, 4, 7, 8, 9, 0, 1,-1,-1,-1,-1},
//    { 2, 7,11, 2, 1, 7, 1, 4, 7, 1, 9, 4,-1},
//    { 8, 4, 7,11,10, 1,11, 1, 3,-1,-1,-1,-1},
//    {11, 4, 7, 1, 4,11, 1,11,10, 1, 0, 4,-1},
//    { 3, 8, 0, 7,11, 4,11, 9, 4,11,10, 9,-1},
//    { 7,11, 4, 4,11, 9,11,10, 9,-1,-1,-1,-1},
//    { 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 3, 0, 8, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1},
//    { 5, 4, 0, 5, 0, 1,-1,-1,-1,-1,-1,-1,-1},
//    { 4, 8, 5, 8, 3, 5, 5, 3, 1,-1,-1,-1,-1},
//    { 2,10, 1, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1},
//    { 0, 8, 3, 5, 4, 9,10, 1, 2,-1,-1,-1,-1},
//    {10, 5, 2, 5, 4, 2, 2, 4, 0,-1,-1,-1,-1},
//    { 3, 4, 8, 3, 2, 4, 2, 5, 4, 2,10, 5,-1},
//    {11, 2, 3, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 5, 4, 8,11, 2, 8, 2, 0,-1,-1,-1,-1},
//    { 3,11, 2, 1, 5, 4, 1, 4, 0,-1,-1,-1,-1},
//    { 8, 5, 4, 2, 5, 8, 2, 8,11, 2, 1, 5,-1},
//    { 5, 4, 9, 1, 3,11, 1,11,10,-1,-1,-1,-1},
//    { 0, 9, 1, 4, 8, 5, 8,10, 5, 8,11,10,-1},
//    { 3, 4, 0, 3,10, 4, 4,10, 5, 3,11,10,-1},
//    { 4, 8, 5, 5, 8,10, 8,11,10,-1,-1,-1,-1},
//    { 9, 5, 7, 9, 7, 8,-1,-1,-1,-1,-1,-1,-1},
//    { 0, 9, 3, 9, 5, 3, 3, 5, 7,-1,-1,-1,-1},
//    { 8, 0, 7, 0, 1, 7, 7, 1, 5,-1,-1,-1,-1},
//    { 1, 7, 3, 1, 5, 7,-1,-1,-1,-1,-1,-1,-1},
//    { 1, 2,10, 5, 7, 8, 5, 8, 9,-1,-1,-1,-1},
//    { 9, 1, 0,10, 5, 2, 5, 3, 2, 5, 7, 3,-1},
//    { 5, 2,10, 8, 2, 5, 8, 5, 7, 8, 0, 2,-1},
//    {10, 5, 2, 2, 5, 3, 5, 7, 3,-1,-1,-1,-1},
//    {11, 2, 3, 8, 9, 5, 8, 5, 7,-1,-1,-1,-1},
//    { 9, 2, 0, 9, 7, 2, 2, 7,11, 9, 5, 7,-1},
//    { 0, 3, 8, 2, 1,11, 1, 7,11, 1, 5, 7,-1},
//    { 2, 1,11,11, 1, 7, 1, 5, 7,-1,-1,-1,-1},
//    { 3, 9, 1, 3, 8, 9, 7,11,10, 7,10, 5,-1},
//    { 9, 1, 0,10, 7,11,10, 5, 7,-1,-1,-1,-1},
//    { 3, 8, 0, 7,10, 5, 7,11,10,-1,-1,-1,-1},
//    {11, 5, 7,11,10, 5,-1,-1,-1,-1,-1,-1,-1},
//    {10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 3, 0,10, 6, 5,-1,-1,-1,-1,-1,-1,-1},
//    { 0, 1, 9, 5,10, 6,-1,-1,-1,-1,-1,-1,-1},
//    {10, 6, 5, 9, 8, 3, 9, 3, 1,-1,-1,-1,-1},
//    { 1, 2, 6, 1, 6, 5,-1,-1,-1,-1,-1,-1,-1},
//    { 0, 8, 3, 2, 6, 5, 2, 5, 1,-1,-1,-1,-1},
//    { 5, 9, 6, 9, 0, 6, 6, 0, 2,-1,-1,-1,-1},
//    { 9, 6, 5, 3, 6, 9, 3, 9, 8, 3, 2, 6,-1},
//    { 3,11, 2,10, 6, 5,-1,-1,-1,-1,-1,-1,-1},
//    { 6, 5,10, 2, 0, 8, 2, 8,11,-1,-1,-1,-1},
//    { 1, 9, 0, 6, 5,10,11, 2, 3,-1,-1,-1,-1},
//    { 1,10, 2, 5, 9, 6, 9,11, 6, 9, 8,11,-1},
//    {11, 6, 3, 6, 5, 3, 3, 5, 1,-1,-1,-1,-1},
//    { 0, 5, 1, 0,11, 5, 5,11, 6, 0, 8,11,-1},
//    { 0, 5, 9, 0, 3, 5, 3, 6, 5, 3,11, 6,-1},
//    { 5, 9, 6, 6, 9,11, 9, 8,11,-1,-1,-1,-1},
//    {10, 6, 5, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1},
//    { 5,10, 6, 7, 3, 0, 7, 0, 4,-1,-1,-1,-1},
//    { 5,10, 6, 0, 1, 9, 8, 4, 7,-1,-1,-1,-1},
//    { 4, 5, 9, 6, 7,10, 7, 1,10, 7, 3, 1,-1},
//    { 7, 8, 4, 5, 1, 2, 5, 2, 6,-1,-1,-1,-1},
//    { 4, 1, 0, 4, 5, 1, 6, 7, 3, 6, 3, 2,-1},
//    { 9, 4, 5, 8, 0, 7, 0, 6, 7, 0, 2, 6,-1},
//    { 4, 5, 9, 6, 3, 2, 6, 7, 3,-1,-1,-1,-1},
//    { 7, 8, 4, 2, 3,11,10, 6, 5,-1,-1,-1,-1},
//    {11, 6, 7,10, 2, 5, 2, 4, 5, 2, 0, 4,-1},
//    {11, 6, 7, 8, 0, 3, 1,10, 2, 9, 4, 5,-1},
//    { 6, 7,11, 1,10, 2, 9, 4, 5,-1,-1,-1,-1},
//    { 6, 7,11, 4, 5, 8, 5, 3, 8, 5, 1, 3,-1},
//    { 6, 7,11, 4, 1, 0, 4, 5, 1,-1,-1,-1,-1},
//    { 4, 5, 9, 3, 8, 0,11, 6, 7,-1,-1,-1,-1},
//    { 9, 4, 5, 7,11, 6,-1,-1,-1,-1,-1,-1,-1},
//    {10, 6, 4,10, 4, 9,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 3, 0, 9,10, 6, 9, 6, 4,-1,-1,-1,-1},
//    { 1,10, 0,10, 6, 0, 0, 6, 4,-1,-1,-1,-1},
//    { 8, 6, 4, 8, 1, 6, 6, 1,10, 8, 3, 1,-1},
//    { 9, 1, 4, 1, 2, 4, 4, 2, 6,-1,-1,-1,-1},
//    { 1, 0, 9, 3, 2, 8, 2, 4, 8, 2, 6, 4,-1},
//    { 2, 4, 0, 2, 6, 4,-1,-1,-1,-1,-1,-1,-1},
//    { 3, 2, 8, 8, 2, 4, 2, 6, 4,-1,-1,-1,-1},
//    { 2, 3,11, 6, 4, 9, 6, 9,10,-1,-1,-1,-1},
//    { 0,10, 2, 0, 9,10, 4, 8,11, 4,11, 6,-1},
//    {10, 2, 1,11, 6, 3, 6, 0, 3, 6, 4, 0,-1},
//    {10, 2, 1,11, 4, 8,11, 6, 4,-1,-1,-1,-1},
//    { 1, 4, 9,11, 4, 1,11, 1, 3,11, 6, 4,-1},
//    { 0, 9, 1, 4,11, 6, 4, 8,11,-1,-1,-1,-1},
//    {11, 6, 3, 3, 6, 0, 6, 4, 0,-1,-1,-1,-1},
//    { 8, 6, 4, 8,11, 6,-1,-1,-1,-1,-1,-1,-1},
//    { 6, 7,10, 7, 8,10,10, 8, 9,-1,-1,-1,-1},
//    { 9, 3, 0, 6, 3, 9, 6, 9,10, 6, 7, 3,-1},
//    { 6, 1,10, 6, 7, 1, 7, 0, 1, 7, 8, 0,-1},
//    { 6, 7,10,10, 7, 1, 7, 3, 1,-1,-1,-1,-1},
//    { 7, 2, 6, 7, 9, 2, 2, 9, 1, 7, 8, 9,-1},
//    { 1, 0, 9, 3, 6, 7, 3, 2, 6,-1,-1,-1,-1},
//    { 8, 0, 7, 7, 0, 6, 0, 2, 6,-1,-1,-1,-1},
//    { 2, 7, 3, 2, 6, 7,-1,-1,-1,-1,-1,-1,-1},
//    { 7,11, 6, 3, 8, 2, 8,10, 2, 8, 9,10,-1},
//    {11, 6, 7,10, 0, 9,10, 2, 0,-1,-1,-1,-1},
//    { 2, 1,10, 7,11, 6, 8, 0, 3,-1,-1,-1,-1},
//    { 1,10, 2, 6, 7,11,-1,-1,-1,-1,-1,-1,-1},
//    { 7,11, 6, 3, 9, 1, 3, 8, 9,-1,-1,-1,-1},
//    { 9, 1, 0,11, 6, 7,-1,-1,-1,-1,-1,-1,-1},
//    { 0, 3, 8,11, 6, 7,-1,-1,-1,-1,-1,-1,-1},
//    {11, 6, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    {11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 0, 8, 3,11, 7, 6,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 0, 1,11, 7, 6,-1,-1,-1,-1,-1,-1,-1},
//    { 7, 6,11, 3, 1, 9, 3, 9, 8,-1,-1,-1,-1},
//    { 1, 2,10, 6,11, 7,-1,-1,-1,-1,-1,-1,-1},
//    { 2,10, 1, 7, 6,11, 8, 3, 0,-1,-1,-1,-1},
//    {11, 7, 6,10, 9, 0,10, 0, 2,-1,-1,-1,-1},
//    { 7, 6,11, 3, 2, 8, 8, 2,10, 8,10, 9,-1},
//    { 2, 3, 7, 2, 7, 6,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 7, 0, 7, 6, 0, 0, 6, 2,-1,-1,-1,-1},
//    { 1, 9, 0, 3, 7, 6, 3, 6, 2,-1,-1,-1,-1},
//    { 7, 6, 2, 7, 2, 9, 2, 1, 9, 7, 9, 8,-1},
//    { 6,10, 7,10, 1, 7, 7, 1, 3,-1,-1,-1,-1},
//    { 6,10, 1, 6, 1, 7, 7, 1, 0, 7, 0, 8,-1},
//    { 9, 0, 3, 6, 9, 3, 6,10, 9, 6, 3, 7,-1},
//    { 6,10, 7, 7,10, 8,10, 9, 8,-1,-1,-1,-1},
//    { 8, 4, 6, 8, 6,11,-1,-1,-1,-1,-1,-1,-1},
//    {11, 3, 6, 3, 0, 6, 6, 0, 4,-1,-1,-1,-1},
//    { 0, 1, 9, 4, 6,11, 4,11, 8,-1,-1,-1,-1},
//    { 1, 9, 4,11, 1, 4,11, 3, 1,11, 4, 6,-1},
//    {10, 1, 2,11, 8, 4,11, 4, 6,-1,-1,-1,-1},
//    {10, 1, 2,11, 3, 6, 6, 3, 0, 6, 0, 4,-1},
//    { 0, 2,10, 0,10, 9, 4,11, 8, 4, 6,11,-1},
//    { 2,11, 3, 6, 9, 4, 6,10, 9,-1,-1,-1,-1},
//    { 3, 8, 2, 8, 4, 2, 2, 4, 6,-1,-1,-1,-1},
//    { 2, 0, 4, 2, 4, 6,-1,-1,-1,-1,-1,-1,-1},
//    { 1, 9, 0, 3, 8, 2, 2, 8, 4, 2, 4, 6,-1},
//    { 9, 4, 1, 1, 4, 2, 4, 6, 2,-1,-1,-1,-1},
//    { 8, 4, 6, 8, 6, 1, 6,10, 1, 8, 1, 3,-1},
//    { 1, 0,10,10, 0, 6, 0, 4, 6,-1,-1,-1,-1},
//    { 8, 0, 3, 9, 6,10, 9, 4, 6,-1,-1,-1,-1},
//    {10, 4, 6,10, 9, 4,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 5, 4, 7, 6,11,-1,-1,-1,-1,-1,-1,-1},
//    { 4, 9, 5, 3, 0, 8,11, 7, 6,-1,-1,-1,-1},
//    { 6,11, 7, 4, 0, 1, 4, 1, 5,-1,-1,-1,-1},
//    { 6,11, 7, 4, 8, 5, 5, 8, 3, 5, 3, 1,-1},
//    { 6,11, 7, 1, 2,10, 9, 5, 4,-1,-1,-1,-1},
//    {11, 7, 6, 8, 3, 0, 1, 2,10, 9, 5, 4,-1},
//    {11, 7, 6,10, 5, 2, 2, 5, 4, 2, 4, 0,-1},
//    { 7, 4, 8, 2,11, 3,10, 5, 6,-1,-1,-1,-1},
//    { 4, 9, 5, 6, 2, 3, 6, 3, 7,-1,-1,-1,-1},
//    { 9, 5, 4, 8, 7, 0, 0, 7, 6, 0, 6, 2,-1},
//    { 4, 0, 1, 4, 1, 5, 6, 3, 7, 6, 2, 3,-1},
//    { 7, 4, 8, 5, 2, 1, 5, 6, 2,-1,-1,-1,-1},
//    { 4, 9, 5, 6,10, 7, 7,10, 1, 7, 1, 3,-1},
//    { 5, 6,10, 0, 9, 1, 8, 7, 4,-1,-1,-1,-1},
//    { 5, 6,10, 7, 0, 3, 7, 4, 0,-1,-1,-1,-1},
//    {10, 5, 6, 4, 8, 7,-1,-1,-1,-1,-1,-1,-1},
//    { 5, 6, 9, 6,11, 9, 9,11, 8,-1,-1,-1,-1},
//    { 0, 9, 5, 0, 5, 3, 3, 5, 6, 3, 6,11,-1},
//    { 0, 1, 5, 0, 5,11, 5, 6,11, 0,11, 8,-1},
//    {11, 3, 6, 6, 3, 5, 3, 1, 5,-1,-1,-1,-1},
//    { 1, 2,10, 5, 6, 9, 9, 6,11, 9,11, 8,-1},
//    { 1, 0, 9, 6,10, 5,11, 3, 2,-1,-1,-1,-1},
//    { 6,10, 5, 2, 8, 0, 2,11, 8,-1,-1,-1,-1},
//    { 3, 2,11,10, 5, 6,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 5, 6, 3, 9, 6, 3, 8, 9, 3, 6, 2,-1},
//    { 5, 6, 9, 9, 6, 0, 6, 2, 0,-1,-1,-1,-1},
//    { 0, 3, 8, 2, 5, 6, 2, 1, 5,-1,-1,-1,-1},
//    { 1, 6, 2, 1, 5, 6,-1,-1,-1,-1,-1,-1,-1},
//    {10, 5, 6, 9, 3, 8, 9, 1, 3,-1,-1,-1,-1},
//    { 0, 9, 1, 5, 6,10,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 0, 3,10, 5, 6,-1,-1,-1,-1,-1,-1,-1},
//    {10, 5, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    {11, 7, 5,11, 5,10,-1,-1,-1,-1,-1,-1,-1},
//    { 3, 0, 8, 7, 5,10, 7,10,11,-1,-1,-1,-1},
//    { 9, 0, 1,10,11, 7,10, 7, 5,-1,-1,-1,-1},
//    { 3, 1, 9, 3, 9, 8, 7,10,11, 7, 5,10,-1},
//    { 2,11, 1,11, 7, 1, 1, 7, 5,-1,-1,-1,-1},
//    { 0, 8, 3, 2,11, 1, 1,11, 7, 1, 7, 5,-1},
//    { 9, 0, 2, 9, 2, 7, 2,11, 7, 9, 7, 5,-1},
//    {11, 3, 2, 8, 5, 9, 8, 7, 5,-1,-1,-1,-1},
//    {10, 2, 5, 2, 3, 5, 5, 3, 7,-1,-1,-1,-1},
//    { 5,10, 2, 8, 5, 2, 8, 7, 5, 8, 2, 0,-1},
//    { 9, 0, 1,10, 2, 5, 5, 2, 3, 5, 3, 7,-1},
//    { 1,10, 2, 5, 8, 7, 5, 9, 8,-1,-1,-1,-1},
//    { 1, 3, 7, 1, 7, 5,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 7, 0, 0, 7, 1, 7, 5, 1,-1,-1,-1,-1},
//    { 0, 3, 9, 9, 3, 5, 3, 7, 5,-1,-1,-1,-1},
//    { 9, 7, 5, 9, 8, 7,-1,-1,-1,-1,-1,-1,-1},
//    { 4, 5, 8, 5,10, 8, 8,10,11,-1,-1,-1,-1},
//    { 3, 0, 4, 3, 4,10, 4, 5,10, 3,10,11,-1},
//    { 0, 1, 9, 4, 5, 8, 8, 5,10, 8,10,11,-1},
//    { 5, 9, 4, 1,11, 3, 1,10,11,-1,-1,-1,-1},
//    { 8, 4, 5, 2, 8, 5, 2,11, 8, 2, 5, 1,-1},
//    { 3, 2,11, 1, 4, 5, 1, 0, 4,-1,-1,-1,-1},
//    { 9, 4, 5, 8, 2,11, 8, 0, 2,-1,-1,-1,-1},
//    {11, 3, 2, 9, 4, 5,-1,-1,-1,-1,-1,-1,-1},
//    { 3, 8, 4, 3, 4, 2, 2, 4, 5, 2, 5,10,-1},
//    {10, 2, 5, 5, 2, 4, 2, 0, 4,-1,-1,-1,-1},
//    { 0, 3, 8, 5, 9, 4,10, 2, 1,-1,-1,-1,-1},
//    { 2, 1,10, 9, 4, 5,-1,-1,-1,-1,-1,-1,-1},
//    { 4, 5, 8, 8, 5, 3, 5, 1, 3,-1,-1,-1,-1},
//    { 5, 0, 4, 5, 1, 0,-1,-1,-1,-1,-1,-1,-1},
//    { 3, 8, 0, 4, 5, 9,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 4, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 7, 4,11, 4, 9,11,11, 9,10,-1,-1,-1,-1},
//    { 3, 0, 8, 7, 4,11,11, 4, 9,11, 9,10,-1},
//    {11, 7, 4, 1,11, 4, 1,10,11, 1, 4, 0,-1},
//    { 8, 7, 4,11, 1,10,11, 3, 1,-1,-1,-1,-1},
//    { 2,11, 7, 2, 7, 1, 1, 7, 4, 1, 4, 9,-1},
//    { 3, 2,11, 4, 8, 7, 9, 1, 0,-1,-1,-1,-1},
//    { 7, 4,11,11, 4, 2, 4, 0, 2,-1,-1,-1,-1},
//    { 2,11, 3, 7, 4, 8,-1,-1,-1,-1,-1,-1,-1},
//    { 2, 3, 7, 2, 7, 9, 7, 4, 9, 2, 9,10,-1},
//    { 4, 8, 7, 0,10, 2, 0, 9,10,-1,-1,-1,-1},
//    { 2, 1,10, 0, 7, 4, 0, 3, 7,-1,-1,-1,-1},
//    {10, 2, 1, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 1, 4, 4, 1, 7, 1, 3, 7,-1,-1,-1,-1},
//    { 1, 0, 9, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1},
//    { 3, 4, 0, 3, 7, 4,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 9,10, 8,10,11,-1,-1,-1,-1,-1,-1,-1},
//    { 0, 9, 3, 3, 9,11, 9,10,11,-1,-1,-1,-1},
//    { 1,10, 0, 0,10, 8,10,11, 8,-1,-1,-1,-1},
//    {10, 3, 1,10,11, 3,-1,-1,-1,-1,-1,-1,-1},
//    { 2,11, 1, 1,11, 9,11, 8, 9,-1,-1,-1,-1},
//    {11, 3, 2, 0, 9, 1,-1,-1,-1,-1,-1,-1,-1},
//    {11, 0, 2,11, 8, 0,-1,-1,-1,-1,-1,-1,-1},
//    {11, 3, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 3, 8, 2, 2, 8,10, 8, 9,10,-1,-1,-1,-1},
//    { 9, 2, 0, 9,10, 2,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 0, 3, 1,10, 2,-1,-1,-1,-1,-1,-1,-1},
//    {10, 2, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 1, 3, 8, 9, 1,-1,-1,-1,-1,-1,-1,-1},
//    { 9, 1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    { 8, 0, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
//    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}
//};
int Fluid::triTable[256][19] =
{
	{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,0,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//1
	{1,0,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//2
	{8,1,3,8,9,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,1,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//4
	{8,0,3,2,1,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//{0,3,2,8,0,1,3,8,10,2,10,3,10,1,8,1,2,0,-1},
	{2,9,10,2,0,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,2,3,10,2,8,9,10,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{11,3,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//8
	{11,0,2,11,8,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,1,0,3,2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{11,1,2,9,1,11,8,9,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,3,1,10,11,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,0,1,8,0,10,11,8,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,3,0,11,3,9,10,11,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,9,10,8,10,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,4,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//16
	{3,4,0,3,7,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,0,9,4,8,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,4,9,7,4,1,3,7,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,1,10,4,8,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{4,3,7,0,3,4,2,1,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,9,10,0,9,2,4,8,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,2,9,9,2,7,7,2,3,9,7,4,-1,-1,-1,-1,-1,-1,-1},
	{4,8,7,11,3,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{4,11,7,2,11,4,0,2,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,9,1,4,8,7,3,2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,4,11,4,9,11,11,9,2,2,9,1,-1,-1,-1,-1,-1,-1,-1},
	{10,3,1,11,3,10,8,7,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{11,1,10,4,1,11,0,1,4,11,7,4,-1,-1,-1,-1,-1,-1,-1},
	{7,4,8,0,9,11,11,9,10,0,11,3,-1,-1,-1,-1,-1,-1,-1},
	{7,4,11,11,4,9,11,9,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,9,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//32
	{5,9,4,8,0,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,0,4,5,1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,8,4,3,8,5,1,3,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,1,10,5,9,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,3,8,2,1,10,9,4,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,5,10,4,5,2,0,4,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,2,5,2,3,5,5,3,4,4,3,8,-1,-1,-1,-1,-1,-1,-1},
	{5,9,4,3,2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{11,0,2,8,0,11,9,4,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,0,4,1,0,5,3,2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,2,5,5,2,8,8,2,11,8,4,5,-1,-1,-1,-1,-1,-1,-1},
	{3,10,11,1,10,3,5,9,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,4,5,8,0,1,10,8,1,11,8,10,-1,-1,-1,-1,-1,-1,-1},
	{4,5,0,0,5,11,11,5,10,0,11,3,-1,-1,-1,-1,-1,-1,-1},
	{4,5,8,8,5,10,8,10,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,9,8,7,5,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,9,0,5,9,3,7,5,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,0,8,1,0,7,5,1,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,1,3,5,3,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,9,8,5,9,7,1,10,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,10,2,5,9,0,3,5,0,7,5,3,-1,-1,-1,-1,-1,-1,-1},
	{0,8,2,2,8,5,5,8,7,5,10,2,-1,-1,-1,-1,-1,-1,-1},
	{10,2,5,5,2,3,5,3,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,7,5,8,7,9,11,3,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,9,7,7,9,2,2,9,0,7,2,11,-1,-1,-1,-1,-1,-1,-1},
	{3,2,11,1,0,8,7,1,8,5,1,7,-1,-1,-1,-1,-1,-1,-1},
	{2,11,1,1,11,7,1,7,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,9,8,5,8,7,1,10,3,3,10,11,-1,-1,-1,-1,-1,-1,-1},
	{7,5,0,0,5,9,11,7,0,0,1,10,10,11,0,-1,-1,-1,-1},
	{10,11,0,0,11,3,5,10,0,0,8,7,7,5,0,-1,-1,-1,-1},
	{10,11,5,11,7,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,10,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//64
	{8,0,3,10,5,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//{0,3,10,8,0,5,3,8,6,10,6,3,6,5,8,5,10,0,-1},
	{0,9,1,10,5,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,1,3,9,1,8,10,5,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,1,5,6,2,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,1,5,2,1,6,0,3,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,9,5,0,9,6,2,0,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,5,8,8,5,2,2,5,6,2,3,8,-1,-1,-1,-1,-1,-1,-1},
	{3,2,11,6,10,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,11,8,2,11,0,6,10,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,0,9,3,2,11,10,5,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,5,6,9,1,2,11,9,2,8,9,11,-1,-1,-1,-1,-1,-1,-1},
	{3,6,11,5,6,3,1,5,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,0,11,11,0,5,5,0,1,11,5,6,-1,-1,-1,-1,-1,-1,-1},
	{11,3,6,3,0,6,6,0,5,5,0,9,-1,-1,-1,-1,-1,-1,-1},
	{5,6,9,9,6,11,9,11,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,5,6,7,4,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,4,0,7,4,3,5,6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,1,0,10,5,6,4,8,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,10,5,9,1,7,7,1,3,9,7,4,-1,-1,-1,-1,-1,-1,-1},
	{1,6,2,5,6,1,7,4,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,1,5,2,5,6,0,3,4,4,3,7,-1,-1,-1,-1,-1,-1,-1},
	{4,8,7,0,9,5,6,0,5,2,0,6,-1,-1,-1,-1,-1,-1,-1},
	{3,7,9,9,7,4,2,3,9,9,5,6,6,2,9,-1,-1,-1,-1},
	{11,3,2,8,7,4,6,10,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,5,6,7,4,2,2,4,0,7,2,11,-1,-1,-1,-1,-1,-1,-1},
	{1,0,9,7,4,8,3,2,11,10,5,6,-1,-1,-1,-1,-1,-1,-1},
	{2,9,1,11,9,2,4,9,11,11,7,4,10,5,6,-1,-1,-1,-1},
	{4,8,7,11,3,5,5,3,1,11,5,6,-1,-1,-1,-1,-1,-1,-1},
	{1,5,11,11,5,6,0,1,11,11,7,4,4,0,11,-1,-1,-1,-1},
	{5,0,9,6,0,5,3,0,6,6,11,3,4,8,7,-1,-1,-1,-1},
	{5,6,9,9,6,11,7,4,9,11,7,9,-1,-1,-1,-1,-1,-1,-1},
	{4,10,9,4,6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,4,6,9,4,10,8,0,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,10,1,6,10,0,4,6,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,8,1,1,8,6,6,8,4,1,6,10,-1,-1,-1,-1,-1,-1,-1},
	{4,1,9,2,1,4,6,2,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,3,8,2,1,9,4,2,9,6,2,4,-1,-1,-1,-1,-1,-1,-1},
	{2,0,4,2,4,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,8,2,2,8,4,2,4,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{4,10,9,6,10,4,2,11,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,0,2,8,2,11,9,4,10,10,4,6,-1,-1,-1,-1,-1,-1,-1},
	{11,3,2,1,0,6,6,0,4,1,6,10,-1,-1,-1,-1,-1,-1,-1},
	{4,6,1,1,6,10,8,4,1,1,2,11,11,8,1,-1,-1,-1,-1},
	{6,9,4,3,9,6,1,9,3,6,11,3,-1,-1,-1,-1,-1,-1,-1},
	{11,8,1,1,8,0,6,11,1,1,9,4,4,6,1,-1,-1,-1,-1},
	{11,3,6,6,3,0,6,0,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{4,6,8,6,11,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,7,6,8,7,10,9,8,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,0,3,10,0,7,9,0,10,7,6,10,-1,-1,-1,-1,-1,-1,-1},
	{6,10,7,10,1,7,7,1,8,8,1,0,-1,-1,-1,-1,-1,-1,-1},
	{6,10,7,7,10,1,7,1,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,1,6,6,1,8,8,1,9,6,8,7,-1,-1,-1,-1,-1,-1,-1},
	{6,2,9,9,2,1,7,6,9,9,0,3,3,7,9,-1,-1,-1,-1},
	{8,7,0,0,7,6,0,6,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,7,2,7,6,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,2,11,6,10,8,8,10,9,6,8,7,-1,-1,-1,-1,-1,-1,-1},
	{0,2,7,7,2,11,9,0,7,7,6,10,10,9,7,-1,-1,-1,-1},
	{8,1,0,7,1,8,10,1,7,7,6,10,3,2,11,-1,-1,-1,-1},
	{2,11,1,1,11,7,6,10,1,7,6,1,-1,-1,-1,-1,-1,-1,-1},
	{9,8,6,6,8,7,1,9,6,6,11,3,3,1,6,-1,-1,-1,-1},
	{9,0,1,6,11,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,7,0,0,7,6,11,3,0,6,11,0,-1,-1,-1,-1,-1,-1,-1},
	{11,7,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,7,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},//128
	{0,3,8,7,11,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,0,9,7,11,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,8,9,3,8,1,7,11,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,10,2,11,6,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,1,10,0,3,8,11,6,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,2,0,10,2,9,11,6,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{11,6,7,10,2,3,8,10,3,9,10,8,-1,-1,-1,-1,-1,-1,-1},
	{2,7,3,2,6,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,7,8,6,7,0,2,6,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,2,6,3,2,7,1,0,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,1,2,8,1,6,9,1,8,7,8,6,-1,-1,-1,-1,-1,-1,-1},
	{7,10,6,1,10,7,3,1,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,10,6,7,1,10,8,1,7,0,1,8,-1,-1,-1,-1,-1,-1,-1},
	{3,0,7,7,0,10,10,0,9,10,6,7,-1,-1,-1,-1,-1,-1,-1},
	{6,7,10,10,7,8,10,8,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,6,4,8,11,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,3,11,0,3,6,4,0,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,8,11,4,8,6,0,9,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{4,9,6,6,9,3,3,9,1,3,11,6,-1,-1,-1,-1,-1,-1,-1},
	{8,6,4,11,6,8,10,2,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,1,10,0,3,11,6,0,11,4,0,6,-1,-1,-1,-1,-1,-1,-1},
	{11,4,8,6,4,11,2,0,9,10,2,9,-1,-1,-1,-1,-1,-1,-1},
	{9,10,3,3,10,2,4,9,3,3,11,6,6,4,3,-1,-1,-1,-1},
	{2,8,3,4,8,2,6,4,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{4,0,2,6,4,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,1,0,3,2,4,4,2,6,3,4,8,-1,-1,-1,-1,-1,-1,-1},
	{9,1,4,4,1,2,4,2,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,8,3,6,8,1,4,8,6,10,6,1,-1,-1,-1,-1,-1,-1,-1},
	{1,10,0,0,10,6,0,6,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,4,3,3,4,8,10,6,3,3,0,9,9,10,3,-1,-1,-1,-1},
	{9,10,4,10,6,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,4,5,6,7,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,0,3,9,4,5,7,11,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,5,1,4,5,0,6,7,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,11,6,3,8,4,5,3,4,1,3,5,-1,-1,-1,-1,-1,-1,-1},
	{5,9,4,1,10,2,6,7,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{11,6,7,2,1,10,8,0,3,9,4,5,-1,-1,-1,-1,-1,-1,-1},
	{6,7,11,4,5,10,2,4,10,0,4,2,-1,-1,-1,-1,-1,-1,-1},
	{4,3,8,5,3,4,2,3,5,5,10,2,7,11,6,-1,-1,-1,-1},
	{2,7,3,6,7,2,4,5,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,9,4,8,0,6,6,0,2,8,6,7,-1,-1,-1,-1,-1,-1,-1},
	{6,3,2,7,3,6,5,1,0,4,5,0,-1,-1,-1,-1,-1,-1,-1},
	{2,6,8,8,6,7,1,2,8,8,4,5,5,1,8,-1,-1,-1,-1},
	{5,9,4,1,10,6,7,1,6,3,1,7,-1,-1,-1,-1,-1,-1,-1},
	{6,1,10,7,1,6,0,1,7,7,8,0,5,9,4,-1,-1,-1,-1},
	{0,4,10,10,4,5,3,0,10,10,6,7,7,3,10,-1,-1,-1,-1},
	{6,7,10,10,7,8,4,5,10,8,4,10,-1,-1,-1,-1,-1,-1,-1},
	{9,6,5,11,6,9,8,11,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{6,3,11,6,0,3,5,0,6,9,0,5,-1,-1,-1,-1,-1,-1,-1},
	{11,0,8,5,0,11,1,0,5,6,5,11,-1,-1,-1,-1,-1,-1,-1},
	{11,6,3,3,6,5,3,5,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,1,10,5,9,11,11,9,8,5,11,6,-1,-1,-1,-1,-1,-1,-1},
	{11,0,3,6,0,11,9,0,6,6,5,9,2,1,10,-1,-1,-1,-1},
	{8,11,5,5,11,6,0,8,5,5,10,2,2,0,5,-1,-1,-1,-1},
	{11,6,3,3,6,5,10,2,3,5,10,3,-1,-1,-1,-1,-1,-1,-1},
	{8,5,9,2,5,8,6,5,2,8,3,2,-1,-1,-1,-1,-1,-1,-1},
	{5,9,6,6,9,0,6,0,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,1,8,8,1,0,6,5,8,8,3,2,2,6,8,-1,-1,-1,-1},
	{5,1,6,1,2,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,1,6,6,1,10,8,3,6,6,5,9,9,8,6,-1,-1,-1,-1},
	{1,10,0,0,10,6,5,9,0,6,5,0,-1,-1,-1,-1,-1,-1,-1},
	{3,0,8,6,5,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,10,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,11,10,5,7,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,11,10,7,11,5,3,8,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{11,5,7,10,5,11,9,1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,10,5,11,10,7,8,9,1,3,8,1,-1,-1,-1,-1,-1,-1,-1},
	{1,11,2,7,11,1,5,7,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,0,3,2,1,7,7,1,5,2,7,11,-1,-1,-1,-1,-1,-1,-1},
	{7,9,5,2,9,7,0,9,2,11,2,7,-1,-1,-1,-1,-1,-1,-1},
	{5,7,2,2,7,11,9,5,2,2,3,8,8,9,2,-1,-1,-1,-1},
	{5,2,10,3,2,5,7,3,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,8,0,5,8,2,7,8,5,2,10,5,-1,-1,-1,-1,-1,-1,-1},
	{0,9,1,10,5,3,3,5,7,10,3,2,-1,-1,-1,-1,-1,-1,-1},
	{8,9,2,2,9,1,7,8,2,2,10,5,5,7,2,-1,-1,-1,-1},
	{3,1,5,7,3,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,0,7,7,0,1,7,1,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,9,3,3,9,5,3,5,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,9,7,9,5,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,5,4,10,5,8,11,10,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,5,4,11,5,0,10,5,11,3,11,0,-1,-1,-1,-1,-1,-1,-1},
	{1,0,9,4,8,10,10,8,11,4,10,5,-1,-1,-1,-1,-1,-1,-1},
	{11,10,4,4,10,5,3,11,4,4,9,1,1,3,4,-1,-1,-1,-1},
	{5,2,1,8,2,5,11,2,8,5,4,8,-1,-1,-1,-1,-1,-1,-1},
	{4,0,11,11,0,3,5,4,11,11,2,1,1,5,11,-1,-1,-1,-1},
	{2,0,5,5,0,9,11,2,5,5,4,8,8,11,5,-1,-1,-1,-1},
	{4,9,5,11,2,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{5,2,10,5,3,2,4,3,5,8,3,4,-1,-1,-1,-1,-1,-1,-1},
	{10,5,2,2,5,4,2,4,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,3,2,5,3,10,8,3,5,5,4,8,1,0,9,-1,-1,-1,-1},
	{10,5,2,2,5,4,9,1,2,4,9,2,-1,-1,-1,-1,-1,-1,-1},
	{4,8,5,5,8,3,5,3,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{4,0,5,0,1,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{4,8,5,5,8,3,0,9,5,3,0,5,-1,-1,-1,-1,-1,-1,-1},
	{4,9,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{11,4,7,9,4,11,10,9,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,0,3,9,4,7,11,9,7,10,9,11,-1,-1,-1,-1,-1,-1,-1},
	{10,1,11,11,1,4,4,1,0,4,7,11,-1,-1,-1,-1,-1,-1,-1},
	{1,3,4,4,3,8,10,1,4,4,7,11,11,10,4,-1,-1,-1,-1},
	{11,4,7,11,9,4,2,9,11,1,9,2,-1,-1,-1,-1,-1,-1,-1},
	{7,9,4,11,9,7,1,9,11,11,2,1,8,0,3,-1,-1,-1,-1},
	{7,11,4,4,11,2,4,2,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{7,11,4,4,11,2,3,8,4,2,3,4,-1,-1,-1,-1,-1,-1,-1},
	{9,2,10,7,2,9,3,2,7,4,7,9,-1,-1,-1,-1,-1,-1,-1},
	{10,9,7,7,9,4,2,10,7,7,8,0,0,2,7,-1,-1,-1,-1},
	{7,3,10,10,3,2,4,7,10,10,1,0,0,4,10,-1,-1,-1,-1},
	{10,1,2,7,8,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,4,1,1,4,7,1,7,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,4,1,1,4,7,8,0,1,7,8,1,-1,-1,-1,-1,-1,-1,-1},
	{0,4,3,4,7,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{8,4,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,9,8,11,10,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,3,9,9,3,11,9,11,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,0,10,10,0,8,10,8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{1,3,10,3,11,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,1,11,11,1,9,11,9,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{0,3,9,9,3,11,2,1,9,11,2,9,-1,-1,-1,-1,-1,-1,-1},
	{2,0,11,0,8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{2,3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,2,8,8,2,10,8,10,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{10,9,2,9,0,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,2,8,8,2,10,1,0,8,10,1,8,-1,-1,-1,-1,-1,-1,-1},
	{10,1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,1,8,1,9,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{9,0,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{3,0,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}
};
//	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//1
//	{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//2
//	{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//4
//	{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},	//	{0, 3, 2, 8, 0, 1, 3, 8, 10, 2, 10, 3, 10, 1, 8, 1, 2, 0, -1},
//	{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//8
//	{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//16
//	{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1, -1, -1, -1},
//	{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1, -1, -1, -1},
//	{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1, -1, -1, -1},
//	{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//32
//	{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1, -1, -1, -1},
//	{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1, -1, -1, -1},
//	{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
//	{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1, -1, -1, -1},
//	{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
//	{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1, -1, -1, -1},
//	{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1, -1, -1, -1},
//	{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
//	{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1, -1, -1, -1},
//	{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1, -1, -1, -1},
//	{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1, -1, -1, -1},
//	{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//64
//	{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//	{0, 3, 10, 8, 0, 5, 3, 8, 6, 10, 6, 3, 6, 5, 8, 5, 10, 0, -1},
//	{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1, -1, -1, -1},
//	{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
//	{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1, -1, -1, -1},
//	{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1, -1, -1, -1},
//	{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1, -1, -1, -1},
//	{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1, -1, -1, -1},
//	{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
//	{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1, -1, -1, -1},
//	{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1, -1, -1, -1},
//	{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
//	{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1, -1, -1, -1},
//	{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1, -1, -1, -1},
//	{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1, -1, -1, -1},
//	{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1, -1, -1, -1},
//	{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1, -1, -1, -1},
//	{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1, -1, -1, -1},
//	{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
//	{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1, -1, -1, -1},
//	{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1, -1, -1, -1},
//	{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1, -1, -1, -1},
//	{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1, -1, -1, -1},
//	{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1, -1, -1, -1},
//	{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1, -1, -1, -1},
//	{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1, -1, -1, -1},
//	{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1},
//	{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1, -1, -1, -1},
//	{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1, -1, -1, -1},
//	{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1, -1, -1, -1},
//	{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1, -1, -1, -1},
//	{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1, -1, -1, -1},
//	{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1, -1, -1, -1},
//	{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1, -1, -1, -1},
//	{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},//128
//	{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
//	{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1, -1, -1, -1},
//	{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1, -1, -1, -1},
//	{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1, -1, -1, -1},
//	{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1, -1, -1, -1},
//	{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
//	{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1, -1, -1, -1},
//	{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1, -1, -1, -1},
//	{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1, -1, -1, -1},
//	{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1, -1, -1, -1},
//	{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
//	{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
//	{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
//	{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1, -1, -1, -1},
//	{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1, -1, -1, -1},
//	{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1, -1, -1, -1},
//	{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1, -1, -1, -1},
//	{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
//	{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1, -1, -1, -1},
//	{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1, -1, -1, -1},
//	{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1, -1, -1, -1},
//	{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1, -1, -1, -1},
//	{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1, -1, -1, -1},
//	{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1, -1, -1, -1},
//	{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1, -1, -1, -1},
//	{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1, -1, -1, -1},
//	{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1, -1, -1, -1},
//	{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1, -1, -1, -1},
//	{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1, -1, -1, -1},
//	{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1, -1, -1, -1},
//	{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1, -1, -1, -1},
//	{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1, -1, -1, -1},
//	{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1, -1, -1, -1},
//	{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1, -1, -1, -1},
//	{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1, -1, -1, -1},
//	{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1, -1, -1, -1},
//	{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1, -1, -1, -1},
//	{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1, -1, -1, -1},
//	{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1, -1, -1, -1},
//	{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1, -1, -1, -1},
//	{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1, -1, -1, -1},
//	{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1, -1, -1, -1},
//	{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1, -1, -1, -1},
//	{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1, -1, -1, -1},
//	{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1, -1, -1, -1},
//	{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1, -1, -1, -1},
//	{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1, -1, -1, -1},
//	{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1, -1, -1, -1},
//	{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
//	{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1, -1, -1, -1},
//	{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1, -1, -1, -1},
//	{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1, -1, -1, -1},
//	{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1, -1, -1, -1},
//	{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1, -1, -1, -1},
//	{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1, -1, -1, -1},
//	{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1, -1, -1, -1},
//	{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1, -1, -1, -1},
//	{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1, -1, -1, -1},
//	{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1, -1, -1, -1},
//	{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
//	{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
//};
#endif /*FLUID_H*/
