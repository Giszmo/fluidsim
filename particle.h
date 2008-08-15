#ifndef PARTICLE_H
#define PARTICLE_H

#include "vektor.h"

//kinds of (control)particles
#define MOVING    0
#define BOUNDARY  1
#define TELEPORT  2
#define SHIFT     3
#define SETSPEED  4
//#define DELETE    5

enum PARTICLE_TYPE
{
	undefined=0,
	moving   =1,///freely moving particles
	boundary =2,///don't move and invert the vertical speed of colliding particles
	teleport =3,///don't move and transport colliding particles to v
	p_remove =4,///don't move and soak up colliding particles
	setspeed =5,///don't move and set the speed of colliding particles to v
	shift    =6,///don't move and add v to colliding particles position
};

class Particle
{
	private:
		Vektor _x; //  x  y  z position
		Vektor _v; // vx vy vz velocity
		Vektor _a; // ax ay az acceleration for particles, outward normal for boundary
		PARTICLE_TYPE _particleType; // moving/boundary/teleport/delete/setspeed
		float _rho; // density
		unsigned long _sortlistindex[4];
		unsigned char _relevantsortlist;
	public:
		Particle()
				:_x ( Vektor ( 0,0,0 ) ),_v ( Vektor ( 0,0,0 ) ),_a ( Vektor ( 0,0,0 ) ),_particleType ( undefined ),_rho ( 0 ), _relevantsortlist ( 0 )
		{
		}
		Particle ( const Vektor x, const Vektor v, const Vektor a, const PARTICLE_TYPE kind, const float rho )
				:_x ( x ),_v ( v ),_a ( a ),_particleType ( kind ),_rho ( rho ), _relevantsortlist ( 0 )
		{
		};
		~Particle() {}
		void set ( const float x, const float y, const float z, const float vx, const float vy, const float vz, const float r, /*const float m, */const float rho, const PARTICLE_TYPE kind )
		{
			_x.set ( x,y,z );
			_v.set ( vx,vy,vz );
			_rho = rho;
			_particleType = kind;
		};
		void set ( const float x, const float y, const float z, const float ax, const float ay, const float az, const float vx, const float vy, const float vz, const PARTICLE_TYPE kind )
		{
			_x.set ( x,y,z );
			_v.set ( vx,vy,vz );
			_a.set ( ax,ay,az );
			_particleType = kind;
		};
		void setsortlistindex ( const unsigned int i, const unsigned long j )
		{
			_sortlistindex[i]=j;
		}
		bool outsidebox (	const float minx, const float miny, const float minz,
		                  const float maxx, const float maxy, const float maxz )
		{
			return ( ( _x.x() < minx ) || ( _x.x() > maxx ) ||
			         ( _x.y() < miny ) || ( _x.y() > maxy ) ||
			         ( _x.z() < minz ) || ( _x.z() > maxz ) );
		}
		unsigned char & relevantsortlist()
		{
			return _relevantsortlist;
		}
		unsigned long get_sortlistindex ( const unsigned char sl )
		{
			return _sortlistindex[sl];
		}
		void set_sortlistindex ( const unsigned char sl, const unsigned long index )
		{
			_sortlistindex[sl]=index;
		}
		void set_relevantsortlist ( const unsigned char rsl )
		{
			_relevantsortlist=rsl;
		}
		void setx ( const Vektor x )
		{
			_x = x;
		};
		void setv ( const Vektor v )
		{
			_v = v;
		};
		void setkind ( const PARTICLE_TYPE kind )
		{
			_particleType = kind;
		};
		Vektor & x()
		{
			return _x;
		};
		Vektor & v()
		{
			return _v;
		};
		Vektor & a()
		{
			return _a;
		};
		float & rho()
		{
			return _rho;
		};
		PARTICLE_TYPE & kind()
		{
			return _particleType;
		};
		void resetrho() {_rho = 0;};
		void setrho ( float & rho ) {_rho =rho;};
		void addrho ( float & rho ) {_rho+=rho;};
		void acc ( const Vektor a )
		{
			_a = a;
		};
		void push ( Vektor f )
		{
			_a = _a + f;
		};
		void move ( const float t )
		{
			_a += _v * t;//.norm(_v.absabs()*(-.01f)/*/_m*/);
			_v += _a * t;
			_x += _v * t;
		};
};
#endif /*PARTICLE_H*/
