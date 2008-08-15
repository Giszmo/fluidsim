#ifndef VEKTOR_H
#define VEKTOR_H

#include <math.h>
#include <float.h>

inline void printspaces ( int n );
class Vektor
{
		float _x, _y, _z;
	public:
		float & x() {return _x;};
		float & y() {return _y;};
		float & z() {return _z;};
		Vektor() {};
		Vektor ( float x, float y, float z ) :_x ( x ),_y ( y ),_z ( z ) {};
		Vektor ( float * v ) :_x ( v[0] ),_y ( v[1] ),_z ( v[2] ) {};
		~Vektor() {};
		Vektor operator+ ( const Vektor b )
		{
			return Vektor (
			           _x + b._x,
			           _y + b._y,
			           _z + b._z );
		}
		void operator+= ( const Vektor b )
		{
			_x += b._x;
			_y += b._y;
			_z += b._z;
		}
		Vektor operator- ( const Vektor b )
		{
			return Vektor (
			           _x - b._x,
			           _y - b._y,
			           _z - b._z );
		}
		Vektor operator% ( const Vektor B )  // Kreuzprodukt
		{
			return Vektor (
			           _y * B._z - _z * B._y,
			           _z * B._x - _x * B._z,
			           _x * B._y - _y * B._x );
		}
		Vektor operator* ( const float b )  // Vektorprodukt
		{
			return Vektor ( _x * b,_y * b,_z * b );
		}
		Vektor operator/ ( const float b )  // Vektorprodukt
		{
			return Vektor ( _x / b,_y / b,_z / b );
		}
		float operator* ( const Vektor B )  // Skalarprodukt
		{
			return _x*B._x+_y*B._y+_z*B._z;
		}
		bool operator== ( const Vektor i ) const
		{
			return ( _x == i._x ) && ( _y == i._y ) && ( _z == i._z );
		}
		bool operator!= ( const Vektor i ) const
		{
			return ( _x != i._x ) || ( _y != i._y ) || ( _z != i._z );
		}
		void write2array ( float * v )
		{
			v[0]=_x;
			v[1]=_y;
			v[2]=_z;
		};
		void set ( const float x, const float y, const float z )
		{
			_x = x;
			_y = y;
			_z = z;
		}
		void setx ( const float x )
		{
			_x = x;
		}
		void sety ( const float y )
		{
			_y = y;
		}
		void setz ( const float z )
		{
			_z = z;
		}
		float abs()
		{
			return sqrtf ( _x*_x + _y*_y + _z*_z );
		}
		float absabs()
		{
			return _x*_x + _y*_y + _z*_z;
		}
		Vektor normed()
		{
			return Vektor ( _x,_y,_z ) * ( 1/sqrtf ( _x*_x + _y*_y + _z*_z ) );
		}
		Vektor norm ( const float n )
		{
			return Vektor ( _x,_y,_z ) * ( n/sqrtf ( _x*_x + _y*_y + _z*_z ) );
		}
};
#endif /*VEKTOR_H*/
