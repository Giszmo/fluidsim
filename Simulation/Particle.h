#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vektor.h"

//kinds of (control)particles
#define MOVING    0
#define BOUNDARY  1
#define TELEPORT  2
#define SHIFT     3
#define SETSPEED  4
//#define DELETE    5

class Particle {
private:
    Vektor _x; //  x  y  z
    Vektor _v; // vx vy vz
    Vektor _a; // ax ay az acceleration for particles, outward normal for boundary
	unsigned char _kind; // moving/boundary/teleport/delete/setspeed
    //float _m; // masse
    float _rho; // dichte ...
    unsigned long _sortlistindex[4];
    unsigned char _relevantsortlist;
public:
    Particle()
			:_x(Vektor(0,0,0)),_v(Vektor(0,0,0)),_a(Vektor(0,0,0)),_kind(100)/*,_m(1)*/,_rho(0), _relevantsortlist(0) {
	}
    Particle(const Vektor x, const Vektor v, const Vektor a, const unsigned char kind/*, const float m*/, const float rho)
            :_x(x),_v(v),_a(a),_kind(kind)/*,_m(m)*/,_rho(rho), _relevantsortlist(0) {
	};
    ~Particle() {}
    void set(const float x, const float y, const float z, const float vx, const float vy, const float vz, const float r, /*const float m, */const float rho, const unsigned char kind) {
        _x.set(x,y,z);
        _v.set(vx,vy,vz);
        //_r = r;
        //_m = m;
        _rho = rho;
		_kind = kind;
    };
    void set(const float x, const float y, const float z, const float ax, const float ay, const float az, const float vx, const float vy, const float vz, const unsigned char kind) {
        _x.set(x,y,z);
        _v.set(vx,vy,vz);
        _a.set(ax,ay,az);
		_kind = kind;
    };
	void setsortlistindex(const unsigned int i, const unsigned long j) {
		_sortlistindex[i]=j;
	}
	bool outsidebox(	const float minx, const float miny, const float minz,
						const float maxx, const float maxy, const float maxz) {
		return ((_x.x() < minx) || (_x.x() > maxx) || 
				(_x.y() < miny) || (_x.y() > maxy) || 
				(_x.z() < minz) || (_x.z() > maxz));
	}
	unsigned char & relevantsortlist() {
		return _relevantsortlist;
	}
	unsigned long get_sortlistindex(const unsigned char sl) {
		return _sortlistindex[sl];
	}
	void set_sortlistindex(const unsigned char sl, const unsigned long index) {
		_sortlistindex[sl]=index;
	}
	void set_relevantsortlist(const unsigned char rsl) {
		_relevantsortlist=rsl;
	}
    void setx(const Vektor x) {
        _x = x;
    };
    void setv(const Vektor v) {
        _v = v;
    };
    void setkind(const unsigned char kind) {
        _kind = kind;
    };
  //  void setm(const float m) {
		//_m = m;
  //  };
    Vektor & x() {
        return _x;
    };
    Vektor & v() {
        return _v;
    };
    Vektor & a() {
        return _a;
    };
    //float & r() {
    //    return _r;
    //};
    //float & m() {
    //    return _m;
    //};
    float & rho() {
        return _rho;
    };
	unsigned char & kind() {
		return _kind;
	};
	void resetrho() {_rho = 0;};
	void setrho(float & rho) {_rho =rho;};
	void addrho(float & rho) {_rho+=rho;};
    void acc(const Vektor a)  {
        _a = a;
    };
    void push(Vektor f) {
        _a = _a + f;// / _rho;
    };
    void move(const float t) {
		_a += _v * t;//.norm(_v.absabs()*(-.01f)/*/_m*/);
//		if(_a.absabs()>400)
//			_a=_a.norm(20);
        _v += _a * t;
		if(_v.absabs()>400)
			_v=_v.norm(20);
        _x += _v * t;
//		_x.sety(0);
    };

};
#endif /*PARTICLE_H*/
