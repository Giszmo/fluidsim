#include "vektor.h"
#include <iostream>
#include <cstring>

class bicubic_bezier_patch  //virtuelle kantenl�nge 1
{
	private:
		float m_ps[4][4];
		float m_sizeX,m_sizeY;
	public:
		bicubic_bezier_patch();
		~bicubic_bezier_patch();
		void set ( float** );
		void setp ( int, int, float );
		float** get();
		float getp ( int, int );
		float getps ( int, int );
		float operator() ( float u,float v );
		float dx ( float u,float v );
		float dy ( float u,float v );
		Vektor normal ( float u,float v );
};
