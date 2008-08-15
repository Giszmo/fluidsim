#ifndef _BICUB_BEZIER_SURFACE_H_
#define _BICUB_BEZIER_SURFACE_H_

#include "bicubic_bezier_patch.h"

class bicubic_bezier_surface/*: public bezier_surface*/
{
	private:
		int mode_; /* edit bzw show */
		int fX_, fY_, oldfX_;				/*  Patches in x / y Richtung / Patches in x Richtung vor letzter Aenderung   (n ist f_calls-1) */
		int nX_, nY_;						/*  Auswertungen in x / y Richtung pro Patch */
		int koi_;							/*  kind of interpolation; */
		int dim_;
		unsigned char * _rawdata;
		bicubic_bezier_patch** bbp_;

	public:
		bicubic_bezier_surface ( unsigned char * rd );
		bicubic_bezier_surface() {};
		~bicubic_bezier_surface();

		float operator() ( float u,float v );
		float dx ( float u,float v );
		float dy ( float u,float v );
		float gets ( int px, int py, int pi, int pj ); /*  Patch_X,Patch_Y,Patch_Index */
		float geth ( int px, int py );
		int fX ( void ) {return fX_;}
		int fY ( void ) {return fY_;}
		void callf ( int nX,int nY );         /* setzt die f_ausgewerteten Punkte 0,3,12 und 15 der (nX_ = nX)*(nY_ = nY) Patches */
		void changepoint ( int x, int y, float h );
		void interpolate();
		void solve ( int nX,int nY,float **A );
		void outputraw();
		void outputpatches();
		void interpolateCardinal();
		Vektor normal ( float u, float v );
};

#endif/* _BICUB_BEZIER_SURFACE_H_ */
