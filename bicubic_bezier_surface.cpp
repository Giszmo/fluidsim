#ifndef NULL
#define NULL 0
#endif
#include <math.h>
#include "vektor.h"

#include "bicubic_bezier_surface.h"

bicubic_bezier_surface::bicubic_bezier_surface ( unsigned char * rd )
{
	_rawdata=rd;
	int i;
	bbp_ = NULL;
	//  f_ = NULL;
	koi_ = 2;
	oldfX_ = 10;
	mode_=1;
	dim_ = 100;
	bbp_ = new bicubic_bezier_patch*[10];
	for ( i=0;i<10;i++ )
		bbp_[i] = new bicubic_bezier_patch[10];
}

bicubic_bezier_surface::~bicubic_bezier_surface()
{
//   int i;
//   if(bbp_)
//     {
	// Speicher wieder freigeben
	//		for(i=0;i<fY_;i++)
	//			delete[] bbp_[i];
	//		delete[] bbp_;

//     }

}

float
bicubic_bezier_surface::operator() ( float u,float v ) //nur werte zwischen 0 und 1
{
	int nrX,nrY;
	nrX = ( int ) floor ( u*fX_ );
	nrY = ( int ) floor ( v*fY_ );
	return    bbp_[nrX][nrY] ( u*fX_-floor ( u*fX_ ),v*fY_-floor ( v*fY_ ) );
}

float
bicubic_bezier_surface::dx ( float u,float v ) //nur werte zwischen 0 und 1
{
	int nrX,nrY;
	nrX = ( int ) floor ( u*fX_ );
	nrY = ( int ) floor ( v*fY_ );
	return bbp_[nrX][nrY].dx ( u*fX_-floor ( u*fX_ ),v*fY_-floor ( v*fY_ ) ); //*fX_;
}

float
bicubic_bezier_surface::dy ( float u,float v ) //nur werte zwischen 0 und 1
{
	int nrX,nrY;
	nrX = ( int ) floor ( u*fX_ );
	nrY = ( int ) floor ( v*fY_ );
	return bbp_[nrX][nrY].dy ( u*fX_-floor ( u*fX_ ),v*fY_-floor ( v*fY_ ) ); //*fY_;
}

void bicubic_bezier_surface::callf ( int nX,int nY )
{
	int i,j;
	float ps;
	fX_ = nX;fY_ = nY;

	if ( bbp_!=NULL )
	{
		for ( i=0;i<oldfX_;i++ )
			if ( bbp_!=NULL )
			{
				delete[] bbp_[i];
			}
		delete[] bbp_;
	}
	bbp_ = new bicubic_bezier_patch*[nY];

	for ( i=0;i<nY;i++ )
		bbp_[i] = new bicubic_bezier_patch[nX];
	oldfX_=nX_;

	// setzen der Eck-Interpolationspunkte
	ps=_rawdata[   0];	bbp_[0   ][0   ].setp ( 0,0,ps );
	ps=_rawdata[  99];	bbp_[0   ][nY-1].setp ( 0,3,ps );
	ps=_rawdata[9900];	bbp_[nX-1][0   ].setp ( 3,0,ps );
	ps=_rawdata[9999];	bbp_[nX-1][nY-1].setp ( 3,3,ps );


	// setzen der linken und rechten Interpolationspunkte
	for ( i=1;i<nY;i++ )
	{
		ps=_rawdata[  0*100+ ( int ) floor ( 99*float ( i ) /float ( nY ) ) ];
		bbp_[0][i].setp ( 0,0,ps );
		bbp_[0][i-1].setp ( 3,0,ps );

		ps=_rawdata[99*100+ ( int ) floor ( 99*float ( i ) /float ( nY ) ) ];
		bbp_[nX-1][i].setp ( 0,3,ps );
		bbp_[nX-1][i-1].setp ( 3,3,ps );
	}

	// setzen der oberen und unteren Interpolationspunkte
	for ( j=1;j<nX;j++ )
	{
		ps=_rawdata[ ( int ) floor ( 99*float ( j ) /float ( nX ) ) *100+0];
		bbp_[j][0].setp ( 0,0,ps );
		bbp_[j-1][0].setp ( 0,3,ps );

		ps=_rawdata[ ( int ) floor ( 99*float ( j ) /float ( nX ) ) *100+99];
		bbp_[j][nY-1].setp ( 3,0,ps );
		bbp_[j-1][nY-1].setp ( 3,3,ps );
	}

	// setzen der inneren Interpolationspunkte
	for ( i=1;i<nY;i++ )
	{
		for ( j=1;j<nX;j++ )
		{
			ps=_rawdata[ ( int ) floor ( 99*float ( j ) /float ( nX ) ) *100+ ( int ) floor ( 99*float ( i ) /float ( nY ) ) ];
			bbp_[j][i].setp ( 0,0,ps );
			bbp_[j][i-1].setp ( 3,0,ps );
			bbp_[j-1][i].setp ( 0,3,ps );
			bbp_[j-1][i-1].setp ( 3,3,ps );
		}
	}
}


void bicubic_bezier_surface::interpolate() //Cardinal
{
	int i,j;//,k;
//  double h;
	float ps;

	// setzen der inneren Randbezierpunkte
	// 					  1  2
	//					4      7	innen
	//					8     11
	//					 13 14
	for ( i=1;i<fY_;i++ )
	{
		for ( j=1;j<fX_;j++ )
		{
			ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j][i].getps ( 0,3 )-bbp_[j-1][i].getps ( 0,0 ) ) /6;
			bbp_[j][i].setp ( 0,1,ps );		bbp_[j][i-1].setp ( 3,1,ps );
			ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j][i].getps ( 3,0 )-bbp_[j][i-1].getps ( 0,0 ) ) /6;
			bbp_[j][i].setp ( 1,0,ps );		bbp_[j-1][i].setp ( 1,3,ps );

			ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j-1][i].getps ( 0,0 )-bbp_[j][i].getps ( 0,3 ) ) /6;
			bbp_[j-1][i].setp ( 0,2,ps );		bbp_[j-1][i-1].setp ( 3,2,ps );
			ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j][i-1].getps ( 0,0 )-bbp_[j][i].getps ( 3,0 ) ) /6;
			bbp_[j][i-1].setp ( 2,0,ps );		bbp_[j-1][i-1].setp ( 2,3,ps );
		}
	}

	// setzen der linken und rechten Randbezierpunkte
	// 					  1  2
	//					4      7	li/re
	//					8     11
	//					 13 14
	for ( i=1;i<fY_;i++ )
	{
		j=0;
		ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j][i].getps ( 3,0 )-bbp_[j][i-1].getps ( 0,0 ) ) /6;
		bbp_[j][i].setp ( 1,0,ps );
		ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j][i-1].getps ( 0,0 )-bbp_[j][i].getps ( 3,0 ) ) /6;
		bbp_[j][i-1].setp ( 2,0,ps );
		ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j][i].getps ( 0,3 )-bbp_[j][i].getps ( 0,0 ) ) /3;
		bbp_[j][i].setp ( 0,1,ps );	bbp_[j][i-1].setp ( 3,1,ps );

		j=fX_-1;
		ps=bbp_[j][i].getps ( 0,3 ) + ( bbp_[j][i].getps ( 3,3 )-bbp_[j][i-1].getps ( 0,3 ) ) /6;
		bbp_[j][i].setp ( 1,3,ps );
		ps=bbp_[j][i].getps ( 0,3 ) + ( bbp_[j][i-1].getps ( 0,3 )-bbp_[j][i].getps ( 3,3 ) ) /6;
		bbp_[j][i-1].setp ( 2,3,ps );
		ps=bbp_[j][i].getps ( 0,3 ) + ( bbp_[j][i].getps ( 0,0 )-bbp_[j][i].getps ( 0,3 ) ) /3;
		bbp_[j][i].setp ( 0,2,ps );	bbp_[j][i-1].setp ( 3,2,ps );
	}

	// setzen der oberen und unteren Randbezierpunkte
	for ( j=1;j<fX_;j++ )
	{
		i=0;
		ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j-1][i].getps ( 0,0 )-bbp_[j][i].getps ( 0,3 ) ) /6;
		bbp_[j-1][i].setp ( 0,2,ps );
		ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j][i].getps ( 0,3 )-bbp_[j-1][i].getps ( 0,0 ) ) /6;
		bbp_[j][i].setp ( 0,1,ps );
		ps=bbp_[j][i].getps ( 0,0 ) + ( bbp_[j][i].getps ( 3,0 )-bbp_[j][i].getps ( 0,0 ) ) /3;
		bbp_[j][i].setp ( 1,0,ps );	bbp_[j-1][i].setp ( 1,3,ps );

		i=fY_-1;
		ps=bbp_[j][i].getps ( 3,0 ) + ( bbp_[j-1][i].getps ( 3,0 )-bbp_[j][i].getps ( 3,3 ) ) /6;
		bbp_[j-1][i].setp ( 3,2,ps );
		ps=bbp_[j][i].getps ( 3,0 ) + ( bbp_[j][i].getps ( 3,3 )-bbp_[j-1][i].getps ( 3,0 ) ) /6;
		bbp_[j][i].setp ( 3,1,ps );
		ps=bbp_[j][i].getps ( 3,0 ) + ( bbp_[j][i].getps ( 0,0 )-bbp_[j][i].getps ( 3,0 ) ) /3;
		bbp_[j][i].setp ( 2,0,ps );	bbp_[j-1][i].setp ( 2,3,ps );
	}

	ps=bbp_[0][0].getps ( 0,0 );			bbp_[0][0].setp ( 0,1,ps );
	ps=bbp_[fX_-1][0].getps ( 0,3 );		bbp_[fX_-1][0].setp ( 0,2,ps );

	ps=bbp_[0][0].getps ( 0,0 );			bbp_[0][0].setp ( 1,0,ps );
	ps=bbp_[fX_-1][0].getps ( 0,3 );		bbp_[fX_-1][0].setp ( 1,3,ps );

	ps=bbp_[0][fY_-1].getps ( 3,0 );		bbp_[0][fY_-1].setp ( 2,0,ps );
	ps=bbp_[fX_-1][fY_-1].getps ( 3,3 );	bbp_[fX_-1][fY_-1].setp ( 2,3,ps );

	ps=bbp_[0][fY_-1].getps ( 3,0 );		bbp_[0][fY_-1].setp ( 3,1,ps );
	ps=bbp_[fX_-1][fY_-1].getps ( 3,3 );	bbp_[fX_-1][fY_-1].setp ( 3,2,ps );

	// setzen der inneren Bezierpunkte
	for ( i=0;i<fY_;i++ )
	{
		for ( j=0;j<fX_;j++ )
		{
			ps=bbp_[j][i].getps ( 0,1 ) +bbp_[j][i].getps ( 1,0 )-bbp_[j][i].getps ( 0,0 );
			bbp_[j][i].setp ( 1,1,ps );
			ps=bbp_[j][i].getps ( 0,2 ) +bbp_[j][i].getps ( 1,3 )-bbp_[j][i].getps ( 0,3 );
			bbp_[j][i].setp ( 1,2,ps );
			ps=bbp_[j][i].getps ( 2,0 ) +bbp_[j][i].getps ( 3,1 )-bbp_[j][i].getps ( 3,0 );
			bbp_[j][i].setp ( 2,1,ps );
			ps=bbp_[j][i].getps ( 3,2 ) +bbp_[j][i].getps ( 2,3 )-bbp_[j][i].getps ( 3,3 );
			bbp_[j][i].setp ( 2,2,ps );
		}
	}
}

float
bicubic_bezier_surface::gets ( int px, int py, int pi, int pj )
{
	return ( bbp_[px][py].getps ( pi, pj ) );
}

float
bicubic_bezier_surface::geth ( int x, int y )
{
	if ( x < fX_ && y < fY_ )
		return ( bbp_[x][y].getps ( 0,0 ) );

	if ( x == fX_ || y == fY_ )
	{
		if ( x > 0 )
			return ( bbp_[x-1][y].getps ( 0,3 ) );
		if ( y > 0 )
			return ( bbp_[x][y-1].getps ( 3,0 ) );
		if ( x == fX_ && y == fY_ )
			return ( bbp_[x-1][y-1].getps ( 3,3 ) );/////////?????? falsche reihenfolge
	}
	return ( 99 );
}

void
bicubic_bezier_surface::changepoint ( int x, int y, float h )
{
	float ps;
	ps=geth ( x,y ) +h;
	if ( x > 0 && y > 0 )
	{
		bbp_[x-1][y-1].setp ( 3,3,ps );
	}
	if ( x < fX_ && y < fY_ )
	{
		bbp_[x  ][y  ].setp ( 0,0,ps );
	}
	if ( x > 0 && y < fY_ )
	{
		bbp_[x-1][y  ].setp ( 0,3,ps );
	}
	if ( x < fX_ && y > 0 )
	{
		bbp_[x  ][y-1].setp ( 3,0,ps );
	}
}

void
bicubic_bezier_surface::solve ( int nX,int nY,float **A )
{
	int i,j;//,k;
	Vektor p;
	for ( i=0;i<nY*fY_;++i )
	{
		for ( j=0;j<nX*fX_;++j )
		{
			A[i][j]=bbp_[int ( floor ( float ( j ) /float ( nX ) ) ) ][int ( floor ( float ( i ) /float ( nY ) ) ) ] ( float ( j%nX ) /float ( nX ),float ( i%nY ) /float ( nY ) );
		}
		A[i][j]=bbp_[fX_-1][int ( floor ( float ( i ) /float ( nY ) ) ) ] ( 1,float ( i%nY ) /float ( nY ) );
	}
	i--;
	j--;
	for ( j=0;j<nX*fX_;++j )
	{
		A[i][nY*fY_-1]=bbp_[int ( floor ( float ( j ) /float ( nX ) ) ) ][fY_-1] ( float ( j%nX ) /float ( nX ),1 );
	}
	A[nX*fX_-1][nY*fY_-1]=bbp_[fX_-1][fY_-1] ( 1,1 );
}

void
bicubic_bezier_surface::outputraw()
{
	for ( int i=0; i<=10; i++ )
	{
		for ( int j=0; j<=10; j++ )
		{
			printf ( "%xh ",_rawdata[i*100+j] );
		}
		printf ( "\n" );
	}
}

void
bicubic_bezier_surface::outputpatches()
{
	for ( int i=0; i<=0; i++ )
	{
		for ( int j=0; j<=4; j++ )
		{
			for ( int k=0; k<4; k++ )
			{
				for ( int l=0; l<4; l++ )
				{
					printf ( "%2.2f ",bbp_[i][j].getps ( k,l ) );
				}
				printf ( "\n" );
			}
			printf ( "\n" );
		}
		printf ( "\n" );
	}
}
Vektor
bicubic_bezier_surface::normal ( float u, float v )
{
	int nrX,nrY;
	nrX = ( int ) floorf ( u*fX_ );
	nrY = ( int ) floorf ( v*fY_ );
	return bbp_[nrX][nrY].normal ( u*fX_-floor ( u*fX_ ),v*fY_-floor ( v*fY_ ) ); //*fX_;
}
