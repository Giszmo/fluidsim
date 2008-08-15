#include "bicubic_bezier_patch.h"

#ifndef NULL
#define NULL 0
#endif


bicubic_bezier_patch::bicubic_bezier_patch()
{
}

bicubic_bezier_patch::~bicubic_bezier_patch()
{
}

void
bicubic_bezier_patch::set ( float** ps )
{
	int i,j;
	for ( i=0;i<=3;i++ )
	{
		for ( j=0;j<=3;j++ )
		{
			m_ps[i][j]=ps[i][j];
		}
	}
}

void
bicubic_bezier_patch::setp ( int i, int j, float ps )
{
	m_ps[i][j]=ps;
}

float** bicubic_bezier_patch::get()
{
	int i,j;
	float** ps;
	ps = new float*[4];
	for ( i=0;i<=3;i++ )
		ps[i] = new float[4];

	for ( i=0;i<=3;i++ )
	{
		for ( j=0;j<=3;j++ )
		{
			ps[i][j]=m_ps[i][j];
		}
	}
	return ps;
}

float bicubic_bezier_patch::getp ( int i, int j )
{
	float ps;
	ps=m_ps[i][j];
	return ps;
}

float bicubic_bezier_patch::getps ( int i, int j )
{
	return m_ps[i][j];
}

float
bicubic_bezier_patch::operator() ( float u,float v )
{
	float h[4][4];
	int i,j,k=0;
	memcpy ( h,m_ps,16*sizeof ( float ) );
	//for(i=0;i<4;i++)		for(j=0;j<4;j++)			h[i][j]=m_ps[i][j];

	for ( k=3;k>=0;k-- )
	{
		for ( i=0;i<k;i++ )
		{
			for ( j=0;j<k;j++ )
			{
				h[j][i]= ( 1-u ) * ( ( 1-v ) *h[j][i]+v*h[j+1][i] ) +u* ( ( 1-v ) *h[j][i+1]+v*h[j+1][i+1] );
			}
		}
	}
	return h[0][0];
}

float
bicubic_bezier_patch::dx ( float u,float v )
{
	float h[4][4];
	int i,j,k=0;
	for ( i=0;i<4;i++ )
		for ( j=0;j<4;j++ )
			h[j][i]=m_ps[i][j];

	for ( k=3;k>=1;k-- )
	{
		for ( i=0;i<k;i++ )
		{
			for ( j=0;j<k;j++ )
			{
				h[i][j]= ( 1-u ) * ( ( 1-v ) *h[i  ][j  ]+
				                     v    *h[i  ][j+1] ) +
				         u    * ( ( 1-v ) *h[i+1][j  ]+
				                  v    *h[i+1][j+1] );
			}
		}
	}
	/*	return    (((1-v)*h[1][1]+
					v    *h[1][0])-
				   ((1-v)*h[0][1]+
				    v    *h[0][0]))*.03;*/
	return ( ( 1-v ) * ( h[1][0]-h[0][0] ) + v * ( h[1][1]-h[0][1] ) ) / 3;
}

float
bicubic_bezier_patch::dy ( float u,float v )
{
	float h[4][4];
	int i,j,k=0;
	for ( i=0;i<4;i++ )
		for ( j=0;j<4;j++ )
			h[j][i]=m_ps[i][j];

	for ( k=3;k>=1;k-- )
	{
		for ( i=0;i<k;i++ )
		{
			for ( j=0;j<k;j++ )
			{
				h[i][j]= ( 1-u ) * ( ( 1-v ) *h[i][j]+v*h[i][j+1] ) +u* ( ( 1-v ) *h[i+1][j]+v*h[i+1][j+1] );
			}
		}
	}
//	return (((1-u)*h[1][1]+u*h[0][1])-((1-u)*h[1][0]+u*h[0][0]))*.03;
	return ( ( 1-u ) * ( h[0][1]-h[0][0] ) + u * ( h[1][1]-h[1][0] ) ) / 3;
}
Vektor
bicubic_bezier_patch::normal ( float u,float v )
{
	float h[4][4];
	memcpy ( h,m_ps,16*sizeof ( float ) );

	h[0][0]= ( 1-u ) * ( ( 1-v ) *h[0][0]+v*h[0][1] ) +u* ( ( 1-v ) *h[1][0]+v*h[1][1] );
	h[0][1]= ( 1-u ) * ( ( 1-v ) *h[0][1]+v*h[0][2] ) +u* ( ( 1-v ) *h[1][1]+v*h[1][2] );
	h[0][2]= ( 1-u ) * ( ( 1-v ) *h[0][2]+v*h[0][3] ) +u* ( ( 1-v ) *h[1][2]+v*h[1][3] );
	h[1][0]= ( 1-u ) * ( ( 1-v ) *h[1][0]+v*h[1][1] ) +u* ( ( 1-v ) *h[2][0]+v*h[2][1] );
	h[1][1]= ( 1-u ) * ( ( 1-v ) *h[1][1]+v*h[1][2] ) +u* ( ( 1-v ) *h[2][1]+v*h[2][2] );
	h[1][2]= ( 1-u ) * ( ( 1-v ) *h[1][2]+v*h[1][3] ) +u* ( ( 1-v ) *h[2][2]+v*h[2][3] );
	h[2][0]= ( 1-u ) * ( ( 1-v ) *h[2][0]+v*h[2][1] ) +u* ( ( 1-v ) *h[3][0]+v*h[3][1] );
	h[2][1]= ( 1-u ) * ( ( 1-v ) *h[2][1]+v*h[2][2] ) +u* ( ( 1-v ) *h[3][1]+v*h[3][2] );
	h[2][2]= ( 1-u ) * ( ( 1-v ) *h[2][2]+v*h[2][3] ) +u* ( ( 1-v ) *h[3][2]+v*h[3][3] );
	h[0][0]= ( 1-u ) * ( ( 1-v ) *h[0][0]+v*h[0][1] ) +u* ( ( 1-v ) *h[1][0]+v*h[1][1] );
	h[0][1]= ( 1-u ) * ( ( 1-v ) *h[0][1]+v*h[0][2] ) +u* ( ( 1-v ) *h[1][1]+v*h[1][2] );
	h[1][0]= ( 1-u ) * ( ( 1-v ) *h[1][0]+v*h[1][1] ) +u* ( ( 1-v ) *h[2][0]+v*h[2][1] );
	h[1][1]= ( 1-u ) * ( ( 1-v ) *h[1][1]+v*h[1][2] ) +u* ( ( 1-v ) *h[2][1]+v*h[2][2] );
	return Vektor (	- ( ( ( ( 1-v ) *h[0][1]+v *h[1][1] )- ( ( 1-v ) *h[0][0]+v *h[1][0] ) ) / 1 ),
	                - ( ( ( ( 1-u ) *h[1][0]+u *h[1][1] )- ( ( 1-u ) *h[0][0]+u *h[0][1] ) ) / 1 ),
	                3 );
}
