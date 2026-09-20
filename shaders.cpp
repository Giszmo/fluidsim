#include "shaders.h"

#define GL_GLEXT_PROTOTYPES 1
#include <GL/gl.h>
#include <GL/glext.h>

#include <cstdio>
#include <cstring>
#include <string>
#include <iostream>

using namespace std;

static bool        g_ok     = false;
static std::string g_error;
static GLuint      g_prog   = 0;
static GLint       g_loc_seg = -1, g_loc_max = -1, g_loc_viewport = -1, g_loc_lit = -1;

//------------------------------------------------------------------------------
// The shaders
//------------------------------------------------------------------------------

//Straight through: the patch is built in object space, so that the control
//points do not have to be recomputed when the camera moves.
static const char * VS = R"(#version 400 compatibility
out vec3 vPos;
out vec3 vNrm;
out vec4 vCol;
void main()
{
	vPos = gl_Vertex.xyz;
	vNrm = gl_Normal;
	vCol = gl_Color;
}
)";

//One PN triangle patch: ten control points for the position, six for the
//normal. Both are computed once per patch and handed to the evaluation stage.
static const char * TCS = R"(#version 400 compatibility
layout(vertices = 3) out;

in  vec3 vPos[];
in  vec3 vNrm[];
in  vec4 vCol[];

out vec3 tcPos[];
out vec3 tcNrm[];
out vec4 tcCol[];

patch out vec3 b210, b120, b021, b012, b102, b201, b111;
patch out vec3 n110, n011, n101;

uniform vec2  uViewport;
uniform float uPixelsPerSegment;
uniform float uMaxLevel;

//Where a point lands on the screen, in pixels. Used only to decide how finely
//to subdivide, so the near clip plane is handled by simply giving up and asking
//for the finest level on anything that crosses it.
vec2 to_pixels ( vec3 p, out bool ok )
{
	vec4 c = gl_ModelViewProjectionMatrix * vec4 ( p, 1.0 );
	ok = c.w > 0.0001;
	return ok ? ( c.xy/c.w ) * 0.5 * uViewport : vec2 ( 0.0 );
}

//Both triangles that own an edge see the same two endpoints, so they ask for
//the same level and the seam cannot crack.
float edge_level ( vec3 a, vec3 b )
{
	bool oka, okb;
	vec2 pa = to_pixels ( a, oka );
	vec2 pb = to_pixels ( b, okb );
	if ( !oka || !okb )
		return uMaxLevel;
	return clamp ( distance ( pa, pb ) / max ( uPixelsPerSegment, 1.0 ), 1.0, uMaxLevel );
}

vec3 edge_point ( vec3 pi, vec3 pj, vec3 ni )
{
	return ( 2.0*pi + pj - dot ( pj-pi, ni ) * ni ) / 3.0;
}

vec3 edge_normal ( vec3 pi, vec3 pj, vec3 ni, vec3 nj )
{
	vec3  d = pj - pi;
	float v = 2.0 * dot ( d, ni+nj ) / max ( dot ( d, d ), 1e-12 );
	return normalize ( ni + nj - v*d );
}

void main()
{
	tcPos[gl_InvocationID] = vPos[gl_InvocationID];
	tcNrm[gl_InvocationID] = vNrm[gl_InvocationID];
	tcCol[gl_InvocationID] = vCol[gl_InvocationID];

	if ( gl_InvocationID == 0 )
	{
		vec3 p0 = vPos[0], p1 = vPos[1], p2 = vPos[2];
		vec3 n0 = normalize ( vNrm[0] ), n1 = normalize ( vNrm[1] ), n2 = normalize ( vNrm[2] );

		b210 = edge_point ( p0, p1, n0 );
		b120 = edge_point ( p1, p0, n1 );
		b021 = edge_point ( p1, p2, n1 );
		b012 = edge_point ( p2, p1, n2 );
		b102 = edge_point ( p2, p0, n2 );
		b201 = edge_point ( p0, p2, n0 );

		vec3 e = ( b210+b120+b021+b012+b102+b201 ) / 6.0;
		vec3 v = ( p0+p1+p2 ) / 3.0;
		b111  = e + ( e - v ) * 0.5;

		n110 = edge_normal ( p0, p1, n0, n1 );
		n011 = edge_normal ( p1, p2, n1, n2 );
		n101 = edge_normal ( p2, p0, n2, n0 );

		//Outer level i belongs to the edge opposite vertex i.
		gl_TessLevelOuter[0] = edge_level ( p1, p2 );
		gl_TessLevelOuter[1] = edge_level ( p2, p0 );
		gl_TessLevelOuter[2] = edge_level ( p0, p1 );
		gl_TessLevelInner[0] = max ( gl_TessLevelOuter[0],
		                       max ( gl_TessLevelOuter[1], gl_TessLevelOuter[2] ) );
	}
}
)";

//Evaluate the cubic patch for the position and the quadratic one for the
//normal, then put the result where the fixed-function pipeline would have.
static const char * TES = R"(#version 400 compatibility
layout(triangles, fractional_odd_spacing, ccw) in;

in  vec3 tcPos[];
in  vec3 tcNrm[];
in  vec4 tcCol[];

patch in vec3 b210, b120, b021, b012, b102, b201, b111;
patch in vec3 n110, n011, n101;

out vec3 fEyePos;
out vec3 fEyeNrm;
out vec4 fCol;

void main()
{
	float u = gl_TessCoord.x, v = gl_TessCoord.y, w = gl_TessCoord.z;
	vec3 p0 = tcPos[0], p1 = tcPos[1], p2 = tcPos[2];
	vec3 n0 = normalize ( tcNrm[0] ), n1 = normalize ( tcNrm[1] ), n2 = normalize ( tcNrm[2] );

	vec3 p = p0*u*u*u + p1*v*v*v + p2*w*w*w
	       + b210*3.0*u*u*v + b120*3.0*u*v*v
	       + b201*3.0*u*u*w + b021*3.0*v*v*w
	       + b102*3.0*u*w*w + b012*3.0*v*w*w
	       + b111*6.0*u*v*w;

	vec3 n = n0*u*u + n1*v*v + n2*w*w
	       + n110*u*v + n011*v*w + n101*u*w;

	vec4 eye = gl_ModelViewMatrix * vec4 ( p, 1.0 );
	fEyePos  = eye.xyz;
	fEyeNrm  = gl_NormalMatrix * normalize ( n );
	fCol     = tcCol[0]*u + tcCol[1]*v + tcCol[2]*w;
	gl_Position = gl_ModelViewProjectionMatrix * vec4 ( p, 1.0 );
}
)";

//The two lights, the material and the colour come from the fixed-function state
//the rest of the program already sets, so tessellated and untessellated water
//are shaded the same.
static const char * FS = R"(#version 400 compatibility
in vec3 fEyePos;
in vec3 fEyeNrm;
in vec4 fCol;

uniform bool uLit;

void main()
{
	if ( !uLit )
	{
		gl_FragColor = fCol;
		return;
	}
	vec3 n = normalize ( fEyeNrm );
	if ( !gl_FrontFacing )
		n = -n;
	vec3 view = normalize ( -fEyePos );
	vec3 rgb  = fCol.rgb * 0.2;
	for ( int i=0;i<2;++i )
	{
		vec4 lp = gl_LightSource[i].position;
		vec3 l  = normalize ( lp.w > 0.0 ? lp.xyz - fEyePos : lp.xyz );
		float d = max ( dot ( n, l ), 0.0 );
		rgb += fCol.rgb * gl_LightSource[i].diffuse.rgb * d;
		if ( d > 0.0 )
		{
			vec3 h = normalize ( l + view );
			rgb += gl_FrontMaterial.specular.rgb * gl_LightSource[i].diffuse.rgb
			     * pow ( max ( dot ( n, h ), 0.0 ), max ( gl_FrontMaterial.shininess, 1.0 ) );
		}
	}
	gl_FragColor = vec4 ( rgb, fCol.a );
}
)";

//------------------------------------------------------------------------------
// Compiling
//------------------------------------------------------------------------------

static bool compile ( GLuint shader, const char * src, const char * what )
{
	glShaderSource ( shader, 1, &src, 0 );
	glCompileShader ( shader );
	GLint ok = 0;
	glGetShaderiv ( shader, GL_COMPILE_STATUS, &ok );
	if ( ok )
		return true;
	char log[4096] = { 0 };
	glGetShaderInfoLog ( shader, sizeof ( log )-1, 0, log );
	g_error = string ( what ) + ": " + log;
	return false;
}

//OpenGL 4.0 or the extension. Asked of the driver rather than assumed, because
//the whole point is that the program still runs where there is neither.
static bool has_tessellation()
{
	GLint major = 0, minor = 0;
	glGetIntegerv ( GL_MAJOR_VERSION, &major );
	glGetIntegerv ( GL_MINOR_VERSION, &minor );
	if ( glGetError() == GL_NO_ERROR && ( major > 4 || ( major == 4 && minor >= 0 ) ) )
		return true;
	GLint n = 0;
	glGetIntegerv ( GL_NUM_EXTENSIONS, &n );
	for ( GLint i=0;i<n;++i )
	{
		const char * e = ( const char * ) glGetStringi ( GL_EXTENSIONS, i );
		if ( e && !strcmp ( e, "GL_ARB_tessellation_shader" ) )
			return true;
	}
	return false;
}

//Called once for every GL context the program makes, so everything is rebuilt
//from scratch each time - a window that has been recreated does not keep a
//program object belonging to a context that is gone.
bool surface_shaders_init()
{
	g_ok = false;
	g_prog = 0;
	g_error.clear();

	if ( !has_tessellation() )
	{
		g_error = "no OpenGL 4 tessellation on this driver";
		return false;
	}

	GLuint vs = glCreateShader ( GL_VERTEX_SHADER );
	GLuint cs = glCreateShader ( GL_TESS_CONTROL_SHADER );
	GLuint es = glCreateShader ( GL_TESS_EVALUATION_SHADER );
	GLuint fs = glCreateShader ( GL_FRAGMENT_SHADER );
	bool built = compile ( vs, VS, "vertex shader" )
	          && compile ( cs, TCS, "tessellation control shader" )
	          && compile ( es, TES, "tessellation evaluation shader" )
	          && compile ( fs, FS, "fragment shader" );
	if ( built )
	{
		g_prog = glCreateProgram();
		glAttachShader ( g_prog, vs );
		glAttachShader ( g_prog, cs );
		glAttachShader ( g_prog, es );
		glAttachShader ( g_prog, fs );
		glLinkProgram ( g_prog );
		GLint ok = 0;
		glGetProgramiv ( g_prog, GL_LINK_STATUS, &ok );
		if ( !ok )
		{
			char log[4096] = { 0 };
			glGetProgramInfoLog ( g_prog, sizeof ( log )-1, 0, log );
			g_error = string ( "link: " ) + log;
			glDeleteProgram ( g_prog );
			g_prog = 0;
			built = false;
		}
	}
	glDeleteShader ( vs );
	glDeleteShader ( cs );
	glDeleteShader ( es );
	glDeleteShader ( fs );

	if ( built )
	{
		g_loc_seg      = glGetUniformLocation ( g_prog, "uPixelsPerSegment" );
		g_loc_max      = glGetUniformLocation ( g_prog, "uMaxLevel" );
		g_loc_viewport = glGetUniformLocation ( g_prog, "uViewport" );
		g_loc_lit      = glGetUniformLocation ( g_prog, "uLit" );
		g_ok = true;
	}
	else
	{
		cout << "surface tessellation off: " << g_error << endl;
	}
	return g_ok;
}

bool surface_tessellation_available()
{
	return g_ok;
}

const char * surface_shaders_error()
{
	return g_error.c_str();
}

void surface_tessellation_bind ( float pixels_per_segment, float max_level )
{
	if ( !g_ok )
		return;
	GLint vp[4] = { 0,0,1,1 };
	glGetIntegerv ( GL_VIEWPORT, vp );
	glUseProgram ( g_prog );
	glUniform1f ( g_loc_seg, pixels_per_segment );
	glUniform1f ( g_loc_max, max_level );
	glUniform2f ( g_loc_viewport, ( float ) vp[2], ( float ) vp[3] );
	glUniform1i ( g_loc_lit, glIsEnabled ( GL_LIGHTING ) ? 1 : 0 );
	glPatchParameteri ( GL_PATCH_VERTICES, 3 );
}

void surface_tessellation_unbind()
{
	if ( g_ok )
		glUseProgram ( 0 );
}
