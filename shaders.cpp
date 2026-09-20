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

static void screenspace_init();

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
	screenspace_init();
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

//==============================================================================
// Screen-space fluid
//
// van der Laan, Green and Sainz, "Screen Space Fluid Rendering with Curvature
// Flow", I3D 2009, with a bilateral blur in place of the curvature flow - it is
// the same idea (smooth the depth buffer, not the geometry) and it does not
// need an iteration count tuned per scene.
//
// Four passes:
//   1  every particle as a sphere, writing its eye-space z and a real depth
//   2  every particle again, additively, for how much water is in front of the
//      pixel - that is what makes thin spray pale and a deep pool dark
//   3  blur the depth, across the silhouette but not across it
//   4  one full-screen triangle: rebuild the eye position from the smoothed
//      depth, take the normal from its screen derivatives, light it
//==============================================================================

static GLuint g_sp_depth = 0, g_sp_thick = 0, g_sp_blur = 0, g_sp_comp = 0;
static bool   g_sp_ok    = false;
static std::string g_sp_error;

static GLuint g_fbo_scene = 0, g_fbo_blur = 0;
static GLuint g_tex_depth = 0, g_tex_thick = 0, g_tex_tmp = 0, g_rb_depth = 0;
static int    g_tw = 0, g_th = 0;

//Nothing is behind the water until something is drawn there; this is the "no
//water here" value the composite pass tests for.
static const float SP_EMPTY = -1.0e9f;

//A sphere, from a square point sprite. Writes eye-space z as colour and the
//sphere's own depth as depth, so the spheres intersect each other properly
//instead of being flat discs.
static const char * SP_SPRITE_VS = R"(#version 400 compatibility
out vec3 vEye;
uniform float uRadius;
uniform vec2  uViewport;
void main()
{
	vEye = ( gl_ModelViewMatrix * gl_Vertex ).xyz;
	gl_Position  = gl_ModelViewProjectionMatrix * gl_Vertex;
	//The projected diameter of a sphere of uRadius at this distance, in pixels.
	gl_PointSize = max ( 1.0, uViewport.y * gl_ProjectionMatrix[1][1] * uRadius / max ( -vEye.z, 0.0001 ) );
}
)";

static const char * SP_DEPTH_FS = R"(#version 400 compatibility
in vec3 vEye;
uniform float uRadius;
out float oEyeZ;
void main()
{
	//gl_PointCoord runs from the top left, so y is flipped to get eye space up.
	vec2  d = vec2 ( gl_PointCoord.x, 1.0-gl_PointCoord.y ) * 2.0 - 1.0;
	float r2 = dot ( d, d );
	if ( r2 > 1.0 )
		discard;
	vec3 p = vEye + uRadius * vec3 ( d, sqrt ( 1.0-r2 ) );
	vec4 c = gl_ProjectionMatrix * vec4 ( p, 1.0 );
	gl_FragDepth = ( c.z/c.w ) * 0.5 + 0.5;
	oEyeZ = p.z;
}
)";

//How much water a ray crosses. No depth test and additive, so everything counts,
//near and far.
static const char * SP_THICK_FS = R"(#version 400 compatibility
in vec3 vEye;
uniform float uRadius;
out float oThickness;
void main()
{
	vec2  d = vec2 ( gl_PointCoord.x, 1.0-gl_PointCoord.y ) * 2.0 - 1.0;
	float r2 = dot ( d, d );
	if ( r2 > 1.0 )
		discard;
	oThickness = 2.0 * uRadius * sqrt ( 1.0-r2 );
}
)";

//One triangle covering the screen, built from the vertex index alone so there is
//no buffer to bind and the current matrices are left untouched for the passes
//that need them.
static const char * SP_FULLSCREEN_VS = R"(#version 400 compatibility
out vec2 vUV;
void main()
{
	vec2 p = vec2 ( ( gl_VertexID << 1 ) & 2, gl_VertexID & 2 );
	vUV = p;
	gl_Position = vec4 ( p*2.0-1.0, 0.0, 1.0 );
}
)";

//Bilateral: a gaussian along one axis, weighted down where the depth jumps, so
//the pool is smoothed and the silhouette against the background is not.
static const char * SP_BLUR_FS = R"(#version 400 compatibility
in vec2 vUV;
uniform sampler2D uDepth;
uniform vec2      uStep;          //one texel along the axis being blurred
uniform float     uProjScale;     //viewport.y * P[1][1]: world units to pixels at z=-1
uniform float     uWorldRadius;   //the splat radius, in world units
uniform float     uMaxRadiusPx;
uniform float     uDepthFalloff;
uniform float     uEmpty;
out float oEyeZ;

//The filter has to be the same size in world units at every distance, or the
//near water is barely smoothed and the far water is smeared into a sheet. So
//the radius is the splat's own projected size, recomputed per pixel from the
//depth that is there.
void main()
{
	float z = texture ( uDepth, vUV ).r;
	if ( z <= uEmpty )
	{
		oEyeZ = z;
		return;
	}
	float n = clamp ( uWorldRadius*uProjScale / max ( -z, 0.0001 ), 1.0, uMaxRadiusPx );
	float sum = 0.0, wsum = 0.0;
	for ( int k=-32; k<=32; ++k )
	{
		float i = float ( k );
		if ( abs ( i ) > n )
			continue;
		float sz = texture ( uDepth, vUV + uStep*i ).r;
		if ( sz <= uEmpty )
			continue;
		float wr = exp ( - ( i*i ) / max ( n*n*0.5, 1.0 ) );
		float dz = ( sz - z ) * uDepthFalloff;
		float wd = exp ( - dz*dz );
		sum  += sz * wr * wd;
		wsum +=      wr * wd;
	}
	oEyeZ = wsum > 0.0 ? sum/wsum : z;
}
)";

//The surface. Everything here comes out of the smoothed depth buffer: the eye
//position by undoing the projection, the normal from that position's screen
//derivatives. The lights and the material are the fixed-function ones, so this
//water is lit by the same two lamps as the rest of the scene.
static const char * SP_COMPOSITE_FS = R"(#version 400 compatibility
in vec2 vUV;
uniform sampler2D uDepth;
uniform sampler2D uThickness;
uniform vec2      uViewport;
uniform float     uWorldRadius;
uniform float     uEmpty;
out vec4 oColour;

vec3 eye_from_depth ( vec2 uv, float z )
{
	vec4 ray = gl_ProjectionMatrixInverse * vec4 ( uv*2.0-1.0, -1.0, 1.0 );
	ray /= ray.w;
	return ray.xyz * ( z / ray.z );
}

void main()
{
	float z = texture ( uDepth, vUV ).r;
	if ( z <= uEmpty )
		discard;

	vec3 p = eye_from_depth ( vUV, z );

	//Derivatives across a silhouette are meaningless, so each axis takes
	//whichever of its two neighbours is nearer in depth. A lone splat - one
	//droplet of spray, a few pixels across - has no valid neighbour on either
	//side, and differencing against the empty background there produced a
	//nonsense normal and a droplet that came out black.
	vec2 texel = 1.0/uViewport;
	float zl = texture ( uDepth, vUV-vec2 ( texel.x,0 ) ).r;
	float zr = texture ( uDepth, vUV+vec2 ( texel.x,0 ) ).r;
	float zd = texture ( uDepth, vUV-vec2 ( 0,texel.y ) ).r;
	float zu = texture ( uDepth, vUV+vec2 ( 0,texel.y ) ).r;
	bool vl = zl > uEmpty, vr = zr > uEmpty, vd = zd > uEmpty, vu = zu > uEmpty;

	vec3 n;
	if ( ( !vl && !vr ) || ( !vd && !vu ) )
	{
		//No surface to take a slope from: face the camera.
		n = vec3 ( 0.0, 0.0, 1.0 );
	}
	else
	{
		vec3 fx = eye_from_depth ( vUV+vec2 ( texel.x,0 ), zr ) - p;
		vec3 bx = p - eye_from_depth ( vUV-vec2 ( texel.x,0 ), zl );
		vec3 dx = ( vr && vl ) ? ( abs ( bx.z ) < abs ( fx.z ) ? bx : fx ) : ( vr ? fx : bx );
		vec3 fy = eye_from_depth ( vUV+vec2 ( 0,texel.y ), zu ) - p;
		vec3 by = p - eye_from_depth ( vUV-vec2 ( 0,texel.y ), zd );
		vec3 dy = ( vu && vd ) ? ( abs ( by.z ) < abs ( fy.z ) ? by : fy ) : ( vu ? fy : by );
		n = normalize ( cross ( dx, dy ) );
		if ( n.z < 0.0 )
			n = -n;
	}

	//How much water the ray crossed, counted in particle diameters so that the
	//look does not change when the scene is scaled or the splat radius tuned.
	float thickness = max ( texture ( uThickness, vUV ).r, 0.0 ) / max ( 2.0*uWorldRadius, 1e-4 );

	//Beer-Lambert: red goes first, which is why deep water is blue-green and a
	//handful of particles of spray is still almost white.
	vec3  absorb = vec3 ( 0.030, 0.010, 0.005 );
	vec3  body   = exp ( -absorb * thickness );
	vec3  deep   = vec3 ( 0.05, 0.25, 0.45 );
	vec3  albedo = mix ( deep, vec3 ( 0.80, 0.92, 1.0 ), body );

	vec3 view = normalize ( -p );

	//Sky light. Without it, deep water goes black in this scene: it absorbs red
	//first, and the only lamp facing the camera here is the red one. Real water
	//is lit mostly by the sky above it, so the ambient is a hemisphere - pale
	//blue from world +z, dark from below - rather than one flat number.
	vec3  up   = normalize ( gl_NormalMatrix * vec3 ( 0.0, 0.0, 1.0 ) );
	float hemi = 0.5 + 0.5*dot ( n, up );
	vec3 rgb   = albedo * mix ( vec3 ( 0.10,0.09,0.08 ), vec3 ( 0.40,0.55,0.75 ), hemi );

	for ( int i=0;i<2;++i )
	{
		vec4 lp = gl_LightSource[i].position;
		vec3 l  = normalize ( lp.w > 0.0 ? lp.xyz - p : lp.xyz );
		float d = max ( dot ( n, l ), 0.0 );
		rgb += albedo * gl_LightSource[i].diffuse.rgb * d;
		if ( d > 0.0 )
		{
			vec3 h = normalize ( l + view );
			rgb += gl_LightSource[i].diffuse.rgb * pow ( max ( dot ( n, h ), 0.0 ), 90.0 );
		}
	}

	//A rim that brightens at grazing angles. Water is mostly reflective there,
	//and without it the silhouette looks like painted plastic.
	float fres = pow ( 1.0 - max ( dot ( n, view ), 0.0 ), 4.0 );
	rgb = mix ( rgb, vec3 ( 0.85,0.92,1.0 ), fres*0.6 );

	//Thin spray is see-through, a pool is not.
	float alpha = clamp ( 1.0 - exp ( -thickness*0.25 ), 0.30, 1.0 );

	vec4 c = gl_ProjectionMatrix * vec4 ( p, 1.0 );
	gl_FragDepth = ( c.z/c.w ) * 0.5 + 0.5;
	oColour = vec4 ( rgb, alpha );
}
)";

static GLuint link_program ( const char * vs_src, const char * fs_src, const char * what )
{
	GLuint vs = glCreateShader ( GL_VERTEX_SHADER );
	GLuint fs = glCreateShader ( GL_FRAGMENT_SHADER );
	GLuint prog = 0;
	if ( compile ( vs, vs_src, what ) && compile ( fs, fs_src, what ) )
	{
		prog = glCreateProgram();
		glAttachShader ( prog, vs );
		glAttachShader ( prog, fs );
		glLinkProgram ( prog );
		GLint ok = 0;
		glGetProgramiv ( prog, GL_LINK_STATUS, &ok );
		if ( !ok )
		{
			char log[4096] = { 0 };
			glGetProgramInfoLog ( prog, sizeof ( log )-1, 0, log );
			g_error = string ( what ) + " link: " + log;
			glDeleteProgram ( prog );
			prog = 0;
		}
	}
	glDeleteShader ( vs );
	glDeleteShader ( fs );
	return prog;
}

static void screenspace_init()
{
	g_sp_ok = false;
	g_sp_error.clear();
	//A fresh context. Whatever the last one held went with it, so the handles
	//are forgotten rather than deleted - deleting them here would name objects
	//in the new context that belong to somebody else.
	g_fbo_scene = g_fbo_blur = 0;
	g_tex_depth = g_tex_thick = g_tex_tmp = g_rb_depth = 0;
	g_tw = g_th = 0;
	g_sp_depth = link_program ( SP_SPRITE_VS, SP_DEPTH_FS, "screen-space depth" );
	g_sp_thick = link_program ( SP_SPRITE_VS, SP_THICK_FS, "screen-space thickness" );
	g_sp_blur  = link_program ( SP_FULLSCREEN_VS, SP_BLUR_FS, "screen-space blur" );
	g_sp_comp  = link_program ( SP_FULLSCREEN_VS, SP_COMPOSITE_FS, "screen-space composite" );
	if ( g_sp_depth && g_sp_thick && g_sp_blur && g_sp_comp )
		g_sp_ok = true;
	else
	{
		g_sp_error = g_error;
		cout << "screen-space fluid off: " << g_sp_error << endl;
	}
}

static GLuint make_r32f ( int w, int h )
{
	GLuint t = 0;
	glGenTextures ( 1, &t );
	glBindTexture ( GL_TEXTURE_2D, t );
	glTexImage2D ( GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, 0 );
	glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
	glTexParameteri ( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
	glBindTexture ( GL_TEXTURE_2D, 0 );
	return t;
}

static bool ensure_targets ( int w, int h )
{
	if ( g_tw == w && g_th == h && g_fbo_scene )
		return true;
	screenspace_free();
	g_tex_depth = make_r32f ( w, h );
	g_tex_thick = make_r32f ( w, h );
	g_tex_tmp   = make_r32f ( w, h );
	glGenRenderbuffers ( 1, &g_rb_depth );
	glBindRenderbuffer ( GL_RENDERBUFFER, g_rb_depth );
	glRenderbufferStorage ( GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h );
	glBindRenderbuffer ( GL_RENDERBUFFER, 0 );

	glGenFramebuffers ( 1, &g_fbo_scene );
	glGenFramebuffers ( 1, &g_fbo_blur );
	glBindFramebuffer ( GL_FRAMEBUFFER, g_fbo_scene );
	glFramebufferTexture2D ( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, g_tex_depth, 0 );
	glFramebufferRenderbuffer ( GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, g_rb_depth );
	const bool ok = glCheckFramebufferStatus ( GL_FRAMEBUFFER ) == GL_FRAMEBUFFER_COMPLETE;
	glBindFramebuffer ( GL_FRAMEBUFFER, 0 );
	if ( !ok )
	{
		g_sp_error = "render targets incomplete";
		screenspace_free();
		return false;
	}
	g_tw = w;
	g_th = h;
	return true;
}

void screenspace_free()
{
	if ( g_fbo_scene ) glDeleteFramebuffers ( 1, &g_fbo_scene );
	if ( g_fbo_blur )  glDeleteFramebuffers ( 1, &g_fbo_blur );
	if ( g_tex_depth ) glDeleteTextures ( 1, &g_tex_depth );
	if ( g_tex_thick ) glDeleteTextures ( 1, &g_tex_thick );
	if ( g_tex_tmp )   glDeleteTextures ( 1, &g_tex_tmp );
	if ( g_rb_depth )  glDeleteRenderbuffers ( 1, &g_rb_depth );
	g_fbo_scene = g_fbo_blur = 0;
	g_tex_depth = g_tex_thick = g_tex_tmp = g_rb_depth = 0;
	g_tw = g_th = 0;
}

bool screenspace_available()
{
	return g_sp_ok;
}

const char * screenspace_error()
{
	return g_sp_error.c_str();
}

//Splat the particles through whichever program is bound.
static void draw_particles ( const float * positions, unsigned long count )
{
	glEnableClientState ( GL_VERTEX_ARRAY );
	glDisableClientState ( GL_NORMAL_ARRAY );
	glDisableClientState ( GL_COLOR_ARRAY );
	glVertexPointer ( 3, GL_FLOAT, 3*sizeof ( float ), positions );
	glDrawArrays ( GL_POINTS, 0, ( GLsizei ) count );
	glDisableClientState ( GL_VERTEX_ARRAY );
}

static void blur_axis ( GLuint src, GLuint dst, float dx, float dy )
{
	glBindFramebuffer ( GL_FRAMEBUFFER, g_fbo_blur );
	glFramebufferTexture2D ( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, dst, 0 );
	glActiveTexture ( GL_TEXTURE0 );
	glBindTexture ( GL_TEXTURE_2D, src );
	glUniform1i ( glGetUniformLocation ( g_sp_blur, "uDepth" ), 0 );
	glUniform2f ( glGetUniformLocation ( g_sp_blur, "uStep" ), dx, dy );
	glDrawArrays ( GL_TRIANGLES, 0, 3 );
}

bool screenspace_render ( const float * positions, unsigned long count,
                          float radius, int smoothing )
{
	if ( !g_sp_ok || count == 0 )
		return false;

	GLint vp[4] = { 0,0,1,1 };
	glGetIntegerv ( GL_VIEWPORT, vp );
	const int w = vp[2], h = vp[3];
	if ( w < 2 || h < 2 || !ensure_targets ( w, h ) )
		return false;

	GLint old_fbo = 0;
	glGetIntegerv ( GL_FRAMEBUFFER_BINDING, &old_fbo );
	glPushAttrib ( GL_ENABLE_BIT | GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT | GL_VIEWPORT_BIT );

	glEnable ( GL_PROGRAM_POINT_SIZE );
	glEnable ( GL_POINT_SPRITE );
	glDisable ( GL_BLEND );
	glDisable ( GL_LIGHTING );

	// 1 - depth
	glBindFramebuffer ( GL_FRAMEBUFFER, g_fbo_scene );
	glFramebufferTexture2D ( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, g_tex_depth, 0 );
	glViewport ( 0,0,w,h );
	glClearColor ( SP_EMPTY, SP_EMPTY, SP_EMPTY, SP_EMPTY );
	glClearDepth ( 1.0 );
	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable ( GL_DEPTH_TEST );
	glDepthMask ( GL_TRUE );
	glUseProgram ( g_sp_depth );
	glUniform1f ( glGetUniformLocation ( g_sp_depth, "uRadius" ), radius );
	glUniform2f ( glGetUniformLocation ( g_sp_depth, "uViewport" ), ( float ) w, ( float ) h );
	draw_particles ( positions, count );

	// 2 - thickness
	glFramebufferTexture2D ( GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, g_tex_thick, 0 );
	glClearColor ( 0,0,0,0 );
	glClear ( GL_COLOR_BUFFER_BIT );
	glDisable ( GL_DEPTH_TEST );
	glDepthMask ( GL_FALSE );
	glEnable ( GL_BLEND );
	glBlendFunc ( GL_ONE, GL_ONE );
	glUseProgram ( g_sp_thick );
	glUniform1f ( glGetUniformLocation ( g_sp_thick, "uRadius" ), radius );
	glUniform2f ( glGetUniformLocation ( g_sp_thick, "uViewport" ), ( float ) w, ( float ) h );
	draw_particles ( positions, count );
	glDisable ( GL_BLEND );

	// 3 - smooth the depth, x then y, as many times as asked
	glUseProgram ( g_sp_blur );
	{
		GLfloat proj[16];
		glGetFloatv ( GL_PROJECTION_MATRIX, proj );
		const float proj_scale = ( float ) h * proj[5];   //P[1][1]
		glUniform1f ( glGetUniformLocation ( g_sp_blur, "uProjScale" ), proj_scale );
		glUniform1f ( glGetUniformLocation ( g_sp_blur, "uWorldRadius" ), radius );
		glUniform1f ( glGetUniformLocation ( g_sp_blur, "uMaxRadiusPx" ), 24.0f );
		glUniform1f ( glGetUniformLocation ( g_sp_blur, "uDepthFalloff" ), 1.0f / ( 2.0f*radius ) );
		glUniform1f ( glGetUniformLocation ( g_sp_blur, "uEmpty" ), SP_EMPTY*0.5f );
	}
	for ( int i=0;i<smoothing;++i )
	{
		blur_axis ( g_tex_depth, g_tex_tmp,   1.0f/w, 0.0f );
		blur_axis ( g_tex_tmp,   g_tex_depth, 0.0f, 1.0f/h );
	}

	// 4 - back to the window, and shade
	glBindFramebuffer ( GL_FRAMEBUFFER, old_fbo );
	glViewport ( vp[0],vp[1],vp[2],vp[3] );
	glUseProgram ( g_sp_comp );
	glActiveTexture ( GL_TEXTURE0 );
	glBindTexture ( GL_TEXTURE_2D, g_tex_depth );
	glActiveTexture ( GL_TEXTURE1 );
	glBindTexture ( GL_TEXTURE_2D, g_tex_thick );
	glActiveTexture ( GL_TEXTURE0 );
	glUniform1i ( glGetUniformLocation ( g_sp_comp, "uDepth" ), 0 );
	glUniform1i ( glGetUniformLocation ( g_sp_comp, "uThickness" ), 1 );
	glUniform2f ( glGetUniformLocation ( g_sp_comp, "uViewport" ), ( float ) w, ( float ) h );
	glUniform1f ( glGetUniformLocation ( g_sp_comp, "uWorldRadius" ), radius );
	glUniform1f ( glGetUniformLocation ( g_sp_comp, "uEmpty" ), SP_EMPTY*0.5f );
	glEnable ( GL_DEPTH_TEST );
	glDepthMask ( GL_TRUE );
	glEnable ( GL_BLEND );
	glBlendFunc ( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glDrawArrays ( GL_TRIANGLES, 0, 3 );

	glUseProgram ( 0 );
	glActiveTexture ( GL_TEXTURE1 );
	glBindTexture ( GL_TEXTURE_2D, 0 );
	glActiveTexture ( GL_TEXTURE0 );
	glBindTexture ( GL_TEXTURE_2D, 0 );
	glPopAttrib();
	return true;
}
