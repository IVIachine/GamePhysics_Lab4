/*
	Copyright 2011-2018 Daniel S. Buckstein

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/

/*
	animal3D SDK: Minimal 3D Animation Framework
	animal3D Math (A3DM) SDK
	By Daniel S. Buckstein

	a3vector.h
	Declarations for vector functions.
*/

#ifndef __ANIMAL3D_A3DM_VECTOR_H
#define __ANIMAL3D_A3DM_VECTOR_H


#include "animal3D/a3/a3config.h"
#include "animal3D/a3/a3macros.h"
#include "animal3D/a3/a3types_integer.h"
#include "animal3D/a3/a3types_real.h"


A3_BEGIN_DECL


//-----------------------------------------------------------------------------
// simple wrapper data structures
// these data structures ultimately serve as containers or aliases for simple 
//	arrays of a3real numbers/floats; all of the vector functions and operations 
//	below can be used to process the variable 'v' contained in each structure, 
//	which is ultimately just an array, or some other user-defined containers 
//	of vector data that hold an array

#ifndef __cplusplus
typedef union a3vec2 a3vec2;
typedef union a3vec3 a3vec3;
typedef union a3vec4 a3vec4;
#endif	// !__cplusplus

// 2D vector
union a3vec2
{
	// components
	struct { a3real x, y; };
	struct { a3real s, t; };
	struct { a3real r, g; };
	struct { a3real v0, v1; };

	// array (INPUT FOR OPERATIONS)
	a3real2 v;
};

// 3D vector
union a3vec3
{
	// components
	struct { a3real x, y, z; };
	struct { a3real s, t, p; };
	struct { a3real r, g, b; };
	struct { a3real v0, v1, v2; };

	// array (INPUT FOR OPERATIONS)
	a3real3 v;

	// partial swizzle
	a3vec2 xy;
	a3vec2 st;
	a3vec2 rg;
};

// 4D vector
union a3vec4
{
	// components
	struct { a3real x, y, z, w; };
	struct { a3real s, t, p, q; };
	struct { a3real r, g, b, a; };
	struct { a3real v0, v1, v2, v3; };

	// array (INPUT FOR OPERATIONS)
	a3real4 v;

	// partial swizzle
	a3vec3 xyz;
	a3vec3 stp;
	a3vec3 rgb;
	a3vec2 xy;
	a3vec2 st;
	a3vec2 rg;
};


//-----------------------------------------------------------------------------
// function declarations using array types
// for all functions that use other inputs to modify a vector, the 
//	result/output appears as the first parameter 'v_out'
// functions that modify a vector directly have a parameter 'v_inout'
// any vector used as input is labelled 'v' and may have some suffix
//	to give it a better description when needed

// general function descriptions (to avoid redundant comments): 
// --set: load values into a vector array
// --length squared: returns the squared length of a vector
// --length: returns the magnitude of a vector
//	*note: do not do length * length to get squared length as this is slow and 
//	redundant; use the squared length function to do this as it is optimized!
// --get unit/normalize: get the unit-length vector that points in the same 
//	direction of the input, or normalize the input to make it unit-length
//	*note: these return inverse length to avoid needing to re-calculate it
// --get negative/negate: get the negated input or negate the input
// --dot: the dot product of two vectors (x0*x1 + y0*y1 ...)
// --distance squared: returns the squared distance between two points
// --distance: returns the distance between two points
//	*note: see note for length, same idea
// --projection ratio: the ratio used to perform vector projection
// --projection: perform vector projection onto the base vector; essentially 
//	flattens the input against the base, so that the result has the same 
//	direction but a different length
//	*note: returns projection ratio to avoid needing to re-calculate
// --sum/add: compute the sum of two vectors; add stores result in left input
// --diff/sub: compute the difference, same idea for sub
// --product/mul: multiply vector by scalar s
// --quotient/div: divide vector by scalar s
// --product/mul components: multiply respective components of vectors
// --quotient/div components: divide respective components of vectors
// --all interpolation algorithms implemented for a3real numbers are also 
//	implemented for vectors (Catmull-Rom, Hermite, tables...)
// --Gram-Schmidt orthonormalization: using a base vector as reference, modify 
//	vector(s) so that they are perpendicular (orthogonal) to all the others
// --triangle area/doubled/squared/normal: calculate the area of a triangle 
//	formed by 3 points with modifiers; 3D has option to calculate normal
// --cross product: calculate a vector that is perpendicular to both inputs
//	*note: flipping the input yields a negated result (non-commutative)

//-----------------------------------------------------------------------------
// 2D VECTOR FUNCTIONS

// create 2D vector as array
///
A3_INLINE a3real2r a3real2Set(a3real2p v_out, const a3real x, const a3real y);
A3_INLINE a3real2r a3real2SetReal2(a3real2p v_out, const a3real2p v);
A3_INLINE a3real2r a3real2SetReal3(a3real2p v_out, const a3real3p v);
A3_INLINE a3real2r a3real2SetReal4(a3real2p v_out, const a3real4p v);

// unary operations
///
A3_INLINE a3real a3real2LengthSquared(const a3real2p v);
A3_INLINE a3real a3real2LengthSquaredInverse(const a3real2p v);
A3_INLINE a3real a3real2Length(const a3real2p v);
A3_INLINE a3real a3real2LengthInverse(const a3real2p v);
A3_INLINE a3real2r a3real2GetUnit(a3real2p v_out, const a3real2p v);
A3_INLINE a3real2r a3real2Normalize(a3real2p v_inout);
A3_INLINE a3real2r a3real2GetUnitInvLength(a3real2p v_out, const a3real2p v, a3real *invLength_out);
A3_INLINE a3real2r a3real2NormalizeGetInvLength(a3real2p v_inout, a3real *invLength_out);
A3_INLINE a3real2r a3real2GetNegative(a3real2p v_out, const a3real2p v);
A3_INLINE a3real2r a3real2Negate(a3real2p v_inout);

// binary operations
///
A3_INLINE a3real a3real2Dot(const a3real2p vL, const a3real2p vR);
A3_INLINE a3real a3real2DistanceSquared(const a3real2p v0, const a3real2p v1);
A3_INLINE a3real a3real2Distance(const a3real2p v0, const a3real2p v1);
A3_INLINE a3real a3real2ProjRatio(const a3real2p v, const a3real2p vBase);
A3_INLINE a3real2r a3real2Projected(a3real2p v_out, const a3real2p v, const a3real2p vBase);
A3_INLINE a3real2r a3real2Proj(a3real2p v_inout, const a3real2p vBase);
A3_INLINE a3real2r a3real2ProjectedGetRatio(a3real2p v_out, const a3real2p v, const a3real2p vBase, a3real *ratio_out);
A3_INLINE a3real2r a3real2ProjGetRatio(a3real2p v_inout, const a3real2p vBase, a3real *ratio_out);
A3_INLINE a3real2r a3real2Sum(a3real2p v_out, const a3real2p vL, const a3real2p vR);
A3_INLINE a3real2r a3real2Diff(a3real2p v_out, const a3real2p vL, const a3real2p vR);
A3_INLINE a3real2r a3real2ProductS(a3real2p v_out, const a3real2p v, const a3real s);
A3_INLINE a3real2r a3real2QuotientS(a3real2p v_out, const a3real2p v, const a3real s);
A3_INLINE a3real2r a3real2ProductComp(a3real2p v_out, const a3real2p vL, const a3real2p vR);
A3_INLINE a3real2r a3real2QuotientComp(a3real2p v_out, const a3real2p vL, const a3real2p vR);
A3_INLINE a3real2r a3real2Add(a3real2p vL_inout, const a3real2p vR);
A3_INLINE a3real2r a3real2Sub(a3real2p vL_inout, const a3real2p vR);
A3_INLINE a3real2r a3real2MulS(a3real2p v_inout, const a3real s);
A3_INLINE a3real2r a3real2DivS(a3real2p v_inout, const a3real s);
A3_INLINE a3real2r a3real2MulComp(a3real2p vL_inout, const a3real2p vR);
A3_INLINE a3real2r a3real2DivComp(a3real2p vL_inout, const a3real2p vR);

// interpolation algorithms
///
A3_INLINE a3real2r a3real2Lerp(a3real2p v_out, const a3real2p v0, const a3real2p v1, const a3real param);
A3_INLINE a3real2r a3real2NLerp(a3real2p v_out, const a3real2p v0, const a3real2p v1, const a3real param);
A3_INLINE a3real2r a3real2Bilerp(a3real2p v_out, const a3real2p v00, const a3real2p v01, const a3real2p v10, const a3real2p v11, const a3real param0, const a3real param1);
A3_INLINE a3real2r a3real2Trilerp(a3real2p v_out, const a3real2p v000, const a3real2p v001, const a3real2p v010, const a3real2p v011, const a3real2p v100, const a3real2p v101, const a3real2p v110, const a3real2p v111, const a3real param0, const a3real param1, const a3real param2);
A3_INLINE a3real2r a3real2CatmullRom(a3real2p v_out, const a3real2p vPrev, const a3real2p v0, const a3real2p v1, const a3real2p vNext, const a3real param);
A3_INLINE a3real2r a3real2HermiteControl(a3real2p v_out, const a3real2p v0, const a3real2p v1, const a3real2p vControl0, const a3real2p vControl1, const a3real param);
A3_INLINE a3real2r a3real2HermiteTangent(a3real2p v_out, const a3real2p v0, const a3real2p v1, const a3real2p vTangent0, const a3real2p vTangent1, const a3real param);
A3_INLINE a3real2r a3real2Bezier0(a3real2p v_out, const a3real2p v0, const a3real param);
A3_INLINE a3real2r a3real2Bezier1(a3real2p v_out, const a3real2p v0, const a3real2p v1, const a3real param);
A3_INLINE a3real2r a3real2Bezier2(a3real2p v_out, const a3real2p v0, const a3real2p v1, const a3real2p v2, const a3real param);
A3_INLINE a3real2r a3real2Bezier3(a3real2p v_out, const a3real2p v0, const a3real2p v1, const a3real2p v2, const a3real2p v3, const a3real param);
A3_INLINE a3real2r a3real2BezierN(a3real2p v_out, const a3real2 *vValues, const a3real param, a3index order_N);
A3_INLINE a3real2r a3real2Slerp(a3real2p v_out, const a3real2p v0, const a3real2p v1, const a3real param);
A3_INLINE a3real2r a3real2SlerpUnit(a3real2p v_out, const a3real2p v0, const a3real2p v1, const a3real param);
A3_INLINE a3real a3real2CalculateArcLengthCatmullRom(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2p vPrev, const a3real2p v0, const a3real2p v1, const a3real2p vNext);
A3_INLINE a3real a3real2CalculateArcLengthHermiteControl(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2p v0, const a3real2p v1, const a3real2p vControl0, const a3real2p vControl1);
A3_INLINE a3real a3real2CalculateArcLengthHermiteTangent(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2p v0, const a3real2p v1, const a3real2p vTangent0, const a3real2p vTangent1);
A3_INLINE a3real a3real2CalculateArcLengthBezier0(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2p v0);
A3_INLINE a3real a3real2CalculateArcLengthBezier1(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2p v0, const a3real2p v1);
A3_INLINE a3real a3real2CalculateArcLengthBezier2(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2p v0, const a3real2p v1, const a3real2p v2);
A3_INLINE a3real a3real2CalculateArcLengthBezier3(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2p v0, const a3real2p v1, const a3real2p v2, const a3real2p v3);
A3_INLINE a3real a3real2CalculateArcLengthBezierN(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2 *vValues, a3index order_N);
A3_INLINE a3real a3real2CalculateArcLengthSlerp(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2p v0, const a3real2p v1);
A3_INLINE a3real a3real2CalculateArcLengthSlerpUnit(a3real2 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real2p v0, const a3real2p v1);

// Gram-Schmidt orthonormalization
///
A3_INLINE a3real2r a3real2GramSchmidtOrtho(a3real2p v_out, const a3real2p v, const a3real2p vBase);
A3_INLINE a3real2r a3real2GramSchmidt(a3real2p v_inout, const a3real2p vBase);

// 2D cross product: get the signed length of the perpendicular axis (Z)
///
A3_INLINE a3real a3real2CrossZ(const a3real2p vL, const a3real2p vR);

// area of a triangle
///
A3_INLINE a3real a3real2TriangleAreaDoubled(const a3real2p v0, const a3real2p v1, const a3real2p v2);
A3_INLINE a3real a3real2TriangleAreaSquared(const a3real2p v0, const a3real2p v1, const a3real2p v2);
A3_INLINE a3real a3real2TriangleArea(const a3real2p v0, const a3real2p v1, const a3real2p v2);

// point is in triangle (barycentric coordinates)
///
A3_INLINE a3boolean a3real2PointIsInTriangle(const a3real2p p, const a3real2p v0, const a3real2p v1, const a3real2p v2);
A3_INLINE a3boolean a3real2PointIsInTriangleBarycentric(const a3real2p p, const a3real2p v0, const a3real2p v1, const a3real2p v2, a3real *param0, a3real *param1, a3real *param2);


//-----------------------------------------------------------------------------
// 3D VECTOR FUNCTIONS

// create 3D vector as array
///
A3_INLINE a3real3r a3real3Set(a3real3p v_out, const a3real x, const a3real y, const a3real z);
A3_INLINE a3real3r a3real3SetReal2(a3real3p v_out, const a3real2p v);
A3_INLINE a3real3r a3real3SetReal2Z(a3real3p v_out, const a3real2p v, const a3real z);
A3_INLINE a3real3r a3real3SetReal3(a3real3p v_out, const a3real3p v);
A3_INLINE a3real3r a3real3SetReal4(a3real3p v_out, const a3real4p v);

// unary operations
///
A3_INLINE a3real a3real3LengthSquared(const a3real3p v);
A3_INLINE a3real a3real3LengthSquaredInverse(const a3real3p v);
A3_INLINE a3real a3real3Length(const a3real3p v);
A3_INLINE a3real a3real3LengthInverse(const a3real3p v);
A3_INLINE a3real3r a3real3GetUnit(a3real3p v_out, const a3real3p v);
A3_INLINE a3real3r a3real3Normalize(a3real3p v_inout);
A3_INLINE a3real3r a3real3GetUnitInvLength(a3real3p v_out, const a3real3p v, a3real *invLength_out);
A3_INLINE a3real3r a3real3NormalizeGetInvLength(a3real3p v_inout, a3real *invLength_out);
A3_INLINE a3real3r a3real3GetNegative(a3real3p v_out, const a3real3p v);
A3_INLINE a3real3r a3real3Negate(a3real3p v_inout);

// binary operations
///
A3_INLINE a3real a3real3Dot(const a3real3p vL, const a3real3p vR);
A3_INLINE a3real a3real3DistanceSquared(const a3real3p v0, const a3real3p v1);
A3_INLINE a3real a3real3Distance(const a3real3p v0, const a3real3p v1);
A3_INLINE a3real a3real3ProjRatio(const a3real3p v, const a3real3p vBase);
A3_INLINE a3real3r a3real3Projected(a3real3p v_out, const a3real3p v, const a3real3p vBase);
A3_INLINE a3real3r a3real3Proj(a3real3p v_inout, const a3real3p vBase);
A3_INLINE a3real3r a3real3ProjectedGetRatio(a3real3p v_out, const a3real3p v, const a3real3p vBase, a3real *ratio_out);
A3_INLINE a3real3r a3real3ProjGetRatio(a3real3p v_inout, const a3real3p vBase, a3real *ratio_out);
A3_INLINE a3real3r a3real3Sum(a3real3p v_out, const a3real3p vL, const a3real3p vR);
A3_INLINE a3real3r a3real3Diff(a3real3p v_out, const a3real3p vL, const a3real3p vR);
A3_INLINE a3real3r a3real3ProductS(a3real3p v_out, const a3real3p v, const a3real s);
A3_INLINE a3real3r a3real3QuotientS(a3real3p v_out, const a3real3p v, const a3real s);
A3_INLINE a3real3r a3real3ProductComp(a3real3p v_out, const a3real3p vL, const a3real3p vR);
A3_INLINE a3real3r a3real3QuotientComp(a3real3p v_out, const a3real3p vL, const a3real3p vR);
A3_INLINE a3real3r a3real3Add(a3real3p vL_inout, const a3real3p vR);
A3_INLINE a3real3r a3real3Sub(a3real3p vL_inout, const a3real3p vR);
A3_INLINE a3real3r a3real3MulS(a3real3p v_inout, const a3real s);
A3_INLINE a3real3r a3real3DivS(a3real3p v_inout, const a3real s);
A3_INLINE a3real3r a3real3MulComp(a3real3p vL_inout, const a3real3p vR);
A3_INLINE a3real3r a3real3DivComp(a3real3p vL_inout, const a3real3p vR);

// interpolation algorithms
///
A3_INLINE a3real3r a3real3Lerp(a3real3p v_out, const a3real3p v0, const a3real3p v1, const a3real param);
A3_INLINE a3real3r a3real3NLerp(a3real3p v_out, const a3real3p v0, const a3real3p v1, const a3real param);
A3_INLINE a3real3r a3real3Bilerp(a3real3p v_out, const a3real3p v00, const a3real3p v01, const a3real3p v10, const a3real3p v11, const a3real param0, const a3real param1);
A3_INLINE a3real3r a3real3Trilerp(a3real3p v_out, const a3real3p v000, const a3real3p v001, const a3real3p v010, const a3real3p v011, const a3real3p v100, const a3real3p v101, const a3real3p v110, const a3real3p v111, const a3real param0, const a3real param1, const a3real param2);
A3_INLINE a3real3r a3real3CatmullRom(a3real3p v_out, const a3real3p vPrev, const a3real3p v0, const a3real3p v1, const a3real3p vNext, const a3real param);
A3_INLINE a3real3r a3real3HermiteControl(a3real3p v_out, const a3real3p v0, const a3real3p v1, const a3real3p vControl0, const a3real3p vControl1, const a3real param);
A3_INLINE a3real3r a3real3HermiteTangent(a3real3p v_out, const a3real3p v0, const a3real3p v1, const a3real3p vTangent0, const a3real3p vTangent1, const a3real param);
A3_INLINE a3real3r a3real3Bezier0(a3real3p v_out, const a3real3p v0, const a3real param);
A3_INLINE a3real3r a3real3Bezier1(a3real3p v_out, const a3real3p v0, const a3real3p v1, const a3real param);
A3_INLINE a3real3r a3real3Bezier2(a3real3p v_out, const a3real3p v0, const a3real3p v1, const a3real3p v2, const a3real param);
A3_INLINE a3real3r a3real3Bezier3(a3real3p v_out, const a3real3p v0, const a3real3p v1, const a3real3p v2, const a3real3p v3, const a3real param);
A3_INLINE a3real3r a3real3BezierN(a3real3p v_out, const a3real3 *vValues, const a3real param, a3index order_N);
A3_INLINE a3real3r a3real3Slerp(a3real3p v_out, const a3real3p v0, const a3real3p v1, const a3real param);
A3_INLINE a3real3r a3real3SlerpUnit(a3real3p v_out, const a3real3p v0, const a3real3p v1, const a3real param);
A3_INLINE a3real a3real3CalculateArcLengthCatmullRom(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3p vPrev, const a3real3p v0, const a3real3p v1, const a3real3p vNext);
A3_INLINE a3real a3real3CalculateArcLengthHermiteControl(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3p v0, const a3real3p v1, const a3real3p vControl0, const a3real3p vControl1);
A3_INLINE a3real a3real3CalculateArcLengthHermiteTangent(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3p v0, const a3real3p v1, const a3real3p vTangent0, const a3real3p vTangent1);
A3_INLINE a3real a3real3CalculateArcLengthBezier0(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3p v0);
A3_INLINE a3real a3real3CalculateArcLengthBezier1(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3p v0, const a3real3p v1);
A3_INLINE a3real a3real3CalculateArcLengthBezier2(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real a3real3CalculateArcLengthBezier3(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3p v0, const a3real3p v1, const a3real3p v2, const a3real3p v3);
A3_INLINE a3real a3real3CalculateArcLengthBezierN(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3 *vValues, a3index order_N);
A3_INLINE a3real a3real3CalculateArcLengthSlerp(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3p v0, const a3real3p v1);
A3_INLINE a3real a3real3CalculateArcLengthSlerpUnit(a3real3 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real3p v0, const a3real3p v1);

// Gram-Schmidt orthonormalization
///
A3_INLINE a3real3r a3real3GramSchmidtOrtho(a3real3p v_out, const a3real3p v, const a3real3p vBase);
A3_INLINE a3real3r a3real3GramSchmidtOrtho2(a3real3p v0_out, a3real3p v1_out, const a3real3p v0, const a3real3p v1, const a3real3p vBase);
A3_INLINE a3real3r a3real3GramSchmidt(a3real3p v_inout, const a3real3p vBase);
A3_INLINE a3real3r a3real3GramSchmidt2(a3real3p v0_inout, a3real3p v1_inout, const a3real3p vBase);

// cross product and applications (3D only)
///
A3_INLINE a3real3r a3real3Cross(a3real3p v_out, const a3real3p vL, const a3real3p vR);
A3_INLINE a3real3r a3real3CrossUnit(a3real3p v_out, const a3real3p vL, const a3real3p vR);
A3_INLINE a3real3r a3real3TriangleNormal(a3real3p n_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real3r a3real3TriangleNormalUnit(a3real3p n_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);

// normal and area of triangle
///
A3_INLINE a3real a3real3TriangleAreaDoubledNormal(a3real3p n_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real a3real3TriangleAreaSquaredNormal(a3real3p n_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real a3real3TriangleAreaNormal(a3real3p n_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real a3real3TriangleAreaDoubledNormalUnit(a3real3p n_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real a3real3TriangleAreaSquaredNormalUnit(a3real3p n_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real a3real3TriangleAreaNormalUnit(a3real3p n_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);

// area of a triangle
///
A3_INLINE a3real a3real3TriangleAreaDoubled(const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real a3real3TriangleAreaSquared(const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real a3real3TriangleArea(const a3real3p v0, const a3real3p v1, const a3real3p v2);

// point is in triangle (barycentric coordinates)
///
A3_INLINE a3boolean a3real3PointIsInTriangle(const a3real3p p, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3boolean a3real3PointIsInTriangleBarycentric(const a3real3p p, const a3real3p v0, const a3real3p v1, const a3real3p v2, a3real *param0, a3real *param1, a3real *param2);


//-----------------------------------------------------------------------------
// 4D VECTOR FUNCTIONS

// create 4D vector as array
///
A3_INLINE a3real4r a3real4Set(a3real4p v_out, const a3real x, const a3real y, const a3real z, const a3real w);
A3_INLINE a3real4r a3real4SetReal2(a3real4p v_out, const a3real2p v);
A3_INLINE a3real4r a3real4SetReal2Z(a3real4p v_out, const a3real2p v, const a3real z);
A3_INLINE a3real4r a3real4SetReal2ZW(a3real4p v_out, const a3real2p v, const a3real z, const a3real w);
A3_INLINE a3real4r a3real4SetReal3(a3real4p v_out, const a3real3p v);
A3_INLINE a3real4r a3real4SetReal3W(a3real4p v_out, const a3real3p v, const a3real w);
A3_INLINE a3real4r a3real4SetReal4(a3real4p v_out, const a3real4p v);

// unary operations
///
A3_INLINE a3real a3real4LengthSquared(const a3real4p v);
A3_INLINE a3real a3real4LengthSquaredInverse(const a3real4p v);
A3_INLINE a3real a3real4Length(const a3real4p v);
A3_INLINE a3real a3real4LengthInverse(const a3real4p v);
A3_INLINE a3real4r a3real4GetUnit(a3real4p v_out, const a3real4p v);
A3_INLINE a3real4r a3real4Normalize(a3real4p v_inout);
A3_INLINE a3real4r a3real4GetUnitInvLength(a3real4p v_out, const a3real4p v, a3real *invLen_out);
A3_INLINE a3real4r a3real4NormalizeGetInvLength(a3real4p v_inout, a3real *invLen_out);
A3_INLINE a3real4r a3real4GetNegative(a3real4p v_out, const a3real4p v);
A3_INLINE a3real4r a3real4Negate(a3real4p v_inout);

// binary operations
///
A3_INLINE a3real a3real4Dot(const a3real4p vL, const a3real4p vR);
A3_INLINE a3real a3real4DistanceSquared(const a3real4p v0, const a3real4p v1);
A3_INLINE a3real a3real4Distance(const a3real4p v0, const a3real4p v1);
A3_INLINE a3real a3real4ProjRatio(const a3real4p v, const a3real4p vBase);
A3_INLINE a3real4r a3real4Projected(a3real4p v_out, const a3real4p v, const a3real4p vBase);
A3_INLINE a3real4r a3real4Proj(a3real4p v_inout, const a3real4p vBase);
A3_INLINE a3real4r a3real4ProjectedGetRatio(a3real4p v_out, const a3real4p v, const a3real4p vBase, a3real *ratio_out);
A3_INLINE a3real4r a3real4ProjGetRatio(a3real4p v_inout, const a3real4p vBase, a3real *ratio_out);
A3_INLINE a3real4r a3real4Sum(a3real4p v_out, const a3real4p vL, const a3real4p vR);
A3_INLINE a3real4r a3real4Diff(a3real4p v_out, const a3real4p vL, const a3real4p vR);
A3_INLINE a3real4r a3real4ProductS(a3real4p v_out, const a3real4p v, const a3real s);
A3_INLINE a3real4r a3real4QuotientS(a3real4p v_out, const a3real4p v, const a3real s);
A3_INLINE a3real4r a3real4ProductComp(a3real4p v_out, const a3real4p vL, const a3real4p vR);
A3_INLINE a3real4r a3real4QuotientComp(a3real4p v_out, const a3real4p vL, const a3real4p vR);
A3_INLINE a3real4r a3real4Add(a3real4p vL_inout, const a3real4p vR);
A3_INLINE a3real4r a3real4Sub(a3real4p vL_inout, const a3real4p vR);
A3_INLINE a3real4r a3real4MulS(a3real4p v_inout, const a3real s);
A3_INLINE a3real4r a3real4DivS(a3real4p v_inout, const a3real s);
A3_INLINE a3real4r a3real4MulComp(a3real4p vL_inout, const a3real4p vR);
A3_INLINE a3real4r a3real4DivComp(a3real4p vL_inout, const a3real4p vR);

// interpolation algorithms
///
A3_INLINE a3real4r a3real4Lerp(a3real4p v_out, const a3real4p v0, const a3real4p v1, const a3real param);
A3_INLINE a3real4r a3real4NLerp(a3real4p v_out, const a3real4p v0, const a3real4p v1, const a3real param);
A3_INLINE a3real4r a3real4Bilerp(a3real4p v_out, const a3real4p v00, const a3real4p v01, const a3real4p v10, const a3real4p v11, const a3real param0, const a3real param1);
A3_INLINE a3real4r a3real4Trilerp(a3real4p v_out, const a3real4p v000, const a3real4p v001, const a3real4p v010, const a3real4p v011, const a3real4p v100, const a3real4p v101, const a3real4p v110, const a3real4p v111, const a3real param0, const a3real param1, const a3real param2);
A3_INLINE a3real4r a3real4CatmullRom(a3real4p v_out, const a3real4p vPrev, const a3real4p v0, const a3real4p v1, const a3real4p vNext, const a3real param);
A3_INLINE a3real4r a3real4HermiteControl(a3real4p v_out, const a3real4p v0, const a3real4p v1, const a3real4p vControl0, const a3real4p vControl1, const a3real param);
A3_INLINE a3real4r a3real4HermiteTangent(a3real4p v_out, const a3real4p v0, const a3real4p v1, const a3real4p vTangent0, const a3real4p vTangent1, const a3real param);
A3_INLINE a3real4r a3real4Bezier0(a3real4p v_out, const a3real4p v0, const a3real param);
A3_INLINE a3real4r a3real4Bezier1(a3real4p v_out, const a3real4p v0, const a3real4p v1, const a3real param);
A3_INLINE a3real4r a3real4Bezier2(a3real4p v_out, const a3real4p v0, const a3real4p v1, const a3real4p v2, const a3real param);
A3_INLINE a3real4r a3real4Bezier3(a3real4p v_out, const a3real4p v0, const a3real4p v1, const a3real4p v2, const a3real4p v3, const a3real param);
A3_INLINE a3real4r a3real4BezierN(a3real4p v_out, const a3real4 *vValues, const a3real param, a3index order_N);
A3_INLINE a3real4r a3real4Slerp(a3real4p v_out, const a3real4p v0, const a3real4p v1, const a3real param);
A3_INLINE a3real4r a3real4SlerpUnit(a3real4p v_out, const a3real4p v0, const a3real4p v1, const a3real param);
A3_INLINE a3real a3real4CalculateArcLengthCatmullRom(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4p vPrev, const a3real4p v0, const a3real4p v1, const a3real4p vNext);
A3_INLINE a3real a3real4CalculateArcLengthHermiteControl(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4p v0, const a3real4p v1, const a3real4p vControl0, const a3real4p vControl1);
A3_INLINE a3real a3real4CalculateArcLengthHermiteTangent(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4p v0, const a3real4p v1, const a3real4p vTangent0, const a3real4p vTangent1);
A3_INLINE a3real a3real4CalculateArcLengthBezier0(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4p v0);
A3_INLINE a3real a3real4CalculateArcLengthBezier1(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4p v0, const a3real4p v1);
A3_INLINE a3real a3real4CalculateArcLengthBezier2(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4p v0, const a3real4p v1, const a3real4p v2);
A3_INLINE a3real a3real4CalculateArcLengthBezier3(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4p v0, const a3real4p v1, const a3real4p v2, const a3real4p v3);
A3_INLINE a3real a3real4CalculateArcLengthBezierN(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4 *vValues, a3index order_N);
A3_INLINE a3real a3real4CalculateArcLengthSlerp(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4p v0, const a3real4p v1);
A3_INLINE a3real a3real4CalculateArcLengthSlerpUnit(a3real4 *sampleTable_out, a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real4p v0, const a3real4p v1);

// Gram-Schmidt orthonormalization
///
A3_INLINE a3real4r a3real4GramSchmidtOrtho(a3real4p v_out, const a3real4p v, const a3real4p vBase);
A3_INLINE a3real4r a3real4GramSchmidtOrtho2(a3real4p v0_out, a3real4p v1_out, const a3real4p v0, const a3real4p v1, const a3real4p vBase);
A3_INLINE a3real4r a3real4GramSchmidtOrtho3(a3real4p v0_out, a3real4p v1_out, a3real4p v2_out, const a3real4p v0, const a3real4p v1, const a3real4p v2, const a3real4p vBase);
A3_INLINE a3real4r a3real4GramSchmidt(a3real4p v_inout, const a3real4p vBase);
A3_INLINE a3real4r a3real4GramSchmidt2(a3real4p v0_inout, a3real4p v1_inout, const a3real4p vBase);
A3_INLINE a3real4r a3real4GramSchmidt3(a3real4p v0_inout, a3real4p v1_inout, a3real4p v2_inout, const a3real4p vBase);


//-----------------------------------------------------------------------------

#ifndef A3_OPEN_SOURCE
// pre-defined constants for vector presets for quick copying
// extern if closed-source, exposed if open-source
///
extern const a3vec2 a3zeroVec2, a3oneVec2, a3xVec2, a3yVec2;
extern const a3vec3 a3zeroVec3, a3oneVec3, a3xVec3, a3yVec3, a3zVec3;
extern const a3vec4 a3zeroVec4, a3oneVec4, a3xVec4, a3yVec4, a3zVec4, a3wVec4;
#endif	// !A3_OPEN_SOURCE


//-----------------------------------------------------------------------------


A3_END_DECL


#ifdef A3_OPEN_SOURCE
#include "_inl/a3vector2_impl.inl"
#include "_inl/a3vector3_impl.inl"
#include "_inl/a3vector4_impl.inl"
#endif	// A3_OPEN_SOURCE

#endif	// !__ANIMAL3D_A3DM_VECTOR_H