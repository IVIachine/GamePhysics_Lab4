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

	a3matrix.h
	Declarations for matrix functions.
*/

#ifndef __ANIMAL3D_A3DM_MATRIX_H
#define __ANIMAL3D_A3DM_MATRIX_H


#include "animal3D/a3/a3config.h"
#include "animal3D/a3/a3macros.h"
#include "animal3D/a3/a3types_integer.h"
#include "animal3D/a3/a3types_real.h"
#include "a3vector.h"


A3_BEGIN_DECL


//-----------------------------------------------------------------------------
// NOTE: mind your handedness!
// matrices representing coordinate frames are different for right-handed 
//	and left-handed systems: 
//	
//	righty (default):						lefty (not traditional): 
//	
//		y										y
//		|										|  z
//		|										| /
//		|										|/
//		o--------x								o--------x
//	   /
//	  /
//	 z
//

#ifdef A3_MAT_LEFTY
//-----------------------------------------------------------------------------
// unsupported feature message
#pragma message ("\
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// \n\
/// ****A3DM WARNING**** User specified A3_MAT_LEFTY not supported in current version; undefining LEFTY macro.              /// \n\
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   ")
#undef A3_MAT_LEFTY
//-----------------------------------------------------------------------------
#endif	// A3_MAT_LEFTY


//-----------------------------------------------------------------------------
// NOTE: there is a difference between column-major and row-major matrices: 
//	for example, if you instantiate a 4x4 matrix with numbers 1-16 like so: 
//		a3mat4 m = { 1, 2, ..., 15, 16) };
//	then a column-major matrix (OpenGL) has the form (represented in memory): 
//
//		r0	>	|  1   5   9  13 |
//		r1	>	|  2   6  10  14 |
//		r2	>	|  3   7  11  15 |
//		r3	>	|  4   8  12  16 |
//	  rows
//				   ^   ^   ^   ^
//				  c0  c1  c2  c3
//				columns
//
//	and a row-major matrix (DirectX, or "the way it's written") has the form: 
//
//		r0	>	|  1   2   3   4 |
//		r1	>	|  5   6   7   8 |
//		r2	>	|  9  10  11  12 |
//		r3	>	| 13  14  15  16 |
//	  rows
//				   ^   ^   ^   ^
//				  c0  c1  c2  c3
//				columns
//
//	rows are always written horizontally and columns are written vertically, 
//	but as a data structure, the numbers are stored in a different order when 
//	using column- or row-major; be aware that sending a matrix that is set up 
//	as row-major into OpenGL would have the effect of sending the transpose of 
//	what was actually intended... careful!

#ifdef A3_MAT_ROWMAJOR
//-----------------------------------------------------------------------------
// unsupported feature message
#pragma message ("\
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// \n\
/// ****A3DM WARNING**** User specified A3_MAT_ROWMAJOR not supported in current version; undefining ROWMAJOR macro.        /// \n\
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   ")
#undef A3_MAT_ROWMAJOR
//-----------------------------------------------------------------------------
#endif	// A3_MAT_ROWMAJOR


//-----------------------------------------------------------------------------
// simple wrapper data structures
// like vector data structures, these serve as basic containers or wrappers 
//	for square matrix data

#ifndef __cplusplus
typedef union a3mat2 a3mat2;
typedef union a3mat3 a3mat3;
typedef union a3mat4 a3mat4;
#endif	// !__cplusplus

// 2D (2x2) matrix
union a3mat2
{
	// components
	struct { a3real x0, y0, x1, y1; };
	struct { a3real m00, m01, m10, m11; };
	struct { a3vec2 v0, v1; };
	a3vec2 v[2];

	// array (INPUT FOR OPERATIONS)
	a3real2x2 m;	// behaves as type 'a3real[][2]' or 'a3real*[2]'
	a3real2 mm;		// behaves as type 'a3real[]' or 'a3real*'
};

// 3D (3x3) matrix
union a3mat3
{
	// components
	struct { a3real x0, y0, z0  A3PAD(0,0),  x1, y1, z1  A3PAD(0,1),  x2, y2, z2  A3PAD(0,2); };
	struct { a3real m00, m01, m02  A3PAD(1,0),  m10, m11, m12  A3PAD(1,1),  m20, m21, m22  A3PAD(1,2); };
	struct { a3vec3 v0, v1, v2; };
	a3vec3 v[3];

	// array (INPUT FOR OPERATIONS)
	a3real3x3 m;	// behaves as type 'a3real[][3]' or 'a3real*[3]'
	a3real3 mm;		// behaves as type 'a3real[]' or 'a3real*'
};

// 4D (4x4) matrix
union a3mat4
{
	// components
	struct { a3real x0, y0, z0, w0, x1, y1, z1, w1, x2, y2, z2, w2, x3, y3, z3, w3; };
	struct { a3real m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33; };
	struct { a3vec4 v0, v1, v2, v3; };
	a3vec4 v[4];

	// array (INPUT FOR OPERATIONS)
	a3real4x4 m;	// behaves as type 'a3real[][4]' or 'a3real*[4]'
	a3real4 mm;		// behaves as type 'a3real[]' or 'a3real*'
};


//-----------------------------------------------------------------------------
// function declarations using array types
// for all functions that use other inputs to modify a matrix, the 
//	result/output appears as the first parameter 'm_out'
// functions that modify a matrix directly have a parameter 'm_inout'
// any matrix used as input is labelled 'm' and may have some suffix
//	to give it a better description when needed

// general function descriptions: 
// --set: fill matrix with data; setting identity yields all zeros except for 
//	ones along the diagonal (effectively the number 1 in matrix form)
//	*note: if using majors and minors, remember that we're using column-major
//	by default, so setting majors means you're providing the columns!
// --determinant: effectively the "length" of a matrix; if the determinant is 
//	exactly 1, you have a rotation matrix
// --negative/negate: negate all values in the input, store in output
// --transpose: flips all values along the diagonal; rows become columns and 
//	columns become rows
// --inverse: the inverse of a matrix M is M^-1 such that the product of a 
//	matrix and its own inverse, i.e. (M)(M^-1) or (M^-1)(M), is identity
//	*fun fact: the transpose of a rotation matrix is its inverse (fast!)
// --add/sub/mul/div: traditional math operators
// --product/concatenation: matrix multiplication; since here we are dealing 
//	with square matrices, only matrices of the same dimension are compatible 
//	to multiply (i.e. a 3x3 matrix cannot multiply with a 4x4 matrix)
// --basic transformations: generate scale or rotation matrices about the axes

//-----------------------------------------------------------------------------
// 2D MATRIX FUNCTIONS

// create 2x2 matrix
///
A3_INLINE a3real2x2r a3real2x2SetIdentity(a3real2x2p m_out);
A3_INLINE a3real2x2r a3real2x2Set(a3real2x2p m_out, const a3real x0, const a3real y0, const a3real x1, const a3real y1);
A3_INLINE a3real2x2r a3real2x2SetMajors(a3real2x2p m_out, const a3real2p v0, const a3real2p v1);
A3_INLINE a3real2x2r a3real2x2SetMinors(a3real2x2p m_out, const a3real2p v0, const a3real2p v1);
A3_INLINE a3real2x2r a3real2x2SetReal2x2(a3real2x2p m_out, const a3real2x2p m);
A3_INLINE a3real2x2r a3real2x2SetReal3x3(a3real2x2p m_out, const a3real3x3p m);
A3_INLINE a3real2x2r a3real2x2SetReal4x4(a3real2x2p m_out, const a3real4x4p m);

// unary operations
///
A3_INLINE a3real a3real2x2Determinant(const a3real2x2p m);
A3_INLINE a3real a3real2x2DeterminantInverse(const a3real2x2p m);
A3_INLINE a3real2x2r a3real2x2GetNegative(a3real2x2p m_out, const a3real2x2p m);
A3_INLINE a3real2x2r a3real2x2GetTransposed(a3real2x2p m_out, const a3real2x2p m);
A3_INLINE a3real2x2r a3real2x2GetInverse(a3real2x2p m_out, const a3real2x2p m);
A3_INLINE a3real2x2r a3real2x2Negate(a3real2x2p m_inout);
A3_INLINE a3real2x2r a3real2x2Transpose(a3real2x2p m_inout);
A3_INLINE a3real2x2r a3real2x2Invert(a3real2x2p m_inout);

// binary operations
///
A3_INLINE a3real2x2r a3real2x2Sum(a3real2x2p m_out, const a3real2x2p mL, const a3real2x2p mR);
A3_INLINE a3real2x2r a3real2x2Diff(a3real2x2p m_out, const a3real2x2p mL, const a3real2x2p mR);
A3_INLINE a3real2x2r a3real2x2ProductS(a3real2x2p m_out, const a3real2x2p m, const a3real s);
A3_INLINE a3real2x2r a3real2x2QuotientS(a3real2x2p m_out, const a3real2x2p m, const a3real s);
A3_INLINE a3real2x2r a3real2x2Add(a3real2x2p mL_inout, const a3real2x2p mR);
A3_INLINE a3real2x2r a3real2x2Sub(a3real2x2p mL_inout, const a3real2x2p mR);
A3_INLINE a3real2x2r a3real2x2MulS(a3real2x2p m_inout, const a3real s);
A3_INLINE a3real2x2r a3real2x2DivS(a3real2x2p m_inout, const a3real s);

// matrix multiplication with compatible matrix or vector
///
A3_INLINE a3real2r a3real2Real2x2ProductL(a3real2p v_out, const a3real2p v, const a3real2x2p m);
A3_INLINE a3real2r a3real2Real2x2ProductR(a3real2p v_out, const a3real2x2p m, const a3real2p v);
A3_INLINE a3real2r a3real2Real2x2MulL(a3real2p v_inout, const a3real2x2p m);
A3_INLINE a3real2r a3real2Real2x2MulR(const a3real2x2p m, a3real2p v_inout);
A3_INLINE a3real2x2r a3real2x2Product(a3real2x2p m_out, const a3real2x2p mL, const a3real2x2p mR);
A3_INLINE a3real2x2r a3real2x2ConcatL(a3real2x2p mL_inout, const a3real2x2p mR);
A3_INLINE a3real2x2r a3real2x2ConcatR(const a3real2x2p mL, a3real2x2p mR_inout);

// basic transformations: uniform scale, non-uniform scale, rotate
///
A3_INLINE a3real2x2r a3real2x2SetScale(a3real2x2p m_out, const a3real s);
A3_INLINE a3real2x2r a3real2x2SetNonUnif(a3real2x2p m_out, const a3real sx, const a3real sy);
A3_INLINE a3real2x2r a3real2x2SetRotate(a3real2x2p m_out, const a3real degrees);


//-----------------------------------------------------------------------------
// 3D MATRIX FUNCTIONS

// create 3x3 matrix
///
A3_INLINE a3real3x3r a3real3x3SetIdentity(a3real3x3p m_out);
A3_INLINE a3real3x3r a3real3x3Set(a3real3x3p m_out, const a3real x0, const a3real y0, const a3real z0, const a3real x1, const a3real y1, const a3real z1, const a3real x2, const a3real y2, const a3real z2);
A3_INLINE a3real3x3r a3real3x3SetMajors(a3real3x3p m_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real3x3r a3real3x3SetMinors(a3real3x3p m_out, const a3real3p v0, const a3real3p v1, const a3real3p v2);
A3_INLINE a3real3x3r a3real3x3SetReal2x2(a3real3x3p m_out, const a3real2x2p m);
A3_INLINE a3real3x3r a3real3x3SetReal3x3(a3real3x3p m_out, const a3real3x3p m);
A3_INLINE a3real3x3r a3real3x3SetReal4x4(a3real3x3p m_out, const a3real4x4p m);

// unary operations
///
A3_INLINE a3real a3real3x3Determinant(const a3real3x3p m);
A3_INLINE a3real a3real3x3DeterminantInverse(const a3real3x3p m);
A3_INLINE a3real3x3r a3real3x3GetNegative(a3real3x3p m_out, const a3real3x3p m);
A3_INLINE a3real3x3r a3real3x3GetTransposed(a3real3x3p m_out, const a3real3x3p m);
A3_INLINE a3real3x3r a3real3x3GetInverse(a3real3x3p m_out, const a3real3x3p m);
A3_INLINE a3real3x3r a3real3x3Negate(a3real3x3p m_inout);
A3_INLINE a3real3x3r a3real3x3Transpose(a3real3x3p m_inout);
A3_INLINE a3real3x3r a3real3x3Invert(a3real3x3p m_inout);

// binary operations
///
A3_INLINE a3real3x3r a3real3x3Sum(a3real3x3p m_out, const a3real3x3p mL, const a3real3x3p mR);
A3_INLINE a3real3x3r a3real3x3Diff(a3real3x3p m_out, const a3real3x3p mL, const a3real3x3p mR);
A3_INLINE a3real3x3r a3real3x3ProductS(a3real3x3p m_out, const a3real3x3p m, const a3real s);
A3_INLINE a3real3x3r a3real3x3QuotientS(a3real3x3p m_out, const a3real3x3p m, const a3real s);
A3_INLINE a3real3x3r a3real3x3Add(a3real3x3p mL_inout, const a3real3x3p mR);
A3_INLINE a3real3x3r a3real3x3Sub(a3real3x3p mL_inout, const a3real3x3p mR);
A3_INLINE a3real3x3r a3real3x3MulS(a3real3x3p m_inout, const a3real s);
A3_INLINE a3real3x3r a3real3x3DivS(a3real3x3p m_inout, const a3real s);

// matrix multiplication with compatible matrix or vector
///
A3_INLINE a3real3r a3real3Real3x3ProductL(a3real3p v_out, const a3real3p v, const a3real3x3p m);
A3_INLINE a3real3r a3real3Real3x3ProductR(a3real3p v_out, const a3real3x3p m, const a3real3p v);
A3_INLINE a3real3r a3real3Real3x3MulL(a3real3p v_inout, const a3real3x3p m);
A3_INLINE a3real3r a3real3Real3x3MulR(const a3real3x3p m, a3real3p v_inout);
A3_INLINE a3real3x3r a3real3x3Product(a3real3x3p m_out, const a3real3x3p mL, const a3real3x3p mR);
A3_INLINE a3real3x3r a3real3x3ConcatL(a3real3x3p mL_inout, const a3real3x3p mR);
A3_INLINE a3real3x3r a3real3x3ConcatR(const a3real3x3p mL, a3real3x3p mR_inout);

// basic transformations
// for anything more than single-axis rotations, use quaternions to handle 
//	3D rotations; quaternions are handedness-agnostic, matrices are not
///
A3_INLINE a3real3x3r a3real3x3SetScale(a3real3x3p m_out, const a3real s);
A3_INLINE a3real3x3r a3real3x3SetNonUnif(a3real3x3p m_out, const a3real sx, const a3real sy, const a3real sz);
A3_INLINE a3real3x3r a3real3x3SetRotateX(a3real3x3p m_out, const a3real degrees);
A3_INLINE a3real3x3r a3real3x3SetRotateY(a3real3x3p m_out, const a3real degrees);
A3_INLINE a3real3x3r a3real3x3SetRotateZ(a3real3x3p m_out, const a3real degrees);
A3_INLINE a3real3x3r a3real3x3SetRotateXYZ(a3real3x3p m_out, const a3real degrees_x, const a3real degrees_y, const a3real degrees_z);
A3_INLINE a3real3x3r a3real3x3SetRotateZYX(a3real3x3p m_out, const a3real degrees_x, const a3real degrees_y, const a3real degrees_z);
A3_INLINE void a3real3x3GetEulerXYZIgnoreScale(const a3real3x3p m, a3real *degrees_x_out, a3real *degrees_y_out, a3real *degrees_z_out);
A3_INLINE void a3real3x3GetEulerZYXIgnoreScale(const a3real3x3p m, a3real *degrees_x_out, a3real *degrees_y_out, a3real *degrees_z_out);


//-----------------------------------------------------------------------------
// 4D MATRIX FUNCTIONS

// create 4x4 matrix
///
A3_INLINE a3real4x4r a3real4x4SetIdentity(a3real4x4p m_out);
A3_INLINE a3real4x4r a3real4x4Set(a3real4x4p m_out, const a3real x0, const a3real y0, const a3real z0, const a3real w0, const a3real x1, const a3real y1, const a3real z1, const a3real w1, const a3real x2, const a3real y2, const a3real z2, const a3real w2, const a3real x3, const a3real y3, const a3real z3, const a3real w3);
A3_INLINE a3real4x4r a3real4x4SetMajors(a3real4x4p m_out, const a3real4p v0, const a3real4p v1, const a3real4p v2, const a3real4p v3);
A3_INLINE a3real4x4r a3real4x4SetMinors(a3real4x4p m_out, const a3real4p v0, const a3real4p v1, const a3real4p v2, const a3real4p v3);
A3_INLINE a3real4x4r a3real4x4SetReal2x2(a3real4x4p m_out, const a3real2x2p m);
A3_INLINE a3real4x4r a3real4x4SetReal3x3(a3real4x4p m_out, const a3real3x3p m);
A3_INLINE a3real4x4r a3real4x4SetReal4x4(a3real4x4p m_out, const a3real4x4p m);

// unary operations
///
A3_INLINE a3real a3real4x4Determinant(const a3real4x4p m);
A3_INLINE a3real a3real4x4DeterminantInverse(const a3real4x4p m);
A3_INLINE a3real4x4r a3real4x4GetNegative(a3real4x4p m_out, const a3real4x4p m);
A3_INLINE a3real4x4r a3real4x4GetTransposed(a3real4x4p m_out, const a3real4x4p m);
A3_INLINE a3real4x4r a3real4x4GetInverse(a3real4x4p m_out, const a3real4x4p m);
A3_INLINE a3real4x4r a3real4x4Negate(a3real4x4p m_inout);
A3_INLINE a3real4x4r a3real4x4Transpose(a3real4x4p m_inout);
A3_INLINE a3real4x4r a3real4x4Invert(a3real4x4p m_inout);

// binary operations
///
A3_INLINE a3real4x4r a3real4x4Sum(a3real4x4p m_out, const a3real4x4p mL, const a3real4x4p mR);
A3_INLINE a3real4x4r a3real4x4Diff(a3real4x4p m_out, const a3real4x4p mL, const a3real4x4p mR);
A3_INLINE a3real4x4r a3real4x4ProductS(a3real4x4p m_out, const a3real4x4p m, const a3real s);
A3_INLINE a3real4x4r a3real4x4QuotientS(a3real4x4p m_out, const a3real4x4p m, const a3real s);
A3_INLINE a3real4x4r a3real4x4Add(a3real4x4p mL_inout, const a3real4x4p mR);
A3_INLINE a3real4x4r a3real4x4Sub(a3real4x4p mL_inout, const a3real4x4p mR);
A3_INLINE a3real4x4r a3real4x4MulS(a3real4x4p m_inout, const a3real s);
A3_INLINE a3real4x4r a3real4x4DivS(a3real4x4p m_inout, const a3real s);

// matrix multiplication with compatible matrix or vector
///
A3_INLINE a3real4r a3real4Real4x4ProductL(a3real4p v_out, const a3real4p v, const a3real4x4p m);
A3_INLINE a3real4r a3real4Real4x4ProductR(a3real4p v_out, const a3real4x4p m, const a3real4p v);
A3_INLINE a3real4r a3real4Real4x4MulL(a3real4p v_inout, const a3real4x4p m);
A3_INLINE a3real4r a3real4Real4x4MulR(const a3real4x4p m, a3real4p v_inout);
A3_INLINE a3real4x4r a3real4x4Product(a3real4x4p m_out, const a3real4x4p mL, const a3real4x4p mR);
A3_INLINE a3real4x4r a3real4x4ConcatL(a3real4x4p mL_inout, const a3real4x4p mR);
A3_INLINE a3real4x4r a3real4x4ConcatR(const a3real4x4p mL, a3real4x4p mR_inout);

// basic transformations
///
A3_INLINE a3real4x4r a3real4x4SetScale(a3real4x4p m_out, const a3real s);
A3_INLINE a3real4x4r a3real4x4SetNonUnif(a3real4x4p m_out, const a3real sx, const a3real sy, const a3real sz);
A3_INLINE a3real4x4r a3real4x4SetRotateX(a3real4x4p m_out, const a3real degrees);
A3_INLINE a3real4x4r a3real4x4SetRotateY(a3real4x4p m_out, const a3real degrees);
A3_INLINE a3real4x4r a3real4x4SetRotateZ(a3real4x4p m_out, const a3real degrees);
A3_INLINE a3real4x4r a3real4x4SetRotateXYZ(a3real4x4p m_out, const a3real degrees_x, const a3real degrees_y, const a3real degrees_z);
A3_INLINE a3real4x4r a3real4x4SetRotateZYX(a3real4x4p m_out, const a3real degrees_x, const a3real degrees_y, const a3real degrees_z);
A3_INLINE void a3real4x4GetEulerXYZIgnoreScale(const a3real4x4p m, a3real *degrees_x_out, a3real *degrees_y_out, a3real *degrees_z_out);
A3_INLINE void a3real4x4GetEulerZYXIgnoreScale(const a3real4x4p m, a3real *degrees_x_out, a3real *degrees_y_out, a3real *degrees_z_out);
A3_INLINE void a3real4x4GetEulerXYZTranslateIgnoreScale(const a3real4x4p m, a3real *degrees_x_out, a3real *degrees_y_out, a3real *degrees_z_out, a3real3p translate_out);
A3_INLINE void a3real4x4GetEulerZYXTranslateIgnoreScale(const a3real4x4p m, a3real *degrees_x_out, a3real *degrees_y_out, a3real *degrees_z_out, a3real3p translate_out);


//-----------------------------------------------------------------------------
// utilities for matrices; extends matrix functionality

// concatenate/multiply/transform vectors as 4x3 (faster than full 4x4)
///
A3_INLINE a3real4r a3real4ProductTransform(a3real4p v_out, const a3real4p v, const a3real4x4p m);
A3_INLINE a3real4r a3real4TransformProduct(a3real4p v_out, const a3real4x4p m, const a3real4p v);
A3_INLINE a3real4r a3real4MulTransform(a3real4p v_inout, const a3real4x4p m);
A3_INLINE a3real4r a3real4TransformMul(const a3real4x4p m, a3real4p v_inout);
A3_INLINE a3real4x4r a3real4x4ProductTransform(a3real4x4p m_out, const a3real4x4p mL, const a3real4x4p mR);
A3_INLINE a3real4x4r a3real4x4MulTransform(a3real4x4p mL_inout, const a3real4x4p mR);
A3_INLINE a3real4x4r a3real4x4TransformMul(const a3real4x4p mL, a3real4x4p mR_inout);

// get inverse as 4x3 (faster than full 4x4)
///
A3_INLINE a3real4x4r a3real4x4TransformInverse(a3real4x4p m_out, const a3real4x4p m);
A3_INLINE a3real4x4r a3real4x4TransformInverseIgnoreScale(a3real4x4p m_out, const a3real4x4p m);
A3_INLINE a3real4x4r a3real4x4TransformInverseUniformScale(a3real4x4p m_out, const a3real4x4p m);
A3_INLINE a3real4x4r a3real4x4TransformInvert(a3real4x4p m_inout);
A3_INLINE a3real4x4r a3real4x4TransformInvertIgnoreScale(a3real4x4p m_inout);
A3_INLINE a3real4x4r a3real4x4TransformInvertUniformScale(a3real4x4p m_inout);

// prepare a "look-at" matrix given "center" and "target"
// eyePos and targetPos cannot be the same, and worldUpVec should not be zero
// optional parameter for outputting inverse
///
A3_INLINE a3real3x3r a3real3x3MakeLookAt(a3real3x3p m_out, a3real3x3p mInv_out_opt, const a3real3p eyePos, const a3real3p targetPos, const a3real3p worldUpVec);
A3_INLINE a3real4x4r a3real4x4MakeLookAt(a3real4x4p m_out, a3real4x4p mInv_out_opt, const a3real3p eyePos, const a3real3p targetPos, const a3real3p worldUpVec);

// perspective projection matrix using vertical field of view (fovyDegrees) or planes
// fovyDegrees (vertical field of view) must be between zero and 180 degrees, exclusive
// aspect represents the aspect ratio of the viewport, cannot be zero
// near must be greater than zero, and far must be greater than near
// for planes, right must be greater than left and top must be greater than bottom
// optional parameter for outputting inverse
///
A3_INLINE a3real4x4r a3real4x4MakePerspectiveProjection(a3real4x4p m_out, a3real4x4p mInv_out_opt, const a3real fovyDegrees, const a3real aspect, const a3real nearDist, const a3real farDist);
A3_INLINE a3real4x4r a3real4x4MakePerspectiveProjectionPlanes(a3real4x4p m_out, a3real4x4p mInv_out_opt, const a3real rightDist, const a3real leftDist, const a3real topDist, const a3real bottomDist, const a3real nearDist, const a3real farDist);

// orthographic projection matrix using viewport width and height
// right and left cannot be the same (width cannot be zero)
// top and bottom cannot be the same (height cannot be zero)
// near and far cannot be the same
// optional parameter for outputting inverse
///
A3_INLINE a3real4x4r a3real4x4MakeOrthographicProjection(a3real4x4p m_out, a3real4x4p mInv_out_opt, const a3real width, const a3real height, const a3real nearDist, const a3real farDist);
A3_INLINE a3real4x4r a3real4x4MakeOrthographicProjectionPlanes(a3real4x4p m_out, a3real4x4p mInv_out_opt, const a3real rightDist, const a3real leftDist, const a3real topDist, const a3real bottomDist, const a3real nearDist, const a3real farDist);

// convert to stereo projections
// interocular distance (distance between eyes) must be non-negative
// convergence distance (distance to convergence/focus) must be greater than zero (perspective)
// copies the stereo matrices to pointers passed as arguments
// optional parameter for outputting inverses
// returns pointer to right matrix
///
A3_INLINE a3real4x4r a3real4x4ConvertProjectionToStereo(a3real4x4p mL_out, a3real4x4p mR_out, a3real4x4p mInvL_out_opt, a3real4x4p mInvR_out_opt, const a3real4x4p mono, const a3real4x4p monoInv_opt, const a3real interocularDist, const a3real convergenceDist);

// check if a matrix could be perspective, orthographic, or stereo
// returns 1 if yes, 0 if no
///
A3_INLINE a3integer a3real4x4CheckPerspective(const a3real4x4p m);
A3_INLINE a3integer a3real4x4CheckOrthographic(const a3real4x4p m);
A3_INLINE a3integer a3real4x4CheckStereo(const a3real4x4p m);

// check if a matrix is perspective-stereo, or orthographic-stereo
// returns 1 if yes, 0 if not persp/ortho, -1 if persp/ortho but not stereo
///
A3_INLINE a3integer a3real4x4CheckPerspectiveStereo(const a3real4x4p m);
A3_INLINE a3integer a3real4x4CheckOrthographicStereo(const a3real4x4p m);


//-----------------------------------------------------------------------------

#ifndef A3_OPEN_SOURCE
// pre-defined constants for identity matrices for quick copying
// extern if closed-source, exposed if open-source
///
extern const a3mat2 a3identityMat2;
extern const a3mat3 a3identityMat3;
extern const a3mat4 a3identityMat4;
#endif	// !A3_OPEN_SOURCE


// define shortcut for right-handed operations
///
#define a3real2Real2x2Product	a3real2Real2x2ProductR
#define a3real3Real3x3Product	a3real3Real3x3ProductR
#define a3real4Real4x4Product	a3real4Real4x4ProductR
#define a3real2Real2x2Mul		a3real2Real2x2MulR
#define a3real3Real3x3Mul		a3real3Real3x3MulR
#define a3real4Real4x4Mul		a3real4Real4x4MulR
#define a3real2x2Concat			a3real2x2ConcatR
#define a3real3x3Concat			a3real3x3ConcatR
#define a3real4x4Concat			a3real4x4ConcatR


//-----------------------------------------------------------------------------


A3_END_DECL


#ifdef A3_OPEN_SOURCE
#include "_inl/a3matrix2_impl.inl"
#include "_inl/a3matrix3_impl.inl"
#include "_inl/a3matrix4_impl.inl"
#endif	// A3_OPEN_SOURCE

#endif	// !__ANIMAL3D_A3DM_MATRIX_H