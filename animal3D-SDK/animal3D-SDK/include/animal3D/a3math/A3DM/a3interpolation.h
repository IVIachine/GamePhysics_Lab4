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

	a3interpolation.h
	Declarations for basic interpolation algorithms.
*/

#ifndef __ANIMAL3D_A3DM_INTERPOLATION_H
#define __ANIMAL3D_A3DM_INTERPOLATION_H


#include "animal3D/a3/a3config.h"
#include "animal3D/a3/a3macros.h"
#include "animal3D/a3/a3types_integer.h"
#include "animal3D/a3/a3types_real.h"


A3_BEGIN_DECL


//-----------------------------------------------------------------------------

// linear interpolation ("LERP"): calculate linear interpolation between 
//	reference values 'n0' and 'n1' using a interpolation parameter 'param'
//	(sometimes referred to as 't' parameter or value)
// t values between 0 and 1 yield interpolation, others yield extrapolation
// result when param = 0 is n0, and result when param = 1 is n1
A3_INLINE a3real a3lerpFunc(const a3real n0, const a3real n1, const a3real param);

// inverse linear interpolation: given a known value 'value' and reference 
//	values 'n0' and 'n1' (same as above), compute the 't' parameter 
//	that would yield 'n' as a result if used in LERP
// alternatively, the safe version makes sure n0 != n1
A3_INLINE a3real a3lerpInverse(const a3real n0, const a3real n1, const a3real value);
A3_INLINE a3real a3lerpSafeInverse(const a3real n0, const a3real n1, const a3real value);


// bi-linear interpolation: calculate lerp between each pair of numbers 
//	(n00, n01 and n10, n11) and then calculate lerp between the results
A3_INLINE a3real a3bilerp(const a3real n00, const a3real n01, const a3real n10, const a3real n11, const a3real param0, const a3real param1);

// tri-linear interpolation: extensions of the above
A3_INLINE a3real a3trilerp(const a3real n000, const a3real n001, const a3real n010, const a3real n011, const a3real n100, const a3real n101, const a3real n110, const a3real n111, const a3real param0, const a3real param1, const a3real param2);

// Catmull-Rom spline interpolation
A3_INLINE a3real a3CatmullRom(const a3real nPrev, const a3real n0, const a3real n1, const a3real nNext, const a3real param);

// two versions of cubic hermite spline (c-spline) interpolation: 
// use control values/handles (q) or tangents/derivatives (m) as inputs
A3_INLINE a3real a3HermiteControl(const a3real n0, const a3real n1, const a3real nControl0, const a3real nControl1, const a3real param);
A3_INLINE a3real a3HermiteTangent(const a3real n0, const a3real n1, const a3real nTangent0, const a3real nTangent1, const a3real param);

// Bezier spline interpolation methods (it's recursive)
// function call is Bezier<order-N>
A3_INLINE a3real a3Bezier0(const a3real n0, const a3real param);
A3_INLINE a3real a3Bezier1(const a3real n0, const a3real n1, const a3real param);
A3_INLINE a3real a3Bezier2(const a3real n0, const a3real n1, const a3real n2, const a3real param);
A3_INLINE a3real a3Bezier3(const a3real n0, const a3real n1, const a3real n2, const a3real n3, const a3real param);
A3_INLINE a3real a3BezierN(const a3real nValues[], const a3real param, a3index order_N);


// utility to calculate the arc length along a curved segment
// returns the total arc length of the segment
// outputs the samples, sample parameters and arc lengths at samples to arrays
// number of elements in each array must be 'numDivisions' + 1
// user can also choose to auto-normalize parameters
// 'paramTable_out' is an optional input; 't' params will not be stored if null
A3_INLINE a3real a3calculateArcLengthCatmullRom(a3real sampleTable_out[], a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real nPrev, const a3real n0, const a3real n1, const a3real nNext);
A3_INLINE a3real a3calculateArcLengthHermiteControl(a3real sampleTable_out[], a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real n0, const a3real n1, const a3real nControl0, const a3real nControl1);
A3_INLINE a3real a3calculateArcLengthHermiteTangent(a3real sampleTable_out[], a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real n0, const a3real n1, const a3real nTangent0, const a3real nTangent1);
A3_INLINE a3real a3calculateArcLengthBezier0(a3real sampleTable_out[], a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real n0);
A3_INLINE a3real a3calculateArcLengthBezier1(a3real sampleTable_out[], a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real n0, const a3real n1);
A3_INLINE a3real a3calculateArcLengthBezier2(a3real sampleTable_out[], a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real n0, const a3real n1, const a3real n2);
A3_INLINE a3real a3calculateArcLengthBezier3(a3real sampleTable_out[], a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real n0, const a3real n1, const a3real n2, const a3real n3);
A3_INLINE a3real a3calculateArcLengthBezierN(a3real sampleTable_out[], a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3real *nValues, a3index order_N);


// utilities to sample from a table using linear interpolation, where 
//	'valueTable' is the table of sampled values, 
//	'paramTable' is the table of sampling parameters used, 
//	'param' is the input value we are searching for, and 
//	'i' is the a3index to start searching from
// assumes all values are valid (does not check)

// find a3index and calculate value using lerp
// assumes that parameter values increase as we step through the table
A3_INLINE a3real a3sampleTableLerpIncrement(const a3real valueTable[], const a3real paramTable[], const a3real param, a3index i);

// same as the above but accounts for interpolation parameters that decrease 
//	while searching instead of increasing
A3_INLINE a3real a3sampleTableLerpDecrement(const a3real valueTable[], const a3real paramTable[], const a3real param, a3index i);

// these are the same but use a step other than 1
A3_INLINE a3real a3sampleTableLerpIncrementStep(const a3real valueTable[], const a3real paramTable[], const a3real param, const a3index step, a3index i);
A3_INLINE a3real a3sampleTableLerpDecrementStep(const a3real valueTable[], const a3real paramTable[], const a3real param, const a3index step, a3index i);

// same as the above but used to only determine at what a3index the search ends
// does not perform interpolation, but returns the a3index that should be used 
//	to grab n1 (n1 is valueTable[i], n0 is valueTable[i-1] or [i-step])
// also returns the interpolation parameter through a pointer
A3_INLINE a3index a3sampleTableLerpIncrementIndex(const a3real paramTable[], const a3real param, a3index i, a3real *param_out);
A3_INLINE a3index a3sampleTableLerpDecrementIndex(const a3real paramTable[], const a3real param, a3index i, a3real *param_out);
A3_INLINE a3index a3sampleTableLerpIncrementStepIndex(const a3real paramTable[], const a3real param, const a3index step, a3index i, a3real *param_out);
A3_INLINE a3index a3sampleTableLerpDecrementStepIndex(const a3real paramTable[], const a3real param, const a3index step, a3index i, a3real *param_out);


// generate a sample table by sampling a function that takes in one argument
// returns the total number of samples (numDivisions*numSubdivisions + 1)
// also a version that will calculate arc length while sampling
A3_INLINE a3index a3sampleTableGenerateNumSamplesRequired(const a3index numDivisions, const a3index numSubdivisions);
A3_INLINE a3index a3sampleTableGenerate(a3real sampleTable_out[], a3real paramTable_out[], const a3index numDivisions, const a3index numSubdivisions, const a3real paramMin, const a3real paramMax, const a3realfunc func);
A3_INLINE a3index a3sampleTableGenerateArcLength(a3real sampleTable_out[], a3real paramTable_out[], a3real arclenTable_out[], const a3boolean autoNormalize, const a3index numDivisions, const a3index numSubdivisions, const a3real paramMin, const a3real paramMax, const a3realfunc func);


//-----------------------------------------------------------------------------


A3_END_DECL


#ifdef A3_OPEN_SOURCE
#include "_inl/a3interpolation_impl.inl"
#endif	// A3_OPEN_SOURCE

#endif	// !__ANIMAL3D_A3DM_INTERPOLATION_H