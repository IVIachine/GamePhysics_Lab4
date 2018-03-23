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

	a3trig.h
	Declarations for fast trigonometry functions.
*/

#ifndef __ANIMAL3D_A3DM_TRIG_H
#define __ANIMAL3D_A3DM_TRIG_H


#include "animal3D/a3/a3config.h"
#include "animal3D/a3/a3macros.h"
#include "animal3D/a3/a3types_integer.h"
#include "animal3D/a3/a3types_real.h"


A3_BEGIN_DECL


//-----------------------------------------------------------------------------

// initialize data tables for fast trigonometry
// 'samplesPerDegree' is the number of divisions in a single degree, 
//	must be 1 or greater
// 'table_out' is an externally-allocated, sufficiently-sized array to 
//	hold the function samples and params
// total number of samples in each function will be samplesPerDegree x 720 + 1
// samples for sin and cos, and the parameters will be retained, so the total 
//	number of elements in the table should be (samplesPerDegree x 720 + 1) x 3
// returns 1 if success, 0 if fail
// can also just check the number of samples that would be required
// alternatively, just set existing tables
A3_INLINE a3index a3trigInit(const a3index samplesPerDegree, a3real table_out[]);
A3_INLINE a3index a3trigInitSetTables(const a3index samplesPerDegree, const a3real table[]);
A3_INLINE a3index a3trigInitSamplesRequired(const a3index samplesPerDegree);

// delete data tables
// returns 1 if success, 0 if fail
A3_INLINE a3index a3trigFree();


//-----------------------------------------------------------------------------

// validate inputs to the functions below
// sin, cos, tan in degrees: input is [-360, 360]
// sin, cos, tan in radians: input is [-2pi, 2pi]
// asin, acos input is [-1, 1]
// atan input is (-pi/2, pi/2)
A3_INLINE a3real a3trigValid_sind(a3real x);
A3_INLINE a3real a3trigValid_sinr(a3real x);
A3_INLINE a3real a3trigValid_asin(a3real x);
A3_INLINE a3real a3trigValid_atan(a3real x);


//-----------------------------------------------------------------------------

// compute more precise trig functions by using interpolation
// sin, cos, tan assume input is [0, 360] degrees or [0, 2pi] radians, 
//	return value in [-1, 1] for sin, cos and (-inf, inf) for tan
A3_INLINE a3real a3sind(const a3real x);
A3_INLINE a3real a3cosd(const a3real x);
A3_INLINE a3real a3tand(const a3real x);
A3_INLINE a3real a3sinr(const a3real x);
A3_INLINE a3real a3cosr(const a3real x);
A3_INLINE a3real a3tanr(const a3real x);

// performs inverse-trig functions and returns value in degrees
// arcsin assumes input is [-1, 1] and returns [-90, 90] degrees or [-pi/2, pi/2] radians
// arccos assumes input is [-1, 1] and returns [0, 180] degrees or [0, pi] radians
// arctan assumes input is (-inf, inf) and returns (-90, 90) degrees or (-pi/4, pi/4) radians
A3_INLINE a3real a3asind(const a3real x);
A3_INLINE a3real a3acosd(const a3real x);
A3_INLINE a3real a3atand(const a3real x);
A3_INLINE a3real a3asinr(const a3real x);
A3_INLINE a3real a3acosr(const a3real x);
A3_INLINE a3real a3atanr(const a3real x);

// reciprocal trig functions (degrees and radians): 
// cosecant (csc) is (1 / sine)
// secant (sec) is (1 / cosine)
// cotangent (cot) is (1 / tan) or (cos / sin)
A3_INLINE a3real a3cscd(const a3real x);
A3_INLINE a3real a3secd(const a3real x);
A3_INLINE a3real a3cotd(const a3real x);
A3_INLINE a3real a3cscr(const a3real x);
A3_INLINE a3real a3secr(const a3real x);
A3_INLINE a3real a3cotr(const a3real x);

// inverse reciprocal trig functions
// same inputs and outputs as above
A3_INLINE a3real a3acscd(const a3real x);
A3_INLINE a3real a3asecd(const a3real x);
A3_INLINE a3real a3acotd(const a3real x);
A3_INLINE a3real a3acscr(const a3real x);
A3_INLINE a3real a3asecr(const a3real x);
A3_INLINE a3real a3acotr(const a3real x);

// atan2: more precise arc tangent that handles atan using two parameters
// since tan(a) = y/x, breaks down inverse solution depending on y and x
A3_INLINE a3real a3atan2d(const a3real y, const a3real x);
A3_INLINE a3real a3atan2r(const a3real y, const a3real x);

// perform long sine and cosine using a Taylor series
// this is way more precise but very slow
A3_INLINE a3index a3trigSetTaylorIterations(const a3index n);
A3_INLINE a3index a3trigGetTaylorIterations();

A3_INLINE a3real a3sinrTaylor(const a3real x);
A3_INLINE a3real a3cosrTaylor(const a3real x);
A3_INLINE a3real a3tanrTaylor(const a3real x);
A3_INLINE a3real a3sindTaylor(const a3real x);
A3_INLINE a3real a3cosdTaylor(const a3real x);
A3_INLINE a3real a3tandTaylor(const a3real x);

// calculate sin and cos at the same time
A3_INLINE void a3trigTaylor_sinr_cosr(const a3real x, a3real *sin_out, a3real *cos_out);
A3_INLINE void a3trigTaylor_sind_cosd(const a3real x, a3real *sin_out, a3real *cos_out);


//-----------------------------------------------------------------------------
// other trig-related operations

// compute the error ratio occurring from discrete geometry sampling
// e.g. since a sphere in graphics is not continuous, face centers 
//	are closer to the center than the vertices; calculate the 
//	uniform scale required to expand or contract geometry to match 
//	true continuous radii (application: volume lighting)
// total azimuth should be between [0, 360] degrees
// total elevation should be between [0, 180] degrees
A3_INLINE a3real a3pointToEdgeRatio(const a3real totalAzimuth, const a3index numSlices);
A3_INLINE a3real a3edgeToPointRatio(const a3real totalAzimuth, const a3index numSlices);
A3_INLINE a3real a3pointToFaceRatio(const a3real totalAzimuth, const a3real totalElevation, const a3index numSlices, const a3index numStacks);
A3_INLINE a3real a3faceToPointRatio(const a3real totalAzimuth, const a3real totalElevation, const a3index numSlices, const a3index numStacks);


//-----------------------------------------------------------------------------


A3_END_DECL


#ifdef A3_OPEN_SOURCE
#include "_inl/a3trig_impl.inl"
#endif	// A3_OPEN_SOURCE

#endif	// !__ANIMAL3D_A3DM_TRIG_H