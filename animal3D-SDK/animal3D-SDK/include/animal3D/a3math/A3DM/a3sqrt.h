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

	a3sqrt.h
	Declarations for fast square root functions.
*/

#ifndef __ANIMAL3D_A3DM_SQRT_H
#define __ANIMAL3D_A3DM_SQRT_H


#include "animal3D/a3/a3config.h"
#include "animal3D/a3/a3macros.h"
#include "animal3D/a3/a3types_integer.h"
#include "animal3D/a3/a3types_real.h"


A3_BEGIN_DECL


//-----------------------------------------------------------------------------

// compute fast sqrt/inverse sqrt using Quake's method 
//	("0x5f3759df") using Newton's approximation method
// this only works with float precision (f32)
// compute fast sqrt using Quake's method
A3_INLINE a3f32 a3sqrtf(const a3f32 x);
A3_INLINE a3f32 a3sqrtfInverse(const a3f32 x);


// apparently the fastest ever (32-bit only, otherwise uses quake method): 
// https://www.codeproject.com/Articles/69941/Best-Square-Root-Method-Algorithm-Function-Precisi
a3f64 A3_FAST a3sqrt(const a3f64 x);
a3f64 A3_FAST a3sqrtInverse(const a3f64 x);


//-----------------------------------------------------------------------------


A3_END_DECL


#ifdef A3_OPEN_SOURCE
#include "_inl/a3sqrt_impl.inl"
#endif	// A3_OPEN_SOURCE

#endif	// !__ANIMAL3D_A3DM_SQRT_H