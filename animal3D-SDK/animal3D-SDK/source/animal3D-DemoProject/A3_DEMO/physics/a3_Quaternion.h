/*
	Copyright 2011-2017 Daniel S. Buckstein

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
	By Daniel S. Buckstein
	
	a3_Quaternion.h
	Quaternion utilities and algorithms.
*/

#ifndef __ANIMAL3D_QUATERNION_SIMPLE_H
#define __ANIMAL3D_QUATERNION_SIMPLE_H


// math library
#include "animal3D/a3math/A3DM.h"


//-----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// create identity quaternion
	inline a3real4r a3quaternionCreateIdentity(a3real4p q_out);

	// create quaternion from normalized axis and angle
	inline a3real4r a3quaternionCreateAxisAngle(a3real4p q_out, const a3real3p axis_unit, const a3real angle_degrees);

	// create quaternion from two normalized end vectors
	inline a3real4r a3quaternionCreateDelta(a3real4p q_out, const a3real3p v0, const a3real3p v1);

	// extract axis-angle from quaternion
	inline a3real3r a3quaternionGetAxisAngle(a3real3p axis_out, a3real *angle_degrees_out, const a3real4p q);

	// conjugate
	inline a3real4r a3quaternionConjugate(a3real4p qConj_out, const a3real4p q);

	// inverse
	inline a3real4r a3quaternionInverse(a3real4p qInv_out, const a3real4p q);

	// concatenate (multiplication)
	inline a3real4r a3quaternionConcat(a3real4p qConcat_out, const a3real4p qL, const a3real4p qR);

	// rotate 3D vector
	inline a3real3r a3quaternionRotateVec3(a3real3p vRot_out, const a3real4p q, const a3real3p v);

	// rotate 4D vector/point
	inline a3real4r a3quaternionRotateVec4(a3real4p vRot_out, const a3real4p q, const a3real4p v);

	// SLERP between two unit quaternions
	inline a3real4r a3quaternionUnitSLERP(a3real4p qSlerp_out, const a3real4p q0_unit, const a3real4p q1_unit, const a3real t);

	// convert to mat3
	inline a3real3x3r a3quaternionConvertToMat3(a3real3x3p m_out, const a3real4p q);

	// convert to mat4 with translation
	inline a3real4x4r a3quaternionConvertToMat4(a3real4x4p m_out, const a3real4p q, const a3real3p translate);


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_QUATERNION_SIMPLE_H