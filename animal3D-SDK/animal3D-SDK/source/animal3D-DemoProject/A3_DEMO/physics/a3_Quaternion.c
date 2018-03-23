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
	
	a3_Quaternion.c
	Quaternion utility implementations.
*/

#include "a3_Quaternion.h"


//-----------------------------------------------------------------------------

// create identity quaternion
extern inline a3real4r a3quaternionCreateIdentity(a3real4p q_out)
{
	if (q_out)
	{
		// xyz = 0
		// w = 1
		q_out[0] = q_out[1] = q_out[2] = a3realZero;
		q_out[3] = a3realOne;
	}
	return q_out;
}

// create quaternion from normalized axis and angle
extern inline a3real4r a3quaternionCreateAxisAngle(a3real4p q_out, const a3real3p axis_unit, const a3real angle_degrees)
{
	if (q_out && axis_unit)
	{
		// ****TO-DO: implement
		// v = sin(angle / 2) * n
		// w = cos(angle / 2)

	}
	return q_out;
}

extern inline a3real4r a3quaternionCreateDelta(a3real4p q_out, const a3real3p v0_unit, const a3real3p v1_unit)
{
	if (q_out && v0_unit && v1_unit)
	{
		// ****TO-DO: implement
		// SUPER PRO TIP for fast quaternion creation: 
		// Here are some fun facts about unit vectors: 
		//	-> a  dot  b = cos(angle)
		//	-> a cross b = sin(angle) * n
		// Since a quaternion uses half angle, we can solve by using 
		//	the unit halfway vector as 'b'!!!

	}
	return q_out;
}

// extract axis-angle from quaternion
extern inline a3real3r a3quaternionGetAxisAngle(a3real3p axis_out, a3real *angle_degrees_out, const a3real4p q)
{
	if (axis_out && angle_degrees_out && q)
	{
		// ****TO-DO: implement
		// if w is between +/- 1, 
		//	-> extract axis by normalizing vector part
		//	-> extract angle by taking inverse cosine of W and double it
		// else
		//	-> return all zeros

	}
	return axis_out;
}

// conjugate
extern inline a3real4r a3quaternionConjugate(a3real4p qConj_out, const a3real4p q)
{
	if (qConj_out && q)
	{
		// ****TO-DO: implement
		// vector part is negative

	}
	return qConj_out;
}

// inverse
extern inline a3real4r a3quaternionInverse(a3real4p qInv_out, const a3real4p q)
{
	if (qInv_out && q)
	{
		// ****TO-DO: implement
		// conjugate divided by squared magnitude

	}
	return qInv_out;
}

// concatenate (multiplication)
extern inline a3real4r a3quaternionConcat(a3real4p qConcat_out, const a3real4p qL, const a3real4p qR)
{
	if (qConcat_out && qL && qR)
	{
		// ****TO-DO: implement
		// use full formula, it's faster: 
		//	x = w0x1 + x0w1 + y0z1 - z0y1
		//	y = w0y1 - x0z1 + y0w1 + z0x1
		//	z = w0z1 + x0y1 - y0x1 + z0w1
		//	w = w0w1 - x0x1 - y0y1 - z0z1

	}
	return qConcat_out;
}

// rotate 3D vector
extern inline a3real3r a3quaternionRotateVec3(a3real3p vRot_out, const a3real4p q, const a3real3p v)
{
	if (vRot_out && q && v)
	{
		// ****TO-DO: implement
		// expand shortened formula: 
		//	v' = v + (r + r)x(r x v + wv)

	}
	return vRot_out;
}

// rotate 4D vector/point
extern inline a3real4r a3quaternionRotateVec4(a3real4p vRot_out, const a3real4p q, const a3real4p v)
{
	if (vRot_out && q && v)
	{
		// ****TO-DO: implement
		// same as above but set w component

	}
	return vRot_out;
}

// SLERP between two unit quaternions
extern inline a3real4r a3quaternionUnitSLERP(a3real4p qSlerp_out, const a3real4p q0_unit, const a3real4p q1_unit, const a3real t)
{
	if (qSlerp_out && q0_unit && q1_unit)
	{
		// ****TO-DO: implement
		// PRO TIP: if "angle" is negative, flip second quaternion
		// PRO TIP: raw SLERP formula is not enough; what if inputs are parallel?

	}
	return qSlerp_out;
}

// convert to mat3
extern inline a3real3x3r a3quaternionConvertToMat3(a3real3x3p m_out, const a3real4p q)
{
	if (m_out && q)
	{
		// ****TO-DO: implement
		// start by writing shortcuts, then apply conversion formula
		// NOTE: matrices are COLUMN-MAJOR; index like this: 
		//	m_out[column][row]
		//	e.g. upper-right would be m_out[2][0]

		// tmp
		a3real3x3SetIdentity(m_out);
	}
	return m_out;
}

// convert to mat4 with translation
extern inline a3real4x4r a3quaternionConvertToMat4(a3real4x4p m_out, const a3real4p q, const a3real3p translate)
{
	if (m_out && q)
	{
		// ****TO-DO: implement
		// same as above but copy translate into fourth column
		//	and setting bottom row to (0, 0, 0, 1)
		// NOTE: matrices are COLUMN-MAJOR

		// tmp
		a3real4x4SetIdentity(m_out);
		a3real3SetReal3(m_out[3], translate);
	}
	return m_out;
}


//-----------------------------------------------------------------------------
