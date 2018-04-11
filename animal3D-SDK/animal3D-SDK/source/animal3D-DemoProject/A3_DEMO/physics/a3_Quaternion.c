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

		a3vec3 v;
		a3real3ProductS(v.v, axis_unit, a3sind(angle_degrees / 2));
		a3real w = a3cosd(angle_degrees / 2);
		a3real4Set(q_out, v.x, v.y, v.z, w);
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

		a3vec3 temp;
		a3real3Sum(temp.v, v0_unit, v1_unit);
		a3real3Normalize(temp.v);
		a3real3Cross(q_out, v0_unit, temp.v);
		q_out[3] = a3real3Dot(v0_unit, temp.v);
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

		if (q[3] >= -1 && q[3] <= 1)
		{
			a3real3 tmp;
			a3real3Set(tmp, q[0], q[1], q[2]);
			a3real3Normalize(tmp);
			*angle_degrees_out = a3acosd(q[3]) * a3realTwo;
			axis_out = tmp;
		}
		else
		{
			*angle_degrees_out = a3realZero;
			a3real3Set(axis_out, a3realZero, a3realZero, a3realZero);
		}
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
		a3real4Set(qConj_out, -q[0], -q[1], -q[2], q[3]);
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
		a3quaternionConjugate(qInv_out, q);
		a3real magSquared = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
		a3real3DivS(qInv_out, magSquared);
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

		a3real4 v0, v1;
		a3real3ProductS(v0, qR, qL[3]);
		a3real3ProductS(v1, qL, qR[3]);
		a3real4Add(v0, v1);

		a3real3Cross(qConcat_out, qR, qL);
		a3real4Add(qConcat_out, v0);

		a3real w = qL[3] * qR[3] - a3real3Dot(qL, qR);
		qConcat_out[3] = w;
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
		a3real3 tmp, tmp2;
		a3real3Cross(tmp, q, v);
		a3real3ProductS(tmp, v, q[3]);
		a3real3Add(tmp, tmp2);

		a3real3 result;
		a3real3ProductS(result, v, a3realTwo);
		a3real3Cross(vRot_out, result, tmp);

		a3real3Add(vRot_out, v);
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
		a3real3 tmp, tmp2;
		a3real3Cross(tmp, q, v);
		a3real3ProductS(tmp, v, q[3]);
		a3real3Add(tmp, tmp2);

		a3real3 result;
		a3real3ProductS(result, v, a3realTwo);
		a3real3Cross(vRot_out, result, tmp);

		a3real3Add(vRot_out, v);
		vRot_out[3] = a3realOne;
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
		a3real dot = a3real3Dot(q1_unit, q0_unit);
		a3real4 q1;
		a3real4Set(q1, q0_unit[0], q0_unit[1], q0_unit[2], q0_unit[3]);

		if (dot < 0)
			a3real4Set(q1, -q1_unit[0], -q1_unit[1], -q1_unit[2], -q1_unit[3]);

		a3real angle = a3acosd(dot);
		a3real theta = angle * t;

		a3real s0 = a3cosd(theta) - dot * a3sind(theta) / a3sind(angle);
		a3real s1 = a3sind(theta) / a3sind(angle);

		a3real4 result1, result2;
		a3real4ProductS(result1, q0_unit, s0);
		a3real4ProductS(result2, q1_unit, s1);

		a3real3Add(result1, result2);
		qSlerp_out = result1;
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
		m_out[0][0] = q[3] * q[3] + q[0] * q[0] - q[1] * q[1] - q[2] * q[2];
		m_out[0][1] = a3realTwo * (q[0] * q[1] + q[3] * q[2]);
		m_out[0][2] = a3realTwo * (q[0] * q[2] - q[3] * q[1]);

		m_out[1][0] = a3realTwo * (q[0] * q[1] - q[3] * q[2]);
		m_out[1][1] = q[3] * q[3] - q[0] * q[0] + q[1] * q[1] - q[2] * q[2];
		m_out[1][2] = a3realTwo * (q[1] * q[2] + q[3] * q[0]);

		m_out[2][0] = a3realTwo * (q[0] * q[2] + q[3] * q[1]);
		m_out[2][1] = a3realTwo * (q[1] * q[2] - q[3] * q[0]);
		m_out[2][2] = q[3] * q[3] - q[0] * q[0] - q[1] * q[1] + q[2] * q[2];
	}
	return m_out;
}

// convert to mat4 with translation
extern inline a3real4x4r a3quaternionConvertToMat4(a3real4x4p m_out, const a3real4p q, const a3real3p translate)
{
	if (m_out && q)
	{
		a3real xx = q[0] * q[0];
		a3real xy = q[0] * q[1];
		a3real xz = q[0] * q[2];
		a3real xw = q[0] * q[3];

		a3real yy = q[1] * q[1];
		a3real yz = q[1] * q[2];
		a3real yw = q[1] * q[3];

		a3real zz = q[2] * q[2];
		a3real zw = q[2] * q[3];

		// ****TO-DO: implement
		// same as above but copy translate into fourth column
		//	and setting bottom row to (0, 0, 0, 1)
		// NOTE: matrices are COLUMN-MAJOR
		a3real m00 = 1 - 2 * (yy + zz);
		a3real m01 = 2 * (xy - zw);
		a3real m02 = 2 * (xz + yw);

		a3real m10 = 2 * (xy + zw);
		a3real m11 = 1 - 2 * (xx + zz);
		a3real m12 = 2 * (yz - xw);

		a3real m20 = 2 * (xz - yw);
		a3real m21 = 2 * (yz + xw);
		a3real m22 = 1 - 2 * (xx + yy);

		m_out[0][0] = m00;
		m_out[0][1] = m01;
		m_out[0][2] = m02;
		m_out[0][3] = a3realZero;

		m_out[1][0] = m10;
		m_out[1][1] = m11;
		m_out[1][2] = m12;
		m_out[1][3] = a3realZero;

		m_out[2][0] = m20;
		m_out[2][1] = m21;
		m_out[2][2] = m22;
		m_out[2][3] = a3realZero;

		m_out[3][0] = translate[0];
		m_out[3][1] = translate[1];
		m_out[3][2] = translate[2];
		m_out[3][3] = a3realOne;
	}

	return m_out;
}


//-----------------------------------------------------------------------------
