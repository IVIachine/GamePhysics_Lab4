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
	By Daniel S. Buckstein
	
	a3_RigidBody.h
	Rigid body for physics simulation.
*/

/*
* IDs: 0955181 and 0967813
* EGP 425-01 Project 3 4/10/18
* We certify that this work is entirely our own.  The assessor of this project may reproduce this project and provide copies to other academic staff, and/or communicate a copy of this project to a plagiarism-checking service, which may retain a copy of the project on its database.
*/

#ifndef __ANIMAL3D_RIGIDBODY_H
#define __ANIMAL3D_RIGIDBODY_H


//-----------------------------------------------------------------------------

#include "a3_Particle.h"
#include "a3_Quaternion.h"


//-----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_RigidBody		a3_RigidBody;
#endif	// __cplusplus

	
//-----------------------------------------------------------------------------

	// rb
	struct a3_RigidBody
	{
		// linear dynamics
		a3vec3 position;
		a3vec3 velocity;
		a3vec3 acceleration;
		a3vec3 momentum;
		a3vec3 force;
		a3real mass, massInv;

		a3vec4 rotation;
		a3vec4 torque;
		a3vec4 velocity_a; //w(omega) dq/dt = (w)q/2
		a3vec4 acceleration_a; //alpha dw/dt = a

		a3mat3 inertiaTensor, intertiaTensorInv;
		a3mat3 inertiaTensor_t, intertiaTensorInv_t;

		a3vec3 centerMass, centerMass_t;
		// ****TO-DO: 
		//	- add angular dynamics vectors
		//	- add inertia tensor and inverse for convenience
		//	- add inertia tensor and inverse as they update
		//	- add center of mass
		//	- add center of mass as it updates

	};


//-----------------------------------------------------------------------------

	// named Euler methods

	// explicit Euler: integrate current velocity
	inline void a3rigidbodyIntegrateEulerExplicit(a3_RigidBody *rb, const a3real dt);
	
	// semi-implicit Euler: integrate next velocity
	inline void a3rigidbodyIntegrateEulerSemiImplicit(a3_RigidBody *rb, const a3real dt);
	
	// kinematic: integrate average of current and next velocities
	inline void a3rigidbodyIntegrateEulerKinematic(a3_RigidBody *rb, const a3real dt);


//-----------------------------------------------------------------------------

	// rigid body helpers

	// set mass
	inline int a3rigidbodySetMass(a3_RigidBody *rb, const a3real mass);

	// calculate inertia tensor for shapes
	//	(set mass first!)
	//	solid sphere
	inline int a3rigidbodySetLocalInertiaTensorSphereSolid(a3_RigidBody *rb, const a3real radius);
	//	hollow sphere
	inline int a3rigidbodySetLocalInertiaTensorSphereHollow(a3_RigidBody *rb, const a3real radius);
	//	solid box
	inline int a3rigidbodySetLocalInertiaTensorBoxSolid(a3_RigidBody *rb, const a3real width, const a3real height, const a3real depth);
	//	hollow box
	inline int a3rigidbodySetLocalInertiaTensorBoxHollow(a3_RigidBody *rb, const a3real width, const a3real height, const a3real depth);
	//	solid cylinder
	inline int a3rigidbodySetLocalInertiaTensorCylinderSolid(a3_RigidBody *rb, const a3real radius, const a3real height, const int bodyAxis);
	//	solid cone about apex
	inline int a3rigidbodySetLocalInertiaTensorConeSolidApex(a3_RigidBody *rb, const a3real radius, const a3real height, const int bodyAxis);
	//	rod about end
	inline int a3rigidbodySetLocalInertiaTensorRodEnd(a3_RigidBody *rb, const a3real length, const int bodyAxis);
	//	rod about center
	inline int a3rigidbodySetLocalInertiaTensorRodCenter(a3_RigidBody *rb, const a3real length, const int bodyAxis);

	// calculate total mass and center of mass from a set of influences
	inline int a3rigidbodyCalculateLocalCenterOfMass(a3_RigidBody *rb, const a3real3 *influenceList, const a3real *massList, const unsigned int count);

	// calculate inertia tensor and inverse from a set of influences
	//	(calculate center of mass first!)
	inline int a3rigidbodyCalculateLocalInertiaTensor(a3_RigidBody *rb, const a3real3 *influenceList, const a3real *massList, const unsigned int count);


	// check if particle is moving
	inline int a3rigidbodyIsMoving(const a3_RigidBody *rb);

	// check if particle is rotating
	inline int a3rigidbodyIsRotating(const a3_RigidBody *rb);

	// apply force at center of mass
	inline int a3rigidbodyApplyForceDirect(a3_RigidBody *rb, const a3real3p f);

	// apply force at other location
	inline int a3rigidbodyApplyForceLocation(a3_RigidBody *rb, const a3real3p f, const a3real3p loc);

	// convert force to acceleration
	inline int a3rigidbodyConvertForce(a3_RigidBody *rb);

	// convert torque to angular acceleration
	inline int a3rigidbodyConvertTorque(a3_RigidBody *rb);

	// reset force
	inline int a3rigidbodyResetForce(a3_RigidBody *rb);

	// reset torque
	inline int a3rigidbodyResetTorque(a3_RigidBody *rb);


	// update center of mass relative to world
	inline int a3rigidbodyUpdateCenterOfMass(a3_RigidBody *rb, const a3real4x4p transform);

	// update inertia tensor relative to world
	inline int a3rigidbodyUpdateInertiaTensor(a3_RigidBody *rb, const a3real4x4p transform);


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_RIGIDBODY_H