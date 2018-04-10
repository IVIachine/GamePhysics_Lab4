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
	
	a3_PhysicsWorld.h
	Interface for physics world.
*/

/*
Tyler Chermely 0967813
Charles McGarey 0955181

EGP-425-01
Lab 4
3/30/2018

I certify that this work is
entirely our own. The assessor of this project may reproduce this project
and provide copies to other academic staff, and/or communicate a copy of
this project to a plagiarism-checking service, which may retain a copy of the
project on its database.
*/

#ifndef __ANIMAL3D_PHYSICSWORLD_H
#define __ANIMAL3D_PHYSICSWORLD_H


//-----------------------------------------------------------------------------
// physics includes

#include "a3_Collision.h"
#include "../BSP.h"

//-----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_PhysicsWorld					a3_PhysicsWorld;
	typedef struct a3_PhysicsWorldState				a3_PhysicsWorldState;
#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// counters
	enum a3_PhysicsWorldMaxCount
	{
		physicsMaxCount_particle = 32,
		physicsMaxCount_rigidbody = 32,
		physicsMaxCount_bsp = 32,
	};


//-----------------------------------------------------------------------------

	// state of a physics world: things that can be used for graphics ONLY
	//	- position and rotation... why not scale? RIGID bodies don't scale
	struct a3_PhysicsWorldState
	{
		//	- add particle position
		a3vec4 position_p[physicsMaxCount_particle];

		//	- add rigid body position and rotation
		//	- add matrix for solving inertia tensor
		//		-> we will use this and its inverse (transpose) 
		//			to resolve inertia tensors
		a3vec4 position_rb[physicsMaxCount_rigidbody];
		a3vec4 rotation_rb[physicsMaxCount_rigidbody];
		a3mat4 transform_rb[physicsMaxCount_rigidbody];
		a3mat4 transformInv_rb[physicsMaxCount_rigidbody];

		// current counts
		unsigned int count_p, count_rb;
	};


//-----------------------------------------------------------------------------

	// persistent physics world data structure
	struct a3_PhysicsWorld
	{
		//---------------------------------------------------------------------
		// SIMPLE MUTEX LOCK
		int lock;


		//---------------------------------------------------------------------
		// initialized flag
		int init;
		int framesSkipped;

		//---------------------------------------------------------------------
		// timer rate
		double rate;


		//---------------------------------------------------------------------
		// the state to store all of the things that need to go to graphics
		a3_PhysicsWorldState state[1];


		//---------------------------------------------------------------------
		// general variables pertinent to the state

		// e.g. set of rigid bodies
		union {
			a3_RigidBody rigidbody[physicsMaxCount_rigidbody];
			struct {
				a3_RigidBody rb_ground[6];
				a3_RigidBody rb_sphere[7];
			};
		};
		unsigned int rigidbodiesActive;
		
		// e.g. set of particles
		union {
			a3_Particle particle[physicsMaxCount_particle];
		};
		unsigned int particlesActive;

		// e.g. set of hulls for rigid bodies
		union {
			a3_ConvexHull hull[physicsMaxCount_rigidbody];
			struct {
				a3_ConvexHull hull_ground[6];
				a3_ConvexHull hull_sphere[7];
			};
		};

		double t;

		BSP bsps[physicsMaxCount_bsp];
		unsigned int numBSPs;
		//---------------------------------------------------------------------
	};

	
//-----------------------------------------------------------------------------

	// threaded simulation
	void a3physicsUpdate(a3_PhysicsWorld *world, double dt);
	long a3physicsThread(a3_PhysicsWorld *world);

	// world utilities
	int a3physicsWorldStateReset(a3_PhysicsWorldState *worldState);

	// mutex handling
	inline int a3physicsLockWorld(a3_PhysicsWorld *world);
	inline int a3physicsUnlockWorld(a3_PhysicsWorld *world);


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_PHYSICSWORLD_H