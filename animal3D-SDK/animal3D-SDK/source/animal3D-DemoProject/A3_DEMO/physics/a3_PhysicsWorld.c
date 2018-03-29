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
	
	a3_PhysicsWorld.c/.cpp
	Physics world function implementations.
*/


#include "a3_PhysicsWorld.h"

// include utilities on an as-needed basis
#include "animal3D/a3utility/a3_Timer.h"


// external
#include <stdio.h>
#include <string.h>


//-----------------------------------------------------------------------------

// internal utility for initializing and terminating physics world
void a3physicsInitialize_internal(a3_PhysicsWorld *world)
{
	unsigned int i, j;

	// e.g. reset all particles and/or rigid bodies
	memset(world->rigidbody, 0, sizeof(world->rigidbody));
	memset(world->particle, 0, sizeof(world->particle));
	world->t = 0.0;

	// using random rotation
	a3randomSetSeed(0);


	// set up rigid bodies to test ray picking
	world->rigidbodiesActive = 13;
	
	// ****TO-DO: 
	//	- add rotation to all
	
	// static shapes
	world->rb_sphere[0].position.y = -10.0f;
	world->rb_sphere[0].position.z = +5.0f;

	world->rb_cylinder[0].position.y = 0.0f;
	world->rb_cylinder[0].position.z = +5.0f;

	world->rb_box[0].position.y = +10.0f;
	world->rb_box[0].position.z = +5.0f;

	// moving shapes
	world->rb_sphere[1].position.x = -10.0f;
	world->rb_sphere[1].position.y = -10.0f;
	world->rb_sphere[1].position.z = +5.0f;
	world->rb_sphere[1].velocity.x = +5.0f;

	world->rb_cylinder[1].position.x = -25.0f;
	world->rb_cylinder[1].position.y = -10.0f;
	world->rb_cylinder[1].position.z = +5.0f;
	world->rb_cylinder[1].velocity.x = +5.0f;

	world->rb_box[1].position.x = -40.0f;
	world->rb_box[1].position.y = -10.0f;
	world->rb_box[1].position.z = +5.0f;
	world->rb_box[1].velocity.x = +5.0f;

	world->rb_sphere[2].position.x = -10.0f;
	world->rb_sphere[2].position.y = 0.0f;
	world->rb_sphere[2].position.z = +5.0f;
	world->rb_sphere[2].velocity.x = +5.0f;

	world->rb_cylinder[2].position.x = -25.0f;
	world->rb_cylinder[2].position.y = 0.0f;
	world->rb_cylinder[2].position.z = +5.0f;
	world->rb_cylinder[2].velocity.x = +5.0f;

	world->rb_box[2].position.x = -40.0f;
	world->rb_box[2].position.y = 0.0f;
	world->rb_box[2].position.z = +5.0f;
	world->rb_box[2].velocity.x = +5.0f;

	world->rb_sphere[3].position.x = -10.0f;
	world->rb_sphere[3].position.y = +10.0f;
	world->rb_sphere[3].position.z = +5.0f;
	world->rb_sphere[3].velocity.x = +5.0f;

	world->rb_cylinder[3].position.x = -25.0f;
	world->rb_cylinder[3].position.y = +10.0f;
	world->rb_cylinder[3].position.z = +5.0f;
	world->rb_cylinder[3].velocity.x = +5.0f;

	world->rb_box[3].position.x = -40.0f;
	world->rb_box[3].position.y = +10.0f;
	world->rb_box[3].position.z = +5.0f;
	world->rb_box[3].velocity.x = +5.0f;

	// set up hulls
	for (i = j = 0; i < 1; ++i, ++j)
	a3collisionCreateHullPlane(world->hull_ground, world->rb_ground, world->state->transform_rb + j, world->state->transformInv_rb + j,
		(a3real)(16.0), (a3real)(16.0), 1, a3axis_z);

	for (i = 0; i < 4; ++i, ++j)
		a3collisionCreateHullSphere(world->hull_sphere + i, world->rb_sphere + i, world->state->transform_rb + j, world->state->transformInv_rb + j,
			a3randomRange(a3realHalf, a3realTwo));

	for (i = 0; i < 4; ++i, ++j)
		a3collisionCreateHullCylinder(world->hull_cylinder + i, world->rb_cylinder + i, world->state->transform_rb + j, world->state->transformInv_rb + j,
			a3randomRange(a3realHalf, a3realTwo), a3randomRange(a3realOne, a3realFour), a3axis_z);

	for (i = 0; i < 4; ++i, ++j)
	{
		if (i % 2 == 0)
		{
			a3collisionCreateHullBox(world->hull_box + i, world->rb_box + i, world->state->transform_rb + j, world->state->transformInv_rb + j,
				a3randomRange(a3realOne, a3realFour), a3randomRange(a3realOne, a3realFour), a3randomRange(a3realOne, a3realFour), 2);
		}
		else
		{
			////Rotate a small amount
			a3real4Set(world->state->rotation_rb[j].v, 0, 0, -(a3real).45, 1);
			a3collisionCreateHullBox(world->hull_box + i, world->rb_box + i, world->state->transform_rb + j, world->state->transformInv_rb + j,
				a3randomRange(a3realOne, a3realFour), a3randomRange(a3realOne, a3realFour), a3randomRange(a3realOne, a3realFour), 1);

			a3real3Set(world->rigidbody[j].velocity.v, (a3real)1, 0, 0);
			/*a3collisionCreateHullBox(world->hull_box + i, world->rb_box + i, world->state->transform_rb + j, world->state->transformInv_rb + j,
				a3randomRange(a3realOne, a3realFour), a3randomRange(a3realOne, a3realFour), a3randomRange(a3realOne, a3realFour), 2);*/
		}
	}

	// no particles today
	world->particlesActive = 0;

	// raise initialized flag
	world->init = 1;

	// reset state
	a3physicsWorldStateReset(world->state);
}

void a3physicsTerminate_internal(a3_PhysicsWorld *world)
{
	// any term tasks here
}


//-----------------------------------------------------------------------------

// physics simulation
void a3physicsUpdate(a3_PhysicsWorld *world, double dt)
{
	// copy of state to edit before writing to world
	a3_PhysicsWorldState state[1] = { 0 };

	// time as real
	const a3real t_r = (a3real)(world->t);
	const a3real dt_r = (a3real)(dt);

	// generic counter
	unsigned int i;


	// ****TO-DO: 
	//	- write to state
	for (i = 0; i < world->rigidbodiesActive; ++i)
	{
		state->position_rb[i].xyz = world->rigidbody[i].position;
		state->rotation_rb[i] = world->state->rotation_rb[i];

		// rotation
		a3quaternionConvertToMat4(state->transform_rb[i].m, state->rotation_rb[i].v, state->position_rb[i].v);
		a3real4x4TransformInverseIgnoreScale(state->transformInv_rb[i].m, state->transform_rb[i].m);
	}
	state->count_rb = i;
	for (i = 0; i < world->particlesActive; ++i)
	{
		a3real4SetReal3W(state->position_p[i].v, world->particle[i].position.v, a3realOne);
	}
	state->count_p = i;


	// ****TO-DO: 
	//	- apply forces and torques


	// ****TO-DO: 
	//	- convert forces and torques
	//	- integrate to get next step using choice algorithm
	for (i = 0; i < world->rigidbodiesActive; ++i)
	{
		a3rigidbodyIntegrateEulerSemiImplicit(world->rigidbody + i, dt_r);
	}
	for (i = 0; i < world->particlesActive; ++i)
	{
		a3particleIntegrateEulerSemiImplicit(world->particle + i, dt_r);
	}


	// accumulate time
	world->t += dt;

	// write operation is locked
	if (a3physicsLockWorld(world) > 0)
	{
		// copy state to world
		*world->state = *state;
		a3physicsUnlockWorld(world);
	}
}


// physics thread
long a3physicsThread(a3_PhysicsWorld *world)
{
	// physics simulation timer
	a3_Timer physicsTimer[1] = { 0 };

	// second counter for physics (debugging)
	unsigned int currSecond = 0, prevSecond = 0;

	// create world
	a3physicsInitialize_internal(world);

	// start timer
	// rate should be set before beginning thread
	a3timerSet(physicsTimer, world->rate);
	a3timerStart(physicsTimer);

	// if lock is negative, terminate
	while (world->lock >= 0)
	{
		if (a3timerUpdate(physicsTimer))
		{
			// update timer ticked, do the thing
			a3physicsUpdate(world, physicsTimer->previousTick);

			// debug display time in seconds
			currSecond = (unsigned int)(physicsTimer->totalTime);
			if (currSecond > prevSecond)
			{
				prevSecond = currSecond;
				printf("\n physics time: %.4lf;  ticks: %llu \n     ups avg: %.4lf;  dt avg: %.4lf",
					physicsTimer->totalTime,
					physicsTimer->ticks,
					(double)physicsTimer->ticks / physicsTimer->totalTime,
					physicsTimer->totalTime / (double)physicsTimer->ticks
				);
			}
		}
	}

	// terminate world
	a3physicsTerminate_internal(world);
	return 0;
}


//-----------------------------------------------------------------------------

// reset world state
int a3physicsWorldStateReset(a3_PhysicsWorldState *worldState)
{
	unsigned int i;
	if (worldState)
	{
		//	- reset all state data appropriately
		for (i = 0; i < physicsMaxCount_rigidbody; ++i)
		{
			worldState->position_rb[i] = a3wVec4;
			//worldState->rotation_rb[i] = a3wVec4;
			worldState->transform_rb[i] = a3identityMat4;
			worldState->transformInv_rb[i] = a3identityMat4;
		}
		for (i = 0; i < physicsMaxCount_particle; ++i)
		{
			worldState->position_p[i] = a3wVec4;
		}
		return (physicsMaxCount_rigidbody + physicsMaxCount_particle);
	}
	return -1;
}


//-----------------------------------------------------------------------------

// get thread ID
#ifdef _WIN32
#include <Windows.h>
int threadID()
{
	return GetCurrentThreadId();
}
#else
#include <sys/types.h>
int threadID()
{
	return gettid();
}
#endif	// _WIN32

// mutex
extern inline int a3physicsLockWorld(a3_PhysicsWorld *world)
{
	// wait for lock to be released, then set it
	while (world->lock > 0);
	if (world->lock == 0)
	{
		world->lock = threadID();
		return world->lock;
	}
	return -1;
}

extern inline int a3physicsUnlockWorld(a3_PhysicsWorld *world)
{
	const int ret = world->lock;
	if (ret == threadID())
	{
		world->lock = 0;
		return ret;
	}
	return -1;
}


//-----------------------------------------------------------------------------
