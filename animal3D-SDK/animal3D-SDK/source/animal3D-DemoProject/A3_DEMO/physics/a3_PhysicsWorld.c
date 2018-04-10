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

/*
* IDs: 0955181 and 0967813
* EGP 425-01 Project 3 4/10/18
* We certify that this work is entirely our own.  The assessor of this project may reproduce this project and provide copies to other academic staff, and/or communicate a copy of this project to a plagiarism-checking service, which may retain a copy of the project on its database.
*/

#include "a3_PhysicsWorld.h"

// include utilities on an as-needed basis
#include "animal3D/a3utility/a3_Timer.h"


// external
#include <stdio.h>
#include <string.h>


//-----------------------------------------------------------------------------
int setupBSP(BSP * bsp, a3real3p min, a3real3p max)
{
	a3real3Set(bsp->min.v, min[0], min[1], min[2]);
	a3real3Set(bsp->max.v, max[0], max[1], max[2]);

	//printf("Creating BSP with lower bounds (%lf %lf %lf) and upper bounds (%lf %lf %lf)\n", min[0], min[1], min[2], max[0], max[1], max[2]);
	return 0;
}


int setupBSPs(a3_PhysicsWorld * world, a3real3p min, a3real3p max, a3real3p units)
{
	a3vec3 boxUnits, tmp, tmpMin, tmpMax;

	// get the number of BSPs per dimension
	a3real3QuotientComp(boxUnits.v, a3real3Diff(tmp.v, max, min), units);

	int num = 0;
	for (int x = 0; x < (int)boxUnits.x; ++x)
	{
		for (int y = 0; y < (int)boxUnits.y; ++y)
		{
			for (int z = 0; z < (int)boxUnits.z; ++z)
			{
				++num;
				tmpMin.x = min[0] + x * units[0];
				tmpMin.y = min[1] + y * units[1];
				tmpMin.z = min[2] + z * units[2];

				a3real3Sum(tmpMax.v, tmpMin.v, units);

				setupBSP(world->bsps + num, tmpMin.v, tmpMax.v);
			}
		}
	}
	world->numBSPs = num;
	return 0;
}

int updateHulls(a3_PhysicsWorld * world)
{
	for (unsigned int j = 0; j < world->numBSPs; ++j)
	{
		world->bsps[j].numContainedHulls = 0;
	}

	for (unsigned int i = 0; i < world->rigidbodiesActive; ++i)
	{
		for (unsigned int j = 0; j < world->numBSPs; ++j)
		{
			if ((world->bsps[j].min.x <= (world->rigidbody + i)->position.x && world->bsps[j].max.x >= (world->rigidbody + i)->position.x) &&
				(world->bsps[j].min.y <= (world->rigidbody + i)->position.y && world->bsps[j].max.y >= (world->rigidbody + i)->position.y) &&
				(world->bsps[j].min.z <= (world->rigidbody + i)->position.z && world->bsps[j].max.z >= (world->rigidbody + i)->position.z))
			{
				world->bsps[j].containedHulls[world->bsps[j].numContainedHulls] = world->hull + i;
				++world->bsps[j].numContainedHulls;
				break;
			}
		}
	}
	return 0;
}



// internal utility for initializing and terminating physics world
void a3physicsInitialize_internal(a3_PhysicsWorld *world)
{
	//unsigned int i, j;

	// e.g. reset all particles and/or rigid bodies
	memset(world->rigidbody, 0, sizeof(world->rigidbody));
	memset(world->particle, 0, sizeof(world->particle));
	world->t = 0.0;

	// using random rotation
	a3randomSetSeed(0);


	// set up rigid bodies to test ray picking
	
	// ****TO-DO: 
	//	- add rotation to all
	
	// static shapes
	world->rigidbodiesActive = 0;
	
	const a3real PLANE_SIZE = 30.0f;

	world->rb_ground[0].position.x = a3realZero;
	world->rb_ground[0].position.y = a3realZero;
	world->rb_ground[0].position.z = -PLANE_SIZE;
	a3rigidbodySetMass(world->rb_ground, 1.0f);

	a3collisionCreateHullPlane(world->hull_ground + 0, world->rb_ground + 0, world->state->transform_rb + world->rigidbodiesActive, world->state->transformInv_rb + world->rigidbodiesActive,
		(a3real)(PLANE_SIZE), (a3real)(PLANE_SIZE), 2, a3axis_z);
	++world->rigidbodiesActive;

	world->rb_ground[1].position.x = a3realZero;
	world->rb_ground[1].position.y = a3realZero;
	world->rb_ground[1].position.z = PLANE_SIZE;

	a3rigidbodySetMass(world->rb_ground + 1, 0.0f);
	a3vec3 axis;
	axis.x = 0;
	axis.y = 1;
	axis.z = 0;

	a3quaternionCreateAxisAngle(world->state->rotation_rb[world->rigidbodiesActive].v, axis.v, 180.0f);
	a3collisionCreateHullPlane(world->hull_ground + world->rigidbodiesActive, world->rb_ground + world->rigidbodiesActive,
		world->state->transform_rb + world->rigidbodiesActive, world->state->transformInv_rb + world->rigidbodiesActive,
		(a3real)(PLANE_SIZE), (a3real)(PLANE_SIZE), 2, a3axis_z);
	++world->rigidbodiesActive;

	world->rb_ground[world->rigidbodiesActive].position.x = PLANE_SIZE;
	world->rb_ground[world->rigidbodiesActive].position.y = a3realZero;
	world->rb_ground[world->rigidbodiesActive].position.z = a3realZero;

	a3rigidbodySetMass(world->rb_ground + 1, 0.0f);
	axis.x = 0;
	axis.y = 1;
	axis.z = 0;

	a3quaternionCreateAxisAngle(world->state->rotation_rb[world->rigidbodiesActive].v, axis.v, 90.0f);
	a3collisionCreateHullPlane(world->hull_ground + world->rigidbodiesActive, world->rb_ground + world->rigidbodiesActive,
		world->state->transform_rb + world->rigidbodiesActive, world->state->transformInv_rb + world->rigidbodiesActive,
		(a3real)(PLANE_SIZE), (a3real)(PLANE_SIZE), 2, a3axis_z);
	++world->rigidbodiesActive;

	world->rb_ground[world->rigidbodiesActive].position.x = -PLANE_SIZE;
	world->rb_ground[world->rigidbodiesActive].position.y = a3realZero;
	world->rb_ground[world->rigidbodiesActive].position.z = a3realZero;

	a3rigidbodySetMass(world->rb_ground + 1, 0.0f);
	axis.x = 0;
	axis.y = 1;
	axis.z = 0;

	a3quaternionCreateAxisAngle(world->state->rotation_rb[world->rigidbodiesActive].v, axis.v, 270.0f);
	a3collisionCreateHullPlane(world->hull_ground + world->rigidbodiesActive, world->rb_ground + world->rigidbodiesActive,
		world->state->transform_rb + world->rigidbodiesActive, world->state->transformInv_rb + world->rigidbodiesActive,
		(a3real)(PLANE_SIZE), (a3real)(PLANE_SIZE), 2, a3axis_z);
	++world->rigidbodiesActive;

	world->rb_ground[world->rigidbodiesActive].position.x = a3realZero;
	world->rb_ground[world->rigidbodiesActive].position.y = PLANE_SIZE;
	world->rb_ground[world->rigidbodiesActive].position.z = a3realZero;

	a3rigidbodySetMass(world->rb_ground + 1, 0.0f);
	axis.x = 1;
	axis.y = 0;
	axis.z = 0;

	a3quaternionCreateAxisAngle(world->state->rotation_rb[world->rigidbodiesActive].v, axis.v, 270.0f);
	a3collisionCreateHullPlane(world->hull_ground + world->rigidbodiesActive, world->rb_ground + world->rigidbodiesActive,
		world->state->transform_rb + world->rigidbodiesActive, world->state->transformInv_rb + world->rigidbodiesActive,
		(a3real)(PLANE_SIZE), (a3real)(PLANE_SIZE), 2, a3axis_z);
	++world->rigidbodiesActive;

	world->rb_ground[world->rigidbodiesActive].position.x = a3realZero;
	world->rb_ground[world->rigidbodiesActive].position.y = -PLANE_SIZE;
	world->rb_ground[world->rigidbodiesActive].position.z = a3realZero;

	a3rigidbodySetMass(world->rb_ground + 1, 0.0f);
	axis.x = 1;
	axis.y = 0;
	axis.z = 0;

	a3quaternionCreateAxisAngle(world->state->rotation_rb[world->rigidbodiesActive].v, axis.v, 90.0f);
	a3collisionCreateHullPlane(world->hull_ground + world->rigidbodiesActive, world->rb_ground + world->rigidbodiesActive,
		world->state->transform_rb + world->rigidbodiesActive, world->state->transformInv_rb + world->rigidbodiesActive,
		(a3real)(PLANE_SIZE), (a3real)(PLANE_SIZE), 2, a3axis_z);
	++world->rigidbodiesActive;

	world->rb_sphere[0].position.y = -10.0f;
	world->rb_sphere[0].position.z = +5.0f;
	a3real3Set(world->rb_sphere[0].velocity.v, 0, 0, 0);
	a3rigidbodySetMass(world->rb_sphere, 0.5f);

	// moving shapes
	world->rb_sphere[1].position.x = -10.0f;
	world->rb_sphere[1].position.y = -10.0f;
	world->rb_sphere[1].position.z = +5.0f;
	world->rb_sphere[1].velocity.x = +15.0f;
	a3rigidbodySetMass(world->rb_sphere + 1, 1.5f);

	world->rb_sphere[2].position.x = -10.0f;
	world->rb_sphere[2].position.y = 0.0f;
	world->rb_sphere[2].position.z = +5.0f;
	world->rb_sphere[2].velocity.x = +5.0f;
	a3rigidbodySetMass(world->rb_sphere + 2, 1.5f);

	world->rb_sphere[3].position.x = -10.0f;
	world->rb_sphere[3].position.y = +10.0f;
	world->rb_sphere[3].position.z = +5.0f;
	world->rb_sphere[3].velocity.x = +5.0f;
	a3rigidbodySetMass(world->rb_sphere + 3, 0.5f);

	world->rb_sphere[4].position.x = -15.0f;
	world->rb_sphere[4].position.y = +10.0f;
	world->rb_sphere[4].position.z = +5.0f;
	world->rb_sphere[3].velocity.x = +25.0f;
	a3rigidbodySetMass(world->rb_sphere + 3, 0.5f);

	for (int i = 0; i < 5; ++i, ++world->rigidbodiesActive)
		a3collisionCreateHullSphere(world->hull_sphere + i, world->rb_sphere + i, world->state->transform_rb + world->rigidbodiesActive,
			world->state->transformInv_rb + world->rigidbodiesActive,
			a3randomRange(a3realHalf, a3realTwo));


	// no particles today
	world->particlesActive = 0;

	// raise initialized flag
	world->init = 1;
	a3vec3 min, max, units;
	a3real3Set(min.v, -100, -100, -100);
	a3real3Set(max.v, 100, 100, 100);
	a3real3Set(units.v, 100, 100, 100);

	setupBSPs(world, min.v, max.v, units.v);

	world->framesSkipped = 0;
	// reset state
	a3physicsWorldStateReset(world->state);
}

void a3physicsTerminate_internal(a3_PhysicsWorld *world)
{
	// any term tasks here
}


//-----------------------------------------------------------------------------

void a3handleCollision(a3_ConvexHullCollision* collision, a3_ConvexHull* hull_a, a3_ConvexHull* hull_b)
{
	// http://www.chrishecker.com/images/e/e7/Gdmphys3.pdf

	a3vec3 rVel;
	a3real3Diff(rVel.v, hull_a->rb->velocity.v, hull_b->rb->velocity.v);

	a3real j1 = (-a3realTwo * a3real3Dot(rVel.v, collision->normal_a[0].v)) / (a3real3Dot(collision->normal_a[0].v, 
		collision->normal_a[0].v)*(hull_a->rb->massInv + hull_b->rb->massInv));

	a3real3Add(hull_a->rb->velocity.v, a3real3ProductS(rVel.v, collision->normal_a[0].v, (j1 * hull_a->rb->massInv)));
	a3real3Add(hull_b->rb->velocity.v, a3real3ProductS(rVel.v, collision->normal_b[0].v, (j1 * hull_b->rb->massInv)));
}

// physics simulation
void a3physicsUpdate(a3_PhysicsWorld *world, double dt)
{
	// copy of state to edit before writing to world
	a3_PhysicsWorldState state[1] = { 0 };

	// time as real
	const a3real t_r = (a3real)(world->t);
	const a3real dt_r = (a3real)(dt);

	// generic counter
	unsigned int i, j;


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


	a3vec3 tmp;

	for (int i = 0; i < 5; i++)
	{
		a3forceGravity(tmp.v, a3zVec3.v, world->hull_sphere[i].rb->mass);
		a3real3Add(world->hull_sphere[i].rb->force.v, tmp.v);
		a3real3Add(world->hull_sphere[i].rb->force.v,
			a3forceDrag(tmp.v, world->hull_sphere[i].rb->velocity.v, a3zeroVec3.v, 1.2f, world->hull_sphere[i].prop[a3hullProperty_radius], .47f));
	}

	state->count_rb = i;
	for (i = 0; i < world->particlesActive; ++i)
	{
		a3real4SetReal3W(state->position_p[i].v, world->particle[i].position.v, a3realOne);
	}
	state->count_p = i;

	a3_ConvexHullCollision collision[1] = { 0 };

	if (world->framesSkipped > 5)
	{
		for (unsigned int x = 0; x < world->numBSPs; ++x)
		{
			for (i = 0; i < world->bsps[x].numContainedHulls; ++i)
			{
				for (j = 0; j < world->bsps[x].numContainedHulls; ++j)
				{
					if (i == j) continue;
					if (a3collisionTestConvexHulls(collision, world->bsps[x].containedHulls[i], world->bsps[x].containedHulls[j]) > 0)
					{
						a3handleCollision(collision, world->bsps[x].containedHulls[i], world->bsps[x].containedHulls[j]);
					}
				}
			}
		}
	}
	else
		world->framesSkipped++;

	// ****TO-DO: 
	//	- apply forces and torques

	for (i = 0; i < world->rigidbodiesActive; ++i)
	{
		a3rigidbodyIntegrateEulerKinematic(world->rigidbody + i, dt_r);
		a3real3ProductS(world->rigidbody[i].acceleration.v, world->rigidbody[i].force.v, world->rigidbody[i].massInv);
		//Add set to acceleration
		a3real4ProductS(world->rigidbody[i].acceleration_a.v, world->rigidbody[i].torque.v, world->rigidbody[i].massInv);
		a3real4Normalize(world->rigidbody[i].acceleration_a.v);

		//reset force
		a3rigidbodyResetForce(world->rigidbody + i);
	}
	for (i = 0; i < world->particlesActive; ++i)
	{
		a3particleIntegrateEulerSemiImplicit(world->particle + i, dt_r);
	}


	// accumulate time
	world->t += dt;

	updateHulls(world);


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
