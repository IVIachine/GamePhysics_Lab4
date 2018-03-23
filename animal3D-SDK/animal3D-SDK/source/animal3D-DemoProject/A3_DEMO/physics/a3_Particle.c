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
	
	a3_Particle.c/.cpp
	Implementation of particle.
*/

#include "a3_Particle.h"


//-----------------------------------------------------------------------------

// named Euler methods

// explicit Euler: integrate current velocity
extern inline void a3particleIntegrateEulerExplicit(a3_Particle *p, const a3real dt)
{
	a3vec3 d;

	//	x(t+dt) = x(t) + f(t)dt
	//					 f(t) = dx/dt = v(t)
	//	x(t+dt) = x(t) + v(t)dt
	a3real3Add(p->position.v, a3real3ProductS(d.v, p->velocity.v, dt));

	//	v(t+dt) = v(t) + g(t)dt
	//					 g(t) = dv/dt = a(t)
	//	v(t+dt) = v(t) + a(t)dt
	a3real3Add(p->velocity.v, a3real3ProductS(d.v, p->acceleration.v, dt));
}

// semi-implicit Euler: integrate next velocity
extern inline void a3particleIntegrateEulerSemiImplicit(a3_Particle *p, const a3real dt)
{
	a3vec3 d;

	//	v(t+dt) = v(t) + a(t)dt
	a3real3Add(p->velocity.v, a3real3ProductS(d.v, p->acceleration.v, dt));

	//	x(t+dt) = x(t) + v(t+dt)dt
	a3real3Add(p->position.v, a3real3ProductS(d.v, p->velocity.v, dt));
}

// kinematic: integrate average of current and next velocities
extern inline void a3particleIntegrateEulerKinematic(a3_Particle *p, const a3real dt)
{
	a3vec3 d;

	//	x(t+dt) = x(t) + v(t)dt + a(t)dt2 / 2
	a3real3Add(p->position.v, a3real3ProductS(d.v, a3real3Sum(d.v, p->velocity.v, a3real3ProductS(d.v, p->acceleration.v, a3realHalf * dt)), dt));

	//	v(t+dt) = v(t) + a(t)dt
	a3real3Add(p->velocity.v, a3real3ProductS(d.v, p->acceleration.v, dt));
}


//-----------------------------------------------------------------------------

// particle helpers

// set mass
extern inline int a3particleSetMass(a3_Particle *p, const a3real mass)
{
	if (p)
	{
		// correctly set mass: 
		//	- valid or contingency for invalid
		if (mass > a3realZero)
		{
			p->mass = mass;
			p->massInv = a3recip(mass);
			return 1;
		}
		else
		{
			// invalid mass shall describe "static" particle
			p->mass = p->massInv = a3realZero;
			return 0;
		}
	}
	return -1;
}


// check if particle is moving
extern inline int a3particleIsMoving(const a3_Particle *p)
{
	if (p)
	{
		// determine if particle is moving (speed is not zero)
		const a3real v = a3real3LengthSquared(p->velocity.v);
		return a3isNotNearZero(v);	// safe
	//	return (v != a3realZero);	// unsafe
	}
	return -1;
}

// apply force at center of mass
extern inline int a3particleApplyForceDirect(a3_Particle *p, const a3real3p f)
{
	if (p && f)
	{
		// F(t) += F_applied
		a3real3Add(p->force.v, f);
		return 1;
	}
	return 0;
}

// convert force to acceleration
extern inline int a3particleConvertForce(a3_Particle *p)
{
	if (p)
	{
		// a(t) = F(t) / m
		a3real3ProductS(p->acceleration.v, p->force.v, p->massInv);
		return 1;
	}
	return 0;
}

// reset force
extern inline int a3particleResetForce(a3_Particle *p)
{
	if (p)
	{
		a3real3Set(p->force.v, a3realZero, a3realZero, a3realZero);
		return 1;
	}
	return -1;
}


//-----------------------------------------------------------------------------
