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
	
	a3_Force.c/.cpp
	Implementation of force generators.
*/

#include "a3_Force.h"


//-----------------------------------------------------------------------------

// calculate critical damping for a particle
extern inline a3real a3calcCriticalDamping(const a3real mass, const a3real springCoeff)
{
	// a damped oscillator follows quadratic formula, 
	//	discriminant is used to determine damping: 
	// c2 - 4mk > 0		-> overdamped
	// c2 - 4mk = 0		-> critically damped
	// c2 - 4mk < 0		-> underdamped
	// e.g. http://hyperphysics.phy-astr.gsu.edu/hbase/oscda.html
	return (mass + mass)*springCoeff;
}


//-----------------------------------------------------------------------------

// force generators

// calculate force of gravity in freefall
extern inline a3real3r a3forceGravity(a3real3p f_out, const a3real3p unitUpward, const a3real mass)
{
	// standard acceleration due to gravity: 
	//	g = 9.80665 m/s2 "down"

	// ****TO-DO:
	//	- store gravity vector
	//	- convert acceleration to force: F = ma
	//		g = a -> F = mg

	return f_out;
}

// calculate normal force against plane
extern inline a3real3r a3forceNormal(a3real3p f_out, const a3real3p f_gravity, const a3real3p unitNormal)
{
	// ****TO-DO: 
	//	- projection of gravity vector onto surface vector
	//	- result would be negative, so just multiply by negative

	return f_out;
}

// calculate sliding force on a plane
extern inline a3real3r a3forceSliding(a3real3p f_out, const a3real3p f_gravity, const a3real3p f_normal)
{
	// 2D example: http://www.physicsclassroom.com/class/vectors/Lesson-3/Inclined-Planes
	// need to convert this to 3D: 
	//	- force of gravity:		Fg = mg
	//	- normal force:			Fn = Fg cos(a)
	//	- final:				F = Fg + Fn

	// ****TO-DO: 
	//	- calculate sliding force

	return f_out;
}

// calculate static friction (if object is still)
extern inline a3real3r a3forceFrictionStatic(a3real3p f_out, const a3real3p f_normal, const a3real3p f_opposing, const a3real frictionCoeffStatic)
{
	// ****TO-DO: 
	//	- calculate maximum force: k * Fn
	//	- if opposing force is less than max, apply negative opposing
	//	- if opposing force exceeds max, apply maximum in opposite direction
	//	- start with squared forces for fast compare just in case

	return f_out;
}

// calculate kinetic friction (if particle is moving)
extern inline a3real3r a3forceFrictionKinetic(a3real3p f_out, const a3real3p f_normal, const a3real3p particleVelocity, const a3real frictionCoeffKinetic)
{
	// ****TO-DO: 
	//	- calculate force: k * Fn
	//	- apply in opposite direction of travel

	return f_out;
}

// calculate drag force through fluid
extern inline a3real3r a3forceDrag(a3real3p f_out, const a3real3p particleVelocity, const a3real3p fluidVelocity, const a3real fluidDensity, const a3real objectArea, const a3real objectDragCoeff)
{
	// ****TO-DO: 
	//	- implement drag equation: F = (p u2 A C) / 2
	//		- p is the fluid density
	//		- u is the flow velocity relative to object
	//		- A is the cross-sectional area of the object
	//		- C is the object's drag coefficient

	return f_out;
}

// calculate spring force using Hooke's law
extern inline a3real3r a3forceSpring(a3real3p f_out, const a3real3p particlePosition, const a3real3p anchorPosition, const a3real restingLength, const a3real springCoeff)
{
	// ****TO-DO: 
	//	- calculate displacement (vector from anchor to position)
	//	- calculate tension using Hooke's law: F = -k(L - L0)
	//	- apply tension in direction of anchor

	return f_out;
}

// calculate general damping force
extern inline a3real3r a3forceDampingLinear(a3real3p f_out, const a3real3p particleVelocity, const a3real dampingCoeff)
{
	// ****TO-DO: 
	//	- implement simple linear damper: F = -cv

	return f_out;
}


//-----------------------------------------------------------------------------
