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
	
	a3_Force.h
	Force generators and helpers.
*/

#ifndef __ANIMAL3D_FORCE_H
#define __ANIMAL3D_FORCE_H


//-----------------------------------------------------------------------------

// math library
#include "animal3D/a3math/A3DM.h"


//-----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus

#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// calculate critical damping for a particle
	inline a3real a3calcCriticalDamping(const a3real mass, const a3real springCoeff);


//-----------------------------------------------------------------------------

	// force generators

	// calculate force of gravity in freefall (using earth gravity)
	inline a3real3r a3forceGravity(a3real3p f_out, const a3real3p unitUpward, const a3real mass);

	// calculate normal force against plane
	inline a3real3r a3forceNormal(a3real3p f_out, const a3real3p f_gravity, const a3real3p unitNormal);

	// calculate sliding force on a plane
	inline a3real3r a3forceSliding(a3real3p f_out, const a3real3p f_gravity, const a3real3p f_normal);

	// calculate static friction (if particle is still)
	inline a3real3r a3forceFrictionStatic(a3real3p f_out, const a3real3p f_normal, const a3real3p f_opposing, const a3real frictionCoeffStatic);

	// calculate kinetic friction (if particle is moving)
	inline a3real3r a3forceFrictionKinetic(a3real3p f_out, const a3real3p f_normal, const a3real3p particleVelocity, const a3real frictionCoeffKinetic);

	// calculate drag force through fluid
	inline a3real3r a3forceDrag(a3real3p f_out, const a3real3p particleVelocity, const a3real3p fluidVelocity, const a3real fluidDensity, const a3real objectArea, const a3real objectDragCoeff);

	// calculate spring force using Hooke's law
	inline a3real3r a3forceSpring(a3real3p f_out, const a3real3p particlePosition, const a3real3p anchorPosition, const a3real restingLength, const a3real springCoeff);

	// calculate general damping force
	inline a3real3r a3forceDampingLinear(a3real3p f_out, const a3real3p particleVelocity, const a3real dampingCoeff);


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_FORCE_H