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
	
	a3_Ray.h
	Ray casting and picking algorithms.
*/

#ifndef __ANIMAL3D_RAY_H
#define __ANIMAL3D_RAY_H


// math library
#include "animal3D/a3math/A3DM.h"


//-----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_Ray				a3_Ray;
	typedef struct a3_RayHit			a3_RayHit;
	typedef enum a3_Axis				a3_Axis;
#endif	// __cplusplus

	
//-----------------------------------------------------------------------------

	// ray descriptor
	struct a3_Ray
	{
		a3vec4 origin;			// point of origin (w = 1)
		a3vec4 direction;		// direction vector (w = 0)
	};

	// ray hit descriptor
	struct a3_RayHit
	{
		a3vec4 hit0, hit1;		// near and far point of intersection
		a3real param0, param1;	// relative scales along ray
		int hitFlag;
	};

	// axis
	enum a3_Axis
	{
		a3axis_x,	// index 0 in a matrix
		a3axis_y,	// index 1 in a matrix
		a3axis_z,	// index 2 in a matrix
	};


//-----------------------------------------------------------------------------

	// create ray given start and end points
	inline int a3rayCreateLine(a3_Ray *ray_out, const a3vec3 origin, const a3vec3 end);

	// create ray given origin and direction (will be normalized)
	inline int a3rayCreate(a3_Ray *ray_out, const a3vec3 origin, const a3vec3 direction);

	// transform a ray
	inline int a3rayTransform(a3_Ray *ray_out, const a3_Ray *ray, const a3real4x4p transform);

	// create ray given NDC coordinate and an inverse projection matrix
	inline int a3rayCreateUnprojected(a3_Ray *ray_out, const a3vec3 v_ndc, const a3real4x4p invProj);

	// reset ray hit to ray origin
	inline int a3rayHitReset(a3_RayHit *hit_out, const a3_Ray *ray);

	// validate hit: lower flag if near distance is negative
	inline int a3rayHitValidate(a3_RayHit *hit);


//-----------------------------------------------------------------------------

	// pick against infinite plane
	inline int a3rayTestPlane(a3_RayHit *hit_out, const a3_Ray *ray, const a3_Axis normalAxis, const a3real4x4p transform);

	// pick against finite plane
	inline int a3rayTestPlaneFinite(a3_RayHit *hit_out, const a3_Ray *ray, const a3_Axis normalAxis, const a3real width, const a3real height, const a3real4x4p transform);

	// pick against disc
	inline int a3rayTestDisc(a3_RayHit *hit_out, const a3_Ray *ray, const a3_Axis normalAxis, const a3real radius, const a3real4x4p transform);

	// pick against sphere
	inline int a3rayTestSphere(a3_RayHit *hit_out, const a3_Ray *ray, const a3real radius, const a3real4x4p transform);

	// pick against infinite cylinder
	inline int a3rayTestCylinder(a3_RayHit *hit_out, const a3_Ray *ray, const a3_Axis normalAxis, const a3real radius, const a3real4x4p transform);

	// pick against finite cylinder
	inline int a3rayTestCylinderFinite(a3_RayHit *hit_out, const a3_Ray *ray, const a3_Axis normalAxis, const a3real radius, const a3real length, const a3real4x4p transform);

	// pick against axis-aligned box
	inline int a3rayTestAxisAlignedBox(a3_RayHit *hit_out, const a3_Ray *ray, const a3real width, const a3real height, const a3real depth, const a3real4x4p transform);

	// pick against bounding box
	inline int a3rayTestBoundingBox(a3_RayHit *hit_out, const a3_Ray *ray, const a3real width, const a3real height, const a3real depth, const a3real4x4p transform, const a3real4x4p transformInv);
	

//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_RAY_H