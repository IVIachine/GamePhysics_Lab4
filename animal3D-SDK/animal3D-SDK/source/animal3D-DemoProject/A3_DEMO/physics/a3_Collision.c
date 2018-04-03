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

	a3_Collision.c
	Collision detection and response algorithm implementations.
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



#include "a3_Collision.h"
#include <stdio.h>

//-----------------------------------------------------------------------------

// internal utility to reset property list for hull
inline void a3collisionResetHull_internal(a3_ConvexHull *hull_out)
{
	// ****TO-DO: 
	//	- reset new things as needed
	hull_out->rb = 0;
	hull_out->transform = hull_out->transformInv = 0;

	hull_out->type = a3hullType_none;
	hull_out->flag = a3hullFlag_none;
	hull_out->axis = a3axis_z;

	// reset properties
	hull_out->prop[a3hullProperty_width] = a3realOne;
	hull_out->prop[a3hullProperty_height] = a3realOne;
	hull_out->prop[a3hullProperty_depth] = a3realOne;

	hull_out->prop[a3hullProperty_halfwidth] = a3realHalf;
	hull_out->prop[a3hullProperty_halfheight] = a3realHalf;
	hull_out->prop[a3hullProperty_halfdepth] = a3realHalf;

	hull_out->prop[a3hullProperty_halfwidthSq] = a3realQuarter;
	hull_out->prop[a3hullProperty_halfheightSq] = a3realQuarter;
	hull_out->prop[a3hullProperty_halfdepthSq] = a3realQuarter;
}


//-----------------------------------------------------------------------------

// internal collision tests
//	- there are many more, but these will do for now

inline int a3collisionTestPointSphere(
	const a3real3p pointPosition, const a3real3p sphereCenter, const a3real sphereRadiusSq, a3real3p diff_tmp)
{
	// ****TO-DO: 
	//	- implement test
	a3real3Diff(diff_tmp, pointPosition, sphereCenter);
	return (a3real3LengthSquared(diff_tmp) <= sphereRadiusSq);
}

inline int a3collisionTestPointAABB(
	const a3real3p pointPosition_localToAABB, const a3real3p aabbMinExtents, const a3real3p aabbMaxExtents)
{
	// min_x <= x <= max_x
	// min_y <= y <= max_y
	// min_z <= z <= max_z
	// COLLISION!
	if ((aabbMinExtents[0] <= pointPosition_localToAABB[0] && aabbMaxExtents[0] >= pointPosition_localToAABB[0]) &&
		(aabbMinExtents[1] <= pointPosition_localToAABB[1] && aabbMaxExtents[1] >= pointPosition_localToAABB[1]) &&
		(aabbMinExtents[2] <= pointPosition_localToAABB[2] && aabbMaxExtents[2] >= pointPosition_localToAABB[2]))
	{
		return 1;
	}

	return 0;
}

inline int a3collisionTestPlaneSphere(
	const a3real3p planeCenter, const a3real3p planeTangent, const a3real3p planeBitangent, const a3real3p planeNormal, const a3real planeHalfWidth, const a3real planeHalfHeight, const a3real3p sphereCenter, const a3real sphereRadius, a3real3p diff_tmp)
{
	a3real3Projected(diff_tmp, sphereCenter, planeNormal);
	a3real3Sub(diff_tmp, planeCenter);
	if (a3real3LengthSquared(diff_tmp) > sphereRadius * sphereRadius) return 0;

	a3real3Projected(diff_tmp, sphereCenter, planeTangent);
	a3real3Sub(diff_tmp, planeCenter);
	if (a3real3LengthSquared(diff_tmp) > (sphereRadius + planeHalfWidth) * (sphereRadius + planeHalfWidth)) return 0;

	a3real3Projected(diff_tmp, sphereCenter, planeBitangent);
	a3real3Sub(diff_tmp, planeCenter);
	if (a3real3LengthSquared(diff_tmp) > (sphereRadius + planeHalfHeight) * (sphereRadius + planeHalfHeight)) return 0;

	return 1;
}

inline int a3collisionTestPlaneAABB(
	const a3real3p planeCenter_localToAABB, const a3real3p planeTangent_localToAABB, const a3real3p planeBitangent_localToAABB, const a3real3p planeNormal_localToAABB, const a3real planeHalfWidth, const a3real planeHalfHeight, const a3real3p aabbMinExtents, const a3real3p aabbMaxExtents, a3real3p diff_tmp)
{
	// ****TO-DO: 
	//	- implement test
	// treat plane as AABB in 2 dimensions
	// this is not the best way to do it, but it is A way
	a3real3 planeMin, planeMax, tmp;
	a3real3Add(a3real3Sum(planeMin, a3real3ProductS(planeMin, planeTangent_localToAABB, -planeHalfWidth), a3real3ProductS(tmp, planeBitangent_localToAABB, -planeHalfHeight)), planeCenter_localToAABB);
	a3real3Add(a3real3Sum(planeMax, a3real3ProductS(planeMax, planeTangent_localToAABB, planeHalfWidth), a3real3ProductS(tmp, planeBitangent_localToAABB, planeHalfHeight)), planeCenter_localToAABB);

	if ((aabbMinExtents[0] <= planeMax[0] && aabbMaxExtents[0] >= planeMin[0]) &&
		(aabbMinExtents[1] <= planeMax[1] && aabbMaxExtents[1] >= planeMin[1]) &&
		(aabbMinExtents[2] <= planeMax[2] && aabbMaxExtents[2] >= planeMin[2]))
	{
		return 1;
	}

	return 0;
}

inline int a3collisionTestSpheres(a3_ConvexHullCollision *collision_out,
	const a3real3p sphereCenter_a, const a3real sphereRadius_a, const a3real3p sphereCenter_b, const a3real sphereRadius_b, a3real3p diff_tmp)
{
	// ****TO-DO: 
	//	- implement test
	const a3real sumRadii = sphereRadius_a + sphereRadius_b;
	a3real3Diff(diff_tmp, sphereCenter_a, sphereCenter_b);
	if (a3real3LengthSquared(diff_tmp) <= sumRadii * sumRadii)
	{
		// generate contacts
		a3real3Normalize(diff_tmp);
		a3real3Set(collision_out->normal_b[0].v, diff_tmp[0], diff_tmp[1], diff_tmp[2]);
		
		a3real3ProductS(collision_out->contact_b[0].v, diff_tmp, sphereRadius_b);
		a3real3Add(collision_out->contact_b[0].v, sphereCenter_b);
		
		a3real3Negate(diff_tmp);
		a3real3Set(collision_out->normal_a[0].v, diff_tmp[0], diff_tmp[1], diff_tmp[2]);
		a3real3ProductS(collision_out->contact_a[0].v, diff_tmp, sphereRadius_a);
		a3real3Add(collision_out->contact_a[0].v, sphereCenter_a);
		
		collision_out->contactCount_a = collision_out->contactCount_b = 1;
		return 1;
	}
	return 0;
}

inline int a3collisionTestSphereAABB(a3_ConvexHullCollision *collision_out,
	const a3real3p sphereCenter_localToAABB, const a3real sphereRadius, const a3real3p aabbMinExtents, const a3real3p aabbMaxExtents, a3real3p diff_tmp)
{
	a3vec3 tmp;
	// closest point to the sphere center on the box
	tmp.x = (sphereCenter_localToAABB[0] < aabbMinExtents[0]) ? aabbMinExtents[0] : (sphereCenter_localToAABB[0] > aabbMaxExtents[0]) ? aabbMaxExtents[0] : sphereCenter_localToAABB[0];
	tmp.y = (sphereCenter_localToAABB[1] < aabbMinExtents[1]) ? aabbMinExtents[1] : (sphereCenter_localToAABB[1] > aabbMaxExtents[1]) ? aabbMaxExtents[1] : sphereCenter_localToAABB[1];
	tmp.z = (sphereCenter_localToAABB[2] < aabbMinExtents[2]) ? aabbMinExtents[2] : (sphereCenter_localToAABB[2] > aabbMaxExtents[2]) ? aabbMaxExtents[2] : sphereCenter_localToAABB[2];

	a3real3Diff(diff_tmp, tmp.v, sphereCenter_localToAABB);
	if (a3real3LengthSquared(diff_tmp) <= sphereRadius * sphereRadius)
	{
		// sphere is a
		// box is b
		collision_out->contactCount_a = collision_out->contactCount_b = 1;
		collision_out->contact_a[0] = collision_out->contact_b[0] = tmp;
		a3real3Normalize(diff_tmp);
		a3real3Set(collision_out->normal_b[0].v, diff_tmp[0], diff_tmp[1], diff_tmp[2]);
		a3real3Set(collision_out->normal_a[0].v, -diff_tmp[0], -diff_tmp[1], -diff_tmp[2]);
		return 1;
	}


	return 0;
}

inline int a3collisionTestAABBs(a3_ConvexHullCollision *collision_out,
	const a3real3p aabbMinExtents_a, const a3real3p aabbMaxExtents_a, const a3real3p aabbMinExtents_b, const a3real3p aabbMaxExtents_b, a3real3p diff_tmp)
{
	// find the points of the intersection
	
	if ((aabbMinExtents_a[0] <= aabbMaxExtents_b[0] && aabbMaxExtents_a[0] >= aabbMinExtents_b[0]) &&
		(aabbMinExtents_a[1] <= aabbMaxExtents_b[1] && aabbMaxExtents_a[1] >= aabbMinExtents_b[1]) &&
		(aabbMinExtents_a[2] <= aabbMaxExtents_b[2] && aabbMaxExtents_a[2] >= aabbMinExtents_b[2]))
	{
		// min is collision a 5
		// max is collision b 5
		collision_out->contact_a[5].x = aabbMinExtents_a[0] < aabbMinExtents_b[0] ? aabbMinExtents_b[0] : aabbMinExtents_b[0];
		collision_out->contact_a[5].y = aabbMinExtents_a[1] < aabbMinExtents_b[1] ? aabbMinExtents_b[1] : aabbMinExtents_b[1];
		collision_out->contact_a[5].z = aabbMinExtents_a[2] < aabbMinExtents_b[2] ? aabbMinExtents_b[2] : aabbMinExtents_b[2];


		collision_out->contact_b[5].x = aabbMaxExtents_a[0] < aabbMaxExtents_b[0] ? aabbMaxExtents_b[0] : aabbMaxExtents_b[0];
		collision_out->contact_b[5].y = aabbMaxExtents_a[1] < aabbMaxExtents_b[1] ? aabbMaxExtents_b[1] : aabbMaxExtents_b[1];
		collision_out->contact_b[5].z = aabbMaxExtents_a[2] < aabbMaxExtents_b[2] ? aabbMaxExtents_b[2] : aabbMaxExtents_b[2];

		diff_tmp[0] = (collision_out->contact_a[5].x + collision_out->contact_b[5].x) / 2;
		diff_tmp[1] = (collision_out->contact_a[5].y + collision_out->contact_b[5].y) / 2;
		diff_tmp[2] = (collision_out->contact_a[5].z + collision_out->contact_b[5].z) / 2;
		// contact point



		return 1;
	}

	return 0;
}


//-----------------------------------------------------------------------------

// create point hull
extern inline int a3collisionCreateHullPoint(a3_ConvexHull *hull_out, const a3_RigidBody *rb)
{
	if (hull_out && rb)
	{
		// ****TO-DO: 
		//	- set properties
		a3collisionResetHull_internal(hull_out);

		return hull_out->type;
	}
	return -1;
}

// create plane hull
extern inline int a3collisionCreateHullPlane(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const a3real width, const a3real height, const int isAxisAligned, const a3_Axis normalAxis)
{
	if (hull_out && rb)
	{
		// ****TO-DO: 
		//	- set properties
		a3collisionResetHull_internal(hull_out);
		hull_out->rb = rb;
		hull_out->transform = transform;
		hull_out->transformInv = transformInv;

		hull_out->type = a3hullType_plane;
		hull_out->prop[a3hullProperty_width] = width;
		hull_out->prop[a3hullProperty_height] = height;

		hull_out->prop[a3hullProperty_halfwidth] = width * a3realHalf;
		hull_out->prop[a3hullProperty_halfheight] = height * a3realHalf;

		hull_out->axis = normalAxis;
		hull_out->prop[a3hullFlag_isAxisAligned] = (a3real)isAxisAligned;

		return hull_out->type;
	}
	return -1;
}

// create box hull
extern inline int a3collisionCreateHullBox(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const a3real width, const a3real height, const a3real depth, const int isAxisAligned)
{
	if (hull_out && rb)
	{
		// ****TO-DO: 
		//	- set properties
		a3collisionResetHull_internal(hull_out);
		hull_out->rb = rb;
		hull_out->transform = transform;
		hull_out->transformInv = transformInv;

		hull_out->type = a3hullType_box;
		hull_out->prop[a3hullProperty_width] = width;
		hull_out->prop[a3hullProperty_height] = height;
		hull_out->prop[a3hullProperty_depth] = depth;

		hull_out->prop[a3hullProperty_halfwidth] = width * a3realHalf;
		hull_out->prop[a3hullProperty_halfheight] = height * a3realHalf;
		hull_out->prop[a3hullProperty_halfdepth] = depth * a3realHalf;

		hull_out->prop[a3hullFlag_isAxisAligned] = (a3real)isAxisAligned;

		return hull_out->type;
	}
	return -1;
}

// create sphere hull
extern inline int a3collisionCreateHullSphere(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const a3real radius)
{
	if (hull_out && rb && radius > a3realZero)
	{
		// ****TO-DO: 
		//	- set properties
		a3collisionResetHull_internal(hull_out);
		hull_out->rb = rb;
		hull_out->transform = transform;
		hull_out->transformInv = transformInv;

		hull_out->type = a3hullType_sphere;
		hull_out->prop[a3hullProperty_radius] = radius;
		hull_out->prop[a3hullProperty_radiusSq] = radius * radius;
		return hull_out->type;
	}
	return -1;
}

// create cylinder hull
extern inline int a3collisionCreateHullCylinder(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const a3real radius, const a3real length, const a3_Axis normalAxis)
{
	if (hull_out && rb)
	{
		// ****TO-DO: 
		//	- set properties
		a3collisionResetHull_internal(hull_out);
		hull_out->rb = rb;
		hull_out->transform = transform;
		hull_out->transformInv = transformInv;

		hull_out->axis = normalAxis;

		hull_out->type = a3hullType_cylinder;
		hull_out->prop[a3hullProperty_radius] = radius;
		hull_out->prop[a3hullProperty_radiusSq] = radius * radius;
		hull_out->prop[a3hullProperty_length] = length;
		return hull_out->type;
	}
	return -1;
}

// create mesh hull
extern inline int a3collisionCreateHullMesh(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const void *points, const unsigned int pointCount, const int is3D)
{
	if (hull_out && rb)
	{
		// ****TO-DO: 
		//	- set properties
		a3collisionResetHull_internal(hull_out);

		return hull_out->type;
	}
	return -1;
}


//-----------------------------------------------------------------------------

//Rotate vertex around axis
a3vec3 rotatePosition(a3vec3 n, a3real angle, a3vec3 p)
{
	a3mat3 m;
	a3real3x3Set(m.m,
		n.x*n.x * (1.0f - a3cosd(angle)) + a3cosd(angle),
		n.x*n.y * (1.0f - a3cosd(angle)) + n.z * a3sind(angle),
		n.x*n.z * (1.0f - a3cosd(angle)) - n.y * a3sind(angle),
		n.y*n.x * (1.0f - a3cosd(angle)) - n.z * a3sind(angle),
		n.y*n.y * (1.0f - a3cosd(angle)) + a3cosd(angle),
		n.y*n.z * (1.0f - a3cosd(angle)) + n.x * a3sind(angle),
		n.z*n.x * (1.0f - a3cosd(angle)) + n.y * a3sind(angle),
		n.z*n.y * (1.0f - a3cosd(angle)) - n.x * a3sind(angle),
		n.z*n.z * (1.0f - a3cosd(angle)) + a3cosd(angle)
	);

	a3real3Real3x3MulR(m.m, p.v);
	return p;
}

// high-level collision test
extern inline int a3collisionTestConvexHulls(a3_ConvexHullCollision *collision_out, const a3_ConvexHull *hull_a, const a3_ConvexHull *hull_b)
{
	if (collision_out && hull_a && hull_b)
	{
		//Perform collision tests
		//****NOTE:
		//Charlie did AABB/OBB, Plane/OBB
		//Tyler did rest with help from Charlie

		int status = 0;
		a3real3 tmp;

		switch (hull_a->type)
		{
		case a3hullType_sphere:
		{
			switch (hull_b->type)
			{
			case a3hullType_sphere:
			{
				collision_out->hull_a = hull_a;
				collision_out->hull_b = hull_b;
				status = a3collisionTestSpheres(collision_out, hull_a->transform->v3.v, hull_a->prop[a3hullProperty_radius],
					hull_b->transform->v3.v, hull_b->prop[a3hullProperty_radius], tmp);
			}
			break;
			case a3hullType_box:
			{
				//box sphere
				a3vec3 minB, maxB, diff;
				collision_out->hull_a = hull_a;
				collision_out->hull_b = hull_b;
				if (hull_b->prop[a3hullFlag_isAxisAligned] == 2)
				{
					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];

					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];

					minB.z = hull_b->transform->v3.z - hull_b->prop[a3hullProperty_halfdepth];
					maxB.z = hull_b->transform->v3.z + hull_b->prop[a3hullProperty_halfdepth];
					
					status = a3collisionTestSphereAABB(collision_out, hull_a->transform->v3.v, hull_a->prop[a3hullProperty_radius], minB.v, maxB.v, diff.v);
				}
				else
				{
					a3vec3 temp;
					a3real3Diff(temp.v, hull_a->transform->v3.v, hull_b->transform->v3.v);
					a3real4ProductTransform(temp.v, temp.v, hull_a->transformInv->m);
					a3real3Add(temp.v, hull_b->transform->v3.v);

					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];

					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];

					minB.z = hull_b->transform->v3.z - hull_b->prop[a3hullProperty_halfdepth];
					maxB.z = hull_b->transform->v3.z + hull_b->prop[a3hullProperty_halfdepth];

					status = a3collisionTestSphereAABB(collision_out, temp.v, hull_a->prop[a3hullProperty_radius], minB.v, maxB.v, diff.v);
				}
			}
			break;
			case a3hullType_plane:
			{
				//sphere plane
				a3vec3 minB, maxB, diff;
				if (hull_b->prop[a3hullFlag_isAxisAligned] == 2)
				{
					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];

					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];

					minB.z = hull_b->transform->v3.z;
					maxB.z = hull_b->transform->v3.z;
					status = a3collisionTestSphereAABB(collision_out, hull_a->transform->v3.v, hull_a->prop[a3hullProperty_radius], minB.v, maxB.v, diff.v);
				}
				else
				{
					a3vec3 temp;
					a3real3Diff(temp.v, hull_a->transform->v3.v, hull_b->transform->v3.v);
					a3real4ProductTransform(temp.v, temp.v, hull_a->transformInv->m);
					a3real3Add(temp.v, hull_b->transform->v3.v);

					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];

					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];

					minB.z = hull_b->transform->v3.z;
					maxB.z = hull_b->transform->v3.z;

					status = a3collisionTestSphereAABB(collision_out, temp.v, hull_a->prop[a3hullProperty_radius], minB.v, maxB.v, diff.v);
				}
			}
			break;
			default:
				break;
			}
		}
		break;
		case a3hullType_box:
		{
			switch (hull_b->type)
			{
			case a3hullType_sphere:
			{
				//box sphere
				a3vec3 minA, maxA, diff;
				collision_out->hull_a = hull_b;
				collision_out->hull_b = hull_a;
				if (hull_b->prop[a3hullFlag_isAxisAligned] == 2)
				{
					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_halfwidth];

					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_halfheight];

					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_halfdepth];
					status = a3collisionTestSphereAABB(collision_out, hull_b->transform->v3.v, hull_b->prop[a3hullProperty_radius], minA.v, maxA.v, diff.v);
				}
				else
				{
					a3vec3 temp;
					a3real3Diff(temp.v, hull_b->transform->v3.v, hull_a->transform->v3.v);
					a3real4ProductTransform(temp.v, temp.v, hull_b->transform->m);
					a3real3Add(temp.v, hull_a->transform->v3.v);

					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_halfdepth];

					status = a3collisionTestSphereAABB(collision_out, temp.v, hull_b->prop[a3hullProperty_radius], minA.v, maxA.v, diff.v);
				}
			}
			break;
			case a3hullType_box:
			{
				//box box
				a3vec3 minA, maxA, minB, maxB;
				if (hull_a->prop[a3hullFlag_isAxisAligned] == 2 && hull_b->prop[a3hullFlag_isAxisAligned] == 2)
				{
					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];			
					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = hull_b->transform->v3.z - hull_b->prop[a3hullProperty_halfdepth];
					maxB.z = hull_b->transform->v3.z + hull_b->prop[a3hullProperty_halfdepth];

					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
				}
				else if (hull_a->prop[a3hullFlag_isAxisAligned] == 2 && hull_b->prop[a3hullFlag_isAxisAligned] != 2)
				{
					// calc new a transform
					a3vec3 newTransA;
					a3real3Diff(newTransA.v, hull_b->transform->v3.v, hull_a->transform->v3.v);
					a3real4ProductTransform(newTransA.v, newTransA.v, hull_b->transformInv->m);
					a3real3Add(newTransA.v, hull_a->transform->v3.v);

					minA.x = newTransA.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = newTransA.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = newTransA.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = newTransA.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = newTransA.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = newTransA.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];
					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = hull_b->transform->v3.z - hull_b->prop[a3hullProperty_halfdepth];
					maxB.z = hull_b->transform->v3.z + hull_b->prop[a3hullProperty_halfdepth];

					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
				}
				else if (hull_a->prop[a3hullFlag_isAxisAligned] != 2 && hull_b->prop[a3hullFlag_isAxisAligned] == 2)
				{
					a3vec3 newTransB;
					a3real3Diff(newTransB.v, hull_a->transform->v3.v, hull_b->transform->v3.v);
					a3real4ProductTransform(newTransB.v, newTransB.v, hull_a->transformInv->m);
					a3real3Add(newTransB.v, hull_b->transform->v3.v);

					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = newTransB.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = newTransB.x + hull_b->prop[a3hullProperty_halfwidth];
					minB.y = newTransB.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = newTransB.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = newTransB.z - hull_b->prop[a3hullProperty_halfdepth];
					maxB.z = newTransB.z + hull_b->prop[a3hullProperty_halfdepth];

					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
				}
				else
				{
					// calc new a transform
					a3vec3 newTransA;
					a3real3Diff(newTransA.v, hull_b->transform->v3.v, hull_a->transform->v3.v);
					a3real4ProductTransform(newTransA.v, newTransA.v, hull_b->transformInv->m);
					a3real3Add(newTransA.v, hull_a->transform->v3.v);

					minA.x = newTransA.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = newTransA.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = newTransA.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = newTransA.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = newTransA.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = newTransA.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];
					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = hull_b->transform->v3.z - hull_b->prop[a3hullProperty_halfdepth];
					maxB.z = hull_b->transform->v3.z + hull_b->prop[a3hullProperty_halfdepth];

					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
					if (status == 0) break;

					a3vec3 newTransB;
					a3real3Diff(newTransB.v, hull_a->transform->v3.v, hull_b->transform->v3.v);
					a3real4ProductTransform(newTransB.v, newTransB.v, hull_a->transformInv->m);
					a3real3Add(newTransB.v, hull_b->transform->v3.v);

					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = newTransB.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = newTransB.x + hull_b->prop[a3hullProperty_halfwidth];
					minB.y = newTransB.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = newTransB.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = newTransB.z - hull_b->prop[a3hullProperty_halfdepth];
					maxB.z = newTransB.z + hull_b->prop[a3hullProperty_halfdepth];

					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
				}
			}
			break;
			case a3hullType_plane:
			{
				//box plane
				a3vec3 minA, maxA, minB, maxB;
				if (hull_a->prop[a3hullFlag_isAxisAligned] == 2 && hull_b->prop[a3hullFlag_isAxisAligned] == 2)
				{
					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];
					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = hull_b->transform->v3.z;
					maxB.z = hull_b->transform->v3.z;
					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
				}
				else if (hull_a->prop[a3hullFlag_isAxisAligned] == 2 && hull_b->prop[a3hullFlag_isAxisAligned] != 2)
				{
					a3vec3 newTransA;
					a3real3Diff(newTransA.v, hull_b->transform->v3.v, hull_a->transform->v3.v);
					a3real4ProductTransform(newTransA.v, newTransA.v, hull_b->transformInv->m);
					a3real3Add(newTransA.v, hull_a->transform->v3.v);

					minA.x = newTransA.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = newTransA.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = newTransA.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = newTransA.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = newTransA.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = newTransA.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];
					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = hull_b->transform->v3.z;
					maxB.z = hull_b->transform->v3.z;

					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
				}
				else if (hull_a->prop[a3hullFlag_isAxisAligned] != 2 && hull_b->prop[a3hullFlag_isAxisAligned] == 2)
				{
					a3vec3 newTransB;
					a3real3Diff(newTransB.v, hull_a->transform->v3.v, hull_b->transform->v3.v);
					a3real4ProductTransform(newTransB.v, newTransB.v, hull_a->transformInv->m);
					a3real3Add(newTransB.v, hull_b->transform->v3.v);

					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = newTransB.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = newTransB.x + hull_b->prop[a3hullProperty_halfwidth];
					minB.y = newTransB.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = newTransB.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = hull_b->transform->v3.z;
					maxB.z = hull_b->transform->v3.z;

					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
				}
				else
				{
					// calc new a transform
					a3vec3 newTransA;
					a3real3Diff(newTransA.v, hull_b->transform->v3.v, hull_a->transform->v3.v);
					a3real4ProductTransform(newTransA.v, newTransA.v, hull_b->transformInv->m);
					a3real3Add(newTransA.v, hull_a->transform->v3.v);

					minA.x = newTransA.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = newTransA.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = newTransA.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = newTransA.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = newTransA.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = newTransA.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_halfwidth];
					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = hull_b->transform->v3.z;
					maxB.z = hull_b->transform->v3.z;

					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
					if (status == 0) break;

					a3vec3 newTransB;
					a3real3Diff(newTransB.v, hull_a->transform->v3.v, hull_b->transform->v3.v);
					a3real4ProductTransform(newTransB.v, newTransB.v, hull_a->transformInv->m);
					a3real3Add(newTransB.v, hull_b->transform->v3.v);

					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_halfwidth];
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_halfwidth];
					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_halfheight];
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_halfheight];
					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_halfdepth];
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_halfdepth];
					minB.x = newTransB.x - hull_b->prop[a3hullProperty_halfwidth];
					maxB.x = newTransB.x + hull_b->prop[a3hullProperty_halfwidth];
					minB.y = newTransB.y - hull_b->prop[a3hullProperty_halfheight];
					maxB.y = newTransB.y + hull_b->prop[a3hullProperty_halfheight];
					minB.z = hull_b->transform->v3.z;
					maxB.z = hull_b->transform->v3.z;

					status = a3collisionTestAABBs(collision_out, minA.v, maxA.v, minB.v, maxB.v, tmp);
				}
			}
			break;
			default:
				break;
			}
		}
		break;
		default:
			break;
		}

		if (status)
		{
			collision_out->hull_a = hull_a;
			collision_out->hull_b = hull_b;
		}

		return status;
	}
	return -1;
}


//-----------------------------------------------------------------------------
