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
	// ****TO-DO: 
	//	- implement test

	return 0;
}

inline int a3collisionTestPlaneSphere(
	const a3real3p planeCenter, const a3real3p planeTangent, const a3real3p planeBitangent, const a3real3p planeNormal, const a3real planeHalfWidth, const a3real planeHalfHeight, const a3real3p sphereCenter, const a3real sphereRadius, a3real3p diff_tmp)
{
	// ****TO-DO: 
	//	- implement test

	return 0;
}

inline int a3collisionTestPlaneAABB(
	const a3real3p planeCenter_localToAABB, const a3real3p planeTangent_localToAABB, const a3real3p planeBitangent_localToAABB, const a3real3p planeNormal_localToAABB, const a3real planeHalfWidth, const a3real planeHalfHeight, const a3real3p aabbMinExtents, const a3real3p aabbMaxExtents, a3real3p diff_tmp)
{
	// ****TO-DO: 
	//	- implement test

	return 0;
}

inline int a3collisionTestSpheres(
	const a3real3p sphereCenter_a, const a3real sphereRadius_a, const a3real3p sphereCenter_b, const a3real sphereRadius_b, a3real3p diff_tmp)
{
	// ****TO-DO: 
	//	- implement test
	const a3real sumRadii = sphereRadius_a + sphereRadius_b;
	a3real3Diff(diff_tmp, sphereCenter_a, sphereCenter_b);
	return (a3real3LengthSquared(diff_tmp) <= sumRadii * sumRadii);
}

inline int a3collisionTestSphereAABB(
	const a3real3p sphereCenter_localToAABB, const a3real sphereRadius, const a3real3p aabbMinExtents, const a3real3p aabbMaxExtents, a3real3p diff_tmp)
{
	// ****TO-DO: 
	//	- implement test

	return 0;
}

inline int a3collisionTestAABBs(
	const a3real3p aabbMinExtents_a, const a3real3p aabbMaxExtents_a, const a3real3p aabbMinExtents_b, const a3real3p aabbMaxExtents_b, a3real3p diff_tmp)
{
	// ****TO-DO: 
	//	- implement test
	if ((aabbMinExtents_a[0] <= aabbMaxExtents_b[0] && aabbMaxExtents_a[0] >= aabbMinExtents_b[0]) &&
		(aabbMinExtents_a[1] <= aabbMaxExtents_b[1] && aabbMaxExtents_a[1] >= aabbMinExtents_b[1]) &&
		(aabbMinExtents_a[2] <= aabbMaxExtents_b[2] && aabbMaxExtents_a[2] >= aabbMinExtents_b[2]))
	{
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

// high-level collision test
extern inline int a3collisionTestConvexHulls(a3_ConvexHullCollision *collision_out, const a3_ConvexHull *hull_a, const a3_ConvexHull *hull_b)
{
	if (collision_out && hull_a && hull_b)
	{
		// ****TO-DO: 
		//	- select appropriate internal test
		//	- do necessary preparations for selected test (e.g. transformations)
		//	- do multiple tests if necessary (e.g. two arbitrary boxes)

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
				status = a3collisionTestSpheres(hull_a->transform->v3.v, hull_a->prop[a3hullProperty_radius],
					hull_b->transform->v3.v, hull_b->prop[a3hullProperty_radius], tmp);
			}
			break;
			case a3hullType_box:
			{
				//if AA
				//calc extents perform AABB test
				//else
				//transform sphere to box space, perform AABB test
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
				//Same as one above
			}
			break;
			case a3hullType_box:
			{
				a3vec3 minA, maxA, minB, maxB;
				if (hull_a->prop[a3hullFlag_isAxisAligned] == 2 && hull_b->prop[a3hullFlag_isAxisAligned] == 2)
				{
					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_width] * a3realHalf;
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_width] * a3realHalf;

					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_height] * a3realHalf;
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_height] * a3realHalf;

					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_depth] * a3realHalf;
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_depth] * a3realHalf;

					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_width] * a3realHalf;
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_width] * a3realHalf;
					   			  							
					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_height] * a3realHalf;
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_height] * a3realHalf;
					   			  							
					minB.z = hull_b->transform->v3.z - hull_b->prop[a3hullProperty_depth] * a3realHalf;
					maxB.z = hull_b->transform->v3.z + hull_b->prop[a3hullProperty_depth] * a3realHalf;

					status = a3collisionTestAABBs(minA.v, maxA.v, minB.v, maxB.v, tmp);
				}
				else if (hull_a->prop[a3hullFlag_isAxisAligned] == 2 && hull_b->prop[a3hullFlag_isAxisAligned] != 2 ||
					hull_a->prop[a3hullFlag_isAxisAligned] != 2 && hull_b->prop[a3hullFlag_isAxisAligned] == 2)
				{
					a3mat4 transformedBox;
					a3real4x4SetReal4x4(transformedBox.m, hull_a->transform->m);
					a3real4x4MulTransform(transformedBox.m, hull_b->transformInv->m);

					minA.x = hull_a->transform->v3.x - hull_a->prop[a3hullProperty_width] * a3realHalf;
					maxA.x = hull_a->transform->v3.x + hull_a->prop[a3hullProperty_width] * a3realHalf;

					minA.y = hull_a->transform->v3.y - hull_a->prop[a3hullProperty_height] * a3realHalf;
					maxA.y = hull_a->transform->v3.y + hull_a->prop[a3hullProperty_height] * a3realHalf;

					minA.z = hull_a->transform->v3.z - hull_a->prop[a3hullProperty_depth] * a3realHalf;
					maxA.z = hull_a->transform->v3.z + hull_a->prop[a3hullProperty_depth] * a3realHalf;

					a3real4TransformProduct(minA.v, transformedBox.m, minA.v);
					a3real4TransformProduct(maxA.v, transformedBox.m, maxA.v);

					minB.x = hull_b->transform->v3.x - hull_b->prop[a3hullProperty_width] * a3realHalf;
					maxB.x = hull_b->transform->v3.x + hull_b->prop[a3hullProperty_width] * a3realHalf;

					minB.y = hull_b->transform->v3.y - hull_b->prop[a3hullProperty_height] * a3realHalf;
					maxB.y = hull_b->transform->v3.y + hull_b->prop[a3hullProperty_height] * a3realHalf;

					minB.z = hull_b->transform->v3.z - hull_b->prop[a3hullProperty_depth] * a3realHalf;
					maxB.z = hull_b->transform->v3.z + hull_b->prop[a3hullProperty_depth] * a3realHalf;

					status = a3collisionTestAABBs(minA.v, maxA.v, minB.v, maxB.v, tmp);

					if (status == 1)
						printf("HIT\n");
					
					/*if (status == 0)
						break;

					a3real4x4SetReal4x4(transformedBox.m, hull_b->transform->m);
					a3real4x4MulTransform(transformedBox.m, hull_a->transformInv->m);

					minA.x = hull_a->transform->v3.x;
					maxA.x = minA.x + hull_a->prop[a3hullProperty_width];

					minA.y = hull_a->transform->v3.y;
					maxA.y = minA.y + hull_a->prop[a3hullProperty_height];

					minA.z = hull_a->transform->v3.z;
					maxA.z = minA.z + hull_a->prop[a3hullProperty_depth];

					minB.x = transformedBox.v3.x;
					maxB.x = minB.x + hull_b->prop[a3hullProperty_width];

					minB.y = transformedBox.v3.y;
					maxB.y = minB.y + hull_b->prop[a3hullProperty_height];

					minB.z = transformedBox.v3.z;
					maxB.z = minB.z + hull_b->prop[a3hullProperty_depth];

					status = a3collisionTestAABBs(minA.v, maxA.v, minB.v, maxB.v, tmp);*/
				}
				else
				{

				}
				//if both AABB
				//Calc extents perform AABB test
				//if One AABB
				//calc extents on AABB, calc local extents on AABB, perform AABB test
					//transform AABB into non-AABB, calculate local extents of orig, perform AABB test
				//if neither
				//transform first into other, calc first local, perform AABB
					//transform other into first, calculate other local, perform AABB
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
