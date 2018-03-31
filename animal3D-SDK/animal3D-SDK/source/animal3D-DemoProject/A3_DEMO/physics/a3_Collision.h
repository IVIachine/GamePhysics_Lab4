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
	
	a3_Collision.h
	Collision detection and response algorithms.
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

#ifndef __ANIMAL3D_COLLISION_H
#define __ANIMAL3D_COLLISION_H


// rigid body
#include "a3_RigidBody.h"
#include "a3_Ray.h"


//-----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_ConvexHull			a3_ConvexHull;
	typedef struct a3_ConvexHullCollision	a3_ConvexHullCollision;
	typedef enum a3_ConvexHullType			a3_ConvexHullType;
	typedef enum a3_ConvexHullFlag			a3_ConvexHullFlag;
	typedef enum a3_ConvexHullProperty		a3_ConvexHullProperty;
#endif	// __cplusplus

	
//-----------------------------------------------------------------------------

	// constants
	enum a3_ConvexHullLimits
	{
		a3hullProperty_maxCount_preset = 16,
		a3hullProperty_maxCount_user = 0,

		a3hullContact_maxCount = 8,
	};

	// generic convex hull types
	enum a3_ConvexHullType
	{
		a3hullType_none,
		a3hullType_point,
		a3hullType_plane,
		a3hullType_box,
		a3hullType_sphere,
		a3hullType_cylinder,
		a3hullType_mesh,
	};

	// generic convex hull flags
	enum a3_ConvexHullFlag
	{
		a3hullFlag_none,			// 0x00
		a3hullFlag_is3D,			// 0x01
		a3hullFlag_isAxisAligned,	// 0x02
	};

	// generic convex hull properties
	enum a3_ConvexHullProperty
	{
		// boxy shape properties
		a3hullProperty_width = 0,
		a3hullProperty_height,
		a3hullProperty_depth,
		a3hullProperty_halfwidth,
		a3hullProperty_halfheight,
		a3hullProperty_halfdepth,
		a3hullProperty_halfwidthSq,
		a3hullProperty_halfheightSq,
		a3hullProperty_halfdepthSq,

		// round shape properties
		a3hullProperty_radius = 0,
		a3hullProperty_length,
		a3hullProperty_radiusSq,

		// user-defined properties
		a3hullProperty_user = a3hullProperty_maxCount_preset,
	};

	// generic convex hull descriptor
	struct a3_ConvexHull
	{
		// ****TO-DO: 
		//	- add pertinent data
		const a3_RigidBody* rb;
		const a3mat4* transform, *transformInv;

		a3_ConvexHullType type;
		a3_ConvexHullFlag flag;

		a3_Axis axis;

		// hull shape properties
		a3real prop[a3hullProperty_maxCount_preset + a3hullProperty_maxCount_user];
	};

	// collision descriptor
	struct a3_ConvexHullCollision
	{
		// ****TO-DO: 
		//	- add pertinent data

		const a3_ConvexHull *hull_a, *hull_b;

		// list of contact points + normals
		a3vec3 contact_a[a3hullContact_maxCount], contact_b[a3hullContact_maxCount];
		a3vec3 normal_a[a3hullContact_maxCount], normal_b[a3hullContact_maxCount];
		unsigned int contactCount_a, contactCount_b;
	};


//-----------------------------------------------------------------------------

	// create point hull
	inline int a3collisionCreateHullPoint(a3_ConvexHull *hull_out, const a3_RigidBody *rb);

	// create plane hull
	inline int a3collisionCreateHullPlane(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const a3real width, const a3real height, const int isAxisAligned, const a3_Axis normalAxis);

	// create box hull
	inline int a3collisionCreateHullBox(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const a3real width, const a3real height, const a3real depth, const int isAxisAligned);

	// create sphere hull
	inline int a3collisionCreateHullSphere(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const a3real radius);

	// create cylinder hull
	inline int a3collisionCreateHullCylinder(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const a3real radius, const a3real length, const a3_Axis normalAxis);

	// create mesh hull
	inline int a3collisionCreateHullMesh(a3_ConvexHull *hull_out, const a3_RigidBody *rb, const a3mat4 *transform, const a3mat4 *transformInv, const void *points, const unsigned int pointCount, const int is3D);


//-----------------------------------------------------------------------------

	// high-level collision test
	inline int a3collisionTestConvexHulls(a3_ConvexHullCollision *collision_out, const a3_ConvexHull *hull_a, const a3_ConvexHull *hull_b);
	

//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_COLLISION_H