/*
 * IDs: 0955181 and 0967813
 * EGP 425-01 Project 3 4/10/18
 * We certify that this work is entirely our own.  The assessor of this project may reproduce this project and provide copies to other academic staff, and/or communicate a copy of this project to a plagiarism-checking service, which may retain a copy of the project on its database.
 */

#ifndef BSP_H
#define BSP_H

#define RB_MAX 128

#include "physics\a3_Collision.h"


typedef struct BSP BSP;

// its name says bsp but its functionality says collision island
struct BSP
{

	a3_ConvexHull* containedHulls[RB_MAX];
	unsigned int numContainedHulls;
	a3vec3 min, max;

};

#endif // !BSP_H
