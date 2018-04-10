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
