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
	
	a3_Integration.h
	Utility for handling complex general integration algorithms.
*/

#ifndef __ANIMAL3D_INTEGRATION_H
#define __ANIMAL3D_INTEGRATION_H


//-----------------------------------------------------------------------------

// math library
#include "animal3D/a3math/A3DM.h"


//-----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_Function		a3_Function;
#endif	// __cplusplus
	
	
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

	// function
	struct a3_Function
	{
		// ****TO-DO: 
		//	- add value, derivative function and integrator
		//	- add other utilities

		// y(t): value
		a3real *value;
	};


//-----------------------------------------------------------------------------

	// integration evaluators and helpers
	// ****TO-DO: 
	//	- implement some basic example integration methods
	//	- implement utilities to help with more advanced methods


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_INTEGRATION_H