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
	
	a3_GraphicsObjectHandle.h
	Generic graphics object handle with reference counter. Provides the option 
		to auto-release a graphics object when the counter reaches zero.

	**DO NOT MODIFY THIS FILE**
*/

#ifndef __ANIMAL3D_GRAPHICSOBJECTHANDLE_H
#define __ANIMAL3D_GRAPHICSOBJECTHANDLE_H


#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_GraphicsObjectHandle	a3_GraphicsObjectHandle;
#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// A3: Aliases for auto-release function pointer types: 
	//	Custom: programmer-defined release function (good for debugging)
	typedef void(*a3_GraphicsObjectReleaseFunc)(int, unsigned int *);


	// A3: Generic reference-counting graphics object handle.
	//	member releaseFunc: function pointer to auto-release function
	//	member objectCount: how many objects were initialized under this handle
	//	member refCount: reference counter for this handle (i.e. the option to 
	//		keep track of how many other handles are dependent on this one)
	//	members handle, handlePtr: overlapping handle used by OpenGL, 
	//		represented as both an integer and a pointer (because array)
	//	NOTE: because the handle occurs at the end of the structure, any other 
	//		data type that declares a graphics handle instance may declare an 
	//		array of integers immediately following to store multiple handles; 
	//		e.g. a struct that has 4 graphics objects of the same type may 
	//		declare its handles like this: 
	//			a3_GraphicsObjectHandle handle[1];
	//			unsigned int handles_extra[3];
	struct a3_GraphicsObjectHandle
	{
		a3_GraphicsObjectReleaseFunc releaseFunc;
		int objectCount;
		int refCount;

		union {
			unsigned int handle;
			unsigned int handlePtr[1];
		};
	};


//-----------------------------------------------------------------------------
	
	// A3: Create generic handle with an custom object release function.
	//	param handle_out: non-null pointer to graphics object handle container
	//	param releaseFunc_opt: optional pointer to release function with 
	//		return type void, integer argumant and integer pointer argument
	//	param handleValue: non-zero initial value of handle
	//	param objectCount: number of handles to be released (not always used)
	//	return: handleValue if successfully created
	//	return: -1 if invalid params or if already in-use
	int a3handleCreateHandle(a3_GraphicsObjectHandle *handle_out, const a3_GraphicsObjectReleaseFunc releaseFunc_opt, const unsigned int handleValue, const unsigned int objectCount);

	// A3: Update handle delete callback with custom function.
	//	param handle: non-null pointer to graphics object handle container
	//	param releaseFunc_opt: optional pointer to release function with 
	//		return type void and one integer pointer argument
	//	return: 1 if successfully updated
	//	return: -1 if invalid params
	int a3handleSetReleaseFunc(a3_GraphicsObjectHandle *handle, const a3_GraphicsObjectReleaseFunc releaseFunc_opt);

	// A3: Increment handle counter.
	//	param handle: non-null pointer to graphics object handle container
	//	return: new count if success
	//	return: -1 if invalid param or internal handle is zero
	int a3handleIncrementCount(a3_GraphicsObjectHandle *handle);

	// A3: Decrement handle counter; calls release function if one is provided 
	//		and count hits zero
	//	param handle: non-null pointer to graphics object handle container
	//	return: new count (min 0) if success
	//	return: -1 if invalid param or internal handle is zero
	int a3handleDecrementCount(a3_GraphicsObjectHandle *handle);

	// A3: Force-release handle; set value to zero and call release function.
	//	param handle: non-null pointer to graphics object handle container
	//	return: 0 if success
	//	return: -1 if invalid param or internal handle is zero
	int a3handleRelease(a3_GraphicsObjectHandle *handle);


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#include "_inl/a3_GraphicsObjectHandle.inl"


#endif	// !__ANIMAL3D_GRAPHICSOBJECTHANDLE_H