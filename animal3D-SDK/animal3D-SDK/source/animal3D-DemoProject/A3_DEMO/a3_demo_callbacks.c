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
	
	a3_demo_callbacks.c/.cpp
	Main implementation of hotloaded callbacks.

	********************************************
	*** THIS IS THE LIBRARY'S CALLBACKS FILE ***
	*** App hooks your demo via this file.   ***
	********************************************
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

#include "a3_dylib_config_export.h"
#include "a3_DemoState.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//-----------------------------------------------------------------------------
// miscellaneous functions

// get the size of the persistent state to allocate
//	(good idea to keep it relatively constant, so that addresses don't change 
//	when the library is reloaded... that would mess everything up!)
inline unsigned int a3demo_getPersistentStateSize()
{
	const unsigned int minimum = sizeof(a3_DemoState);
	unsigned int size = 1;
	while (size < minimum)
		size += size;
	return size;
}


// consistent text initialization
inline void a3demo_initializeText(a3_TextRenderer *text)
{
	a3textInitialize(text, 20, 1, 0, 0, 0);
}


// consistent physics thread startup and shutdown
inline void a3demo_initializePhysicsThread(a3_DemoState *demoState)
{
	a3real tmp0;
	unsigned int i;

	// reset world and launch thread
	demoState->physicsWorld->lock = 0;
	demoState->physicsWorld->init = 0;
	demoState->physicsWorld->rate = 100.0;
	a3threadLaunch(demoState->physicsThread, (a3_threadfunc)a3physicsThread, demoState->physicsWorld, "a3physicsThread");

	// wait for thread to initialize
	while (!demoState->physicsWorld->init);

	// set object scales based on hull properties
	// for ground, can set by reading first two values of properties list (width, height)
	for (i = 0; i < 6; ++i)
	{
		a3vec3 tmp;
		tmp.x = 60.0f;
		tmp.y = 60.0f;
		tmp.z = 60.0f;

		a3real3SetReal2Z(demoState->groundObject[i].scale.v, tmp.v, a3realOne);
		demoState->groundObject[i].scaleMode = -1;
	}
	for (i = 0; i < 7; ++i)
	{
		// uniform scale: radius
		tmp0 = demoState->physicsWorld->hull_sphere[i].prop[a3hullProperty_radius];
		a3real3Set(demoState->sphereObject[i].scale.v, tmp0, tmp0, tmp0);
		demoState->sphereObject[i].scaleMode = +1;
	}
}

inline void a3demo_terminatePhysicsThread(a3_DemoState *demoState)
{
	demoState->physicsWorld->lock = -1;
	a3threadWait(demoState->physicsThread);
}


//-----------------------------------------------------------------------------
// event handling

// raypick test
inline void a3demo_raypickTest_internal(a3_DemoState *demoState)
{
	//a3_RayHit hit[1] = { 0 };
	//hit->param0 = 100.0f;
	//
	//*demoState->rayHit = *hit;
	//demoState->hitIndex = -1;
	//
	//
	//// test against objects in order, winner is closest
	//
	//
	//if (a3rayTestSphere(hit, demoState->ray, 
	//	demoState->sphereObject->scale.x, 
	//	demoState->sphereObject->modelMat.m)
	//	&& a3rayHitValidate(hit))
	//{
	//	if (hit->param0 < demoState->rayHit->param0)
	//	{
	//		*demoState->rayHit = *hit;
	//		demoState->hitIndex = 1;
	//		printf("Sphere\n");
	//	}
	//}
	//
	//if (a3rayTestCylinderFinite(hit, demoState->ray, a3axis_z, 
	//	demoState->cylinderObject->scale.x, demoState->cylinderObject->scale.z, 
	//	demoState->cylinderObject->modelMat.m)
	//	&& a3rayHitValidate(hit))
	//{
	//	if (hit->param0 < demoState->rayHit->param0)
	//	{
	//		*demoState->rayHit = *hit;
	//		demoState->hitIndex = 2;
	//		printf("Cylinder\n");
	//	}
	//}
	//
	//if (a3rayTestBoundingBox(hit, demoState->ray, 
	//	demoState->boxObject->scale.x, demoState->boxObject->scale.y, demoState->boxObject->scale.z,
	//	demoState->boxObject->modelMat.m, demoState->boxObject->modelMatInv.m)
	//	&& a3rayHitValidate(hit))
	//{
	//	if (hit->param0 < demoState->rayHit->param0)
	//	{
	//		*demoState->rayHit = *hit;
	//		demoState->hitIndex = 3;
	//		printf("Bounding Box\n");
	//	}
	//}

	//if (a3rayTestBoundingBox(hit, demoState->ray,
	//	demoState->boxObject->scale.x, demoState->boxObject->scale.y, demoState->boxObject->scale.z,
	//	demoState->boxObject->modelMat.m, demoState->boxObject->modelMatInv.m)
	//	&& a3rayHitValidate(hit))
	//{
	//	if (hit->param0 < demoState->rayHit->param0)
	//	{
	//		*demoState->rayHit = *hit;
	//		demoState->hitIndex = 3;
	//	}
	//}

	
}



//-----------------------------------------------------------------------------
// callback prototypes
// NOTE: do not move to header; they should be private to this file
// NOTE: you may name these functions whatever you like, just be sure to 
//	update your debug config file: 
//	"<root>/resource/animal3D-data/animal3D-demoinfo-debug.txt"
// copy this config line and the DLL to your main config with a new name when 
//	you're happy with it: 
//	"<root>/resource/animal3D-data/animal3D-demoinfo.txt"


#ifdef __cplusplus
extern "C"
{
#endif	// __cplusplus

	A3API a3_DemoState *a3test_load();
	A3API void a3test_unload(a3_DemoState *demoState, int hotload);
	A3API a3_DemoState *a3test_hotload(a3_DemoState *demoState);
	A3API int a3test_display(a3_DemoState *demoState);
	A3API int a3test_idle(a3_DemoState *demoState);
	A3API void a3test_windowActivate(a3_DemoState *demoState);
	A3API void a3test_windowDeactivate(a3_DemoState *demoState);
	A3API void a3test_windowMove(a3_DemoState *demoState, int newWindowPosX, int newWindowPosY);
	A3API void a3test_windowResize(a3_DemoState *demoState, int newWindowWidth, int newWindowHeight);
	A3API void a3test_keyPress(a3_DemoState *demoState, int virtualKey);
	A3API void a3test_keyHold(a3_DemoState *demoState, int virtualKey);
	A3API void a3test_keyRelease(a3_DemoState *demoState, int virtualKey);
	A3API void a3test_keyCharPress(a3_DemoState *demoState, int asciiKey);
	A3API void a3test_keyCharHold(a3_DemoState *demoState, int asciiKey);
	A3API void a3test_mouseClick(a3_DemoState *demoState, int button, int cursorX, int cursorY);
	A3API void a3test_mouseDoubleClick(a3_DemoState *demoState, int button, int cursorX, int cursorY);
	A3API void a3test_mouseRelease(a3_DemoState *demoState, int button, int cursorX, int cursorY);
	A3API void a3test_mouseWheel(a3_DemoState *demoState, int delta, int cursorX, int cursorY);
	A3API void a3test_mouseMove(a3_DemoState *demoState, int cursorX, int cursorY);
	A3API void a3test_mouseLeave(a3_DemoState *demoState);

#ifdef __cplusplus
}
#endif	// __cplusplus


//-----------------------------------------------------------------------------
// callback implementations

// demo is loaded
A3API a3_DemoState *a3test_load()
{
	// HEAP allocate persistent state
	// stack object will be deleted at the end of the function
	// good idea to set the whole block of memory to zero
	const unsigned int stateSize = a3demo_getPersistentStateSize();
	a3_DemoState *demoState = (a3_DemoState *)malloc(stateSize);
	memset(demoState, 0, stateSize);

	// set up trig table (A3DM)
	a3trigInit(4, demoState->trigTable);

	// initialize state variables
	// e.g. timer, thread, etc.
	a3timerSet(demoState->renderTimer, 30.0);
	a3timerStart(demoState->renderTimer);

	// text
	a3demo_initializeText(demoState->text);
	demoState->textInit = demoState->showText = 1;


	// use Y-"up"
//	demoState->verticalAxis = 1;

	// enable asset streaming between loads
//	demoState->streaming = 1;


	// set default GL state
	a3demo_setDefaultGraphicsState();

	// textures
	a3demo_loadTextures(demoState);

	// geometry
	a3demo_loadGeometry(demoState);

	// shaders
	a3demo_loadShaders(demoState);

	// scene objects
	a3demo_initScene(demoState);


	// launch physics simulation thread
	a3demo_initializePhysicsThread(demoState);


	// return persistent state pointer
	return demoState;
}

// demo is unloaded; option to unload to prep for hotload
A3API void a3test_unload(a3_DemoState *demoState, int hotload)
{
	// release things that need releasing always, whether hotloading or not
	// e.g. kill thread
	// nothing in this example, but then...

	// release persistent state if not hotloading
	// good idea to release in reverse order that things were loaded...
	//	...normally; however, in this case there are reference counters 
	//	to make sure things get deleted when they need to, so there is 
	//	no need to reverse!
	if (!hotload)
	{
		// free fixed objects
		a3textRelease(demoState->text);

		// free graphics objects
		a3demo_unloadTextures(demoState);
		a3demo_unloadGeometry(demoState);
		a3demo_unloadShaders(demoState);

		// validate unload
		a3demo_validateUnload(demoState);

		// erase other stuff
		a3trigFree();

		// wait for physics thread to stop
		a3demo_terminatePhysicsThread(demoState);

		// erase persistent state
		free(demoState);
	}
}

// demo is hotloaded
A3API a3_DemoState *a3test_hotload(a3_DemoState *demoState)
{
	// initial allocation if applicable
	if (!demoState)
		demoState = a3test_load();

	// do any re-allocation tasks
	else
	{
		const unsigned int stateSize = a3demo_getPersistentStateSize();
		a3_DemoState copy = *demoState;

		// wait for physics thread to stop
		a3demo_terminatePhysicsThread(demoState);

		// example 1: copy memory directly
		free(demoState);
		demoState = (a3_DemoState *)malloc(stateSize);
		memset(demoState, 0, stateSize);
		*demoState = copy;

		a3demo_refresh(demoState);
		a3trigInitSetTables(4, demoState->trigTable);

		// restart physics thread
		a3demo_initializePhysicsThread(demoState);
	}

	// return pointer to new persistent state
	return demoState;
}

// window updates display
// NOTE: DO NOT USE FOR RENDERING
A3API int a3test_display(a3_DemoState *demoState)
{
	// do nothing, should just return 1 to indicate that the 
	//	window's display area is controlled by the demo
	return 1;
}

// window idles
A3API int a3test_idle(a3_DemoState *demoState)
{
	// perform any idle tasks, such as rendering
	if (!demoState->exitFlag)
	{
		if (a3timerUpdate(demoState->renderTimer) > 0)
		{
			// render timer ticked, update demo state and draw
			a3demo_update(demoState, demoState->renderTimer->secondsPerTick);
			a3demo_input(demoState, demoState->renderTimer->secondsPerTick);
			a3demo_render(demoState);

			// update input
			a3mouseUpdate(demoState->mouse);
			a3keyboardUpdate(demoState->keyboard);
			a3XboxControlUpdate(demoState->xcontrol);

			// render occurred this idle: return +1
			return +1;
		}

		// nothing happened this idle: return 0
		return 0;
	}

	// demo should exit now: return -1
	return -1;
}

// window gains focus
A3API void a3test_windowActivate(a3_DemoState *demoState)
{
	// nothing really needs to be done here...
	//	but it's here just in case
}

// window loses focus
A3API void a3test_windowDeactivate(a3_DemoState *demoState)
{
	// reset input; it won't track events if the window is inactive, 
	//	active controls will freeze and you'll get strange behaviors
	a3keyboardReset(demoState->keyboard);
	a3mouseReset(demoState->mouse);
	a3XboxControlReset(demoState->xcontrol);
	a3XboxControlSetRumble(demoState->xcontrol, 0, 0);
}

// window moves
A3API void a3test_windowMove(a3_DemoState *demoState, int newWindowPosX, int newWindowPosY)
{
	// nothing needed here
}

// window resizes
A3API void a3test_windowResize(a3_DemoState *demoState, int newWindowWidth, int newWindowHeight)
{
	// account for borders here
	const int frameBorder = 0;
	const unsigned int frameWidth = newWindowWidth + frameBorder + frameBorder;
	const unsigned int frameHeight = newWindowHeight + frameBorder + frameBorder;

	// copy new values to demo state
	demoState->frameBorder = frameBorder;
	demoState->windowWidth = newWindowWidth;
	demoState->windowHeight = newWindowHeight;
	demoState->frameWidth = frameWidth;
	demoState->frameHeight = frameHeight;

	// use framebuffer deactivate utility to set viewport
	a3framebufferDeactivateSetViewport(a3fbo_depthDisable, -frameBorder, -frameBorder, demoState->frameWidth, demoState->frameHeight);

	// viewing info for projection matrix
	demoState->camera->aspect = (a3real)frameWidth / (a3real)frameHeight;

	// initialize cameras dependent on viewport
	a3real4x4MakePerspectiveProjection(demoState->camera->projectionMat.m, demoState->camera->projectionMatInv.m,
		demoState->camera->fovy, demoState->camera->aspect, demoState->camera->znear, demoState->camera->zfar);
}

// any key is pressed
A3API void a3test_keyPress(a3_DemoState *demoState, int virtualKey)
{
	// persistent state update
	a3keyboardSetState(demoState->keyboard, (a3_KeyboardKey)virtualKey, a3input_down);
}

// any key is held
A3API void a3test_keyHold(a3_DemoState *demoState, int virtualKey)
{
	// persistent state update
	a3keyboardSetState(demoState->keyboard, (a3_KeyboardKey)virtualKey, a3input_down);
}

// any key is released
A3API void a3test_keyRelease(a3_DemoState *demoState, int virtualKey)
{
	// persistent state update
	a3keyboardSetState(demoState->keyboard, (a3_KeyboardKey)virtualKey, a3input_up);
}

// ASCII key is pressed (immediately preceded by "any key" pressed call above)
// NOTE: there is no release counterpart
A3API void a3test_keyCharPress(a3_DemoState *demoState, int asciiKey)
{
	// persistent state update
	a3keyboardSetStateASCII(demoState->keyboard, (char)asciiKey);

	// handle special cases immediately
	switch (asciiKey)
	{
		// uncomment to make escape key kill the current demo
		// if disabled, use 'exit demo' menu option
//	case 27: 
//		demoState->exitFlag = 1;
//		break;

		// reload (T) or toggle (t) text
	case 'T': 
		if (!a3textIsInitialized(demoState->text))
		{
			a3demo_initializeText(demoState->text);
			demoState->textInit = 1;
		}
		else
		{
			a3textRelease(demoState->text);
			demoState->textInit = 0;
		}
		break;
	case 't':
		demoState->showText = 1 - demoState->showText;
		break;

		// reload all shaders in real-time
	case 'P': 
		a3demo_unloadShaders(demoState);
		a3demo_loadShaders(demoState);
		break;

		// reload physics thread
	case 'R':
		a3demo_terminatePhysicsThread(demoState);
		a3demo_initializePhysicsThread(demoState);
		break;

		// change modes and toggles
	case '.':
		demoState->demoMode = (demoState->demoMode + 1) % demoState->demoModeCount;
		break;
	case ',':
		demoState->demoMode = (demoState->demoMode + demoState->demoModeCount - 1) % demoState->demoModeCount;
		break;


		// other toggles
	case 'x':
		demoState->displayAxes = 1 - demoState->displayAxes;
		break;
	case 'g':
		demoState->displayGrid = 1 - demoState->displayGrid;
		break;
	case 'b':
		demoState->displaySkybox = 1 - demoState->displaySkybox;
		break;
	case 'c':
		demoState->displayPhysicsText = 1 - demoState->displayPhysicsText;
		break;
	}
}

// ASCII key is held
A3API void a3test_keyCharHold(a3_DemoState *demoState, int asciiKey)
{
	// persistent state update
	a3keyboardSetStateASCII(demoState->keyboard, (char)asciiKey);
}

// mouse button is clicked
A3API void a3test_mouseClick(a3_DemoState *demoState, int button, int cursorX, int cursorY)
{
	// persistent state update
	a3mouseSetState(demoState->mouse, (a3_MouseButton)button, a3input_down);
	a3mouseSetPosition(demoState->mouse, cursorX, cursorY);

	// on right click, cast ray into scene
	if (button == a3mouse_right)
	{
		// position in screen-space converted to NDC
		const a3vec3 ndc = {
			+a3deserialize((float)cursorX / (float)demoState->windowWidth),
			-a3deserialize((float)cursorY / (float)demoState->windowHeight),
			-a3realOne
		};

		// un-project ray into scene
		a3rayCreateUnprojected(demoState->ray, ndc, demoState->camera->projectionMatInv.m);
		a3rayTransform(demoState->ray, demoState->ray, demoState->cameraObject->modelMat.m);

		// perform raypick test
		a3demo_raypickTest_internal(demoState);
	}
}

// mouse button is double-clicked
A3API void a3test_mouseDoubleClick(a3_DemoState *demoState, int button, int cursorX, int cursorY)
{
	// persistent state update
	a3mouseSetState(demoState->mouse, (a3_MouseButton)button, a3input_down);
	a3mouseSetPosition(demoState->mouse, cursorX, cursorY);
}

// mouse button is released
A3API void a3test_mouseRelease(a3_DemoState *demoState, int button, int cursorX, int cursorY)
{
	// persistent state update
	a3mouseSetState(demoState->mouse, (a3_MouseButton)button, a3input_up);
	a3mouseSetPosition(demoState->mouse, cursorX, cursorY);
}

// mouse wheel is turned
A3API void a3test_mouseWheel(a3_DemoState *demoState, int delta, int cursorX, int cursorY)
{
	// persistent state update
	a3mouseSetStateWheel(demoState->mouse, (a3_MouseWheelState)delta);
	a3mouseSetPosition(demoState->mouse, cursorX, cursorY);

	// can use this to change zoom
	// zoom should be faster farther away
	demoState->camera->fovy -= demoState->camera->ctrlZoomSpeed * (demoState->camera->fovy / a3realOneEighty) * (float)delta;
	demoState->camera->fovy = a3clamp(demoState->camera->ctrlZoomSpeed, a3realOneEighty - demoState->camera->ctrlZoomSpeed, demoState->camera->fovy);

	a3real4x4MakePerspectiveProjection(demoState->camera->projectionMat.m, demoState->camera->projectionMatInv.m,
		demoState->camera->fovy, demoState->camera->aspect, demoState->camera->znear, demoState->camera->zfar);
}

// mouse moves
A3API void a3test_mouseMove(a3_DemoState *demoState, int cursorX, int cursorY)
{
	// persistent state update
	a3mouseSetPosition(demoState->mouse, cursorX, cursorY);
}

// mouse leaves window
A3API void a3test_mouseLeave(a3_DemoState *demoState)
{
	// reset mouse state or any buttons pressed will freeze
	a3mouseReset(demoState->mouse);
}


//-----------------------------------------------------------------------------
