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
	
	a3_app_window.h
	Interface for window creation and event loop.

	**DO NOT MODIFY THIS FILE**
*/

#if (defined _WINDOWS || defined _WIN32)
#ifndef __ANIMAL3D_APP_WINDOW_H
#define __ANIMAL3D_APP_WINDOW_H


#include <Windows.h>

#include "a3_app_utils/a3_app_load.h"


//-----------------------------------------------------------------------------
// structures

// A3: More self-explanitory aliases for existing structures.
typedef WNDCLASSEXA			a3_WindowClass;
typedef HWND				a3_WindowHandle;
typedef HDC					a3_DeviceContext;
typedef HGLRC				a3_RenderingContext_OpenGL;
typedef HMENU				a3_WindowMenuHandle;
typedef TRACKMOUSEEVENT		a3_MouseTracker[1];


// A3: Window interface.
// **DO NOT MODIFY THIS STRUCTURE**
typedef struct tagA3_WINDOWINTERFACE
{
	// platform-specifics
	a3_WindowHandle windowHandle;
	a3_DeviceContext deviceContext;
	a3_RenderingContext_OpenGL renderingContext;
	int windowID;
	int isStandalone;

	// platform-specific demo features
	a3_WindowMenuHandle windowMenuHandle, fileMenuHandle;
	a3_MouseTracker mouseTracker;

	// shared
	a3_DemoWindowInterface demo[1];
} a3_WindowInterface;


//-----------------------------------------------------------------------------
// functions

// A3: Create window class for application rendering window.
//	param windowClass_out: pointer to window class object to be created
//	param appInstance: handle to application
//	param windowClassName: name of application's window class
//	param extraBytes: extra space to reserve
//	param iconID: resource ID for a taskbar icon; pass -1 to ignore
//	param smallIconID: resource ID for small window icon; pass -1 to ignore
//	return: 1 if success (window class created)
//	return: 0 if fail (window class not created)
//	return: -1 if invalid params (null pointers)
int a3appCreateDefaultRenderWindowClass(a3_WindowClass *windowClass_out, HINSTANCE appInstance, const char *windowClassName, const unsigned int extraBytes, const int iconID, const int smallIconID);

// A3: Release application window class.
//	param windowClass: pointer to window class object to be released
//	param appInstance: handle to application
//	return: 1 if success (window class released)
//	return: 0 if fail (window class not released)
//	return: -1 if invalid params (null pointers)
int a3appReleaseWindowClass(a3_WindowClass *windowClass, HINSTANCE appInstance);

// A3: Create default OpenGL rendering context.
//	param context_out: pointer to context to be created
//	param windowClass: window class used to describe render window
//	return: 1 if success (created context)
//	return: 0 if fail (context not created)
//	return: -1 if invalid params (null pointers)
int a3appCreateDefaultRenderingContext(a3_RenderingContext_OpenGL *context_out, a3_WindowClass *windowClass);

// A3: Allow two rendering contexts to share resources. This must be done 
//		before starting to actually use either context to create GL resources.
//	param contextA: pointer to first context
//	param contextB: pointer to second context
//	return: 1 if success
//	return: 0 if failed
//	return: -1 if invalid params (null pointers)
int a3appShareRenderingContext(a3_RenderingContext_OpenGL *contextA, a3_RenderingContext_OpenGL *contextB);

// A3: Release OpenGL rendering context properly.
//	param context: pointer to context to be released
//	return: 1 if success (released context)
//	return: 0 if fail (context not released)
//	return: -1 if invalid params (null pointers)
int a3appReleaseRenderingContext(a3_RenderingContext_OpenGL *context);

// A3: Create window.
//	param windowHandle_out: pointer to hold resulting window structure
//	param windowClass: pointer to window class used to describe window
//	param renderContext: pointer to rendering context to use for this window
//	param windowName: name of window (appears at the top of the window)
//	params winWidth, winHeight: the size of the window in pixels (not including border)
//	param standalone: create standalone, single-demo window without menu
//	return: 1 if success (window created)
//	return: 0 if fail (not created)
//	return: -1 if invalid params (null pointers)
int a3appCreateWindow(a3_WindowInterface *window_out, a3_WindowClass *windowClass, a3_RenderingContext_OpenGL *renderContext, const char *windowName, const unsigned int winWidth, const unsigned int winHeight, const int standalone);

// A3: Release window.
//	param windowHandle: handle of window to be released
//	return: 1 if success (window released)
//	return: 0 if fail (not released)
//	return: -1 if invalid params (null pointer)
int a3appReleaseWindow(a3_WindowInterface *window);

// A3: Begin processing window messages.
//	param window: pointer to window object to begin processing messages for
//	return: return code from window when it exits
int a3appBeginWindowMainLoop(a3_WindowInterface *window);

// A3: Begin processing window messages on a new thread.
//	param window: pointer to window object to begin processing messages for
//	param handle: pointer to thread handle
//	return: thread ID if successful
//	return: 0 if thread handle already used
//	return: -1 if invalid params (null pointers)
int a3appBeginWindowMainLoopThreaded(a3_WindowInterface *window, HANDLE *thread);

// A3: Get return value from thread.
//	param handle: pointer to thread handle.
//	return: return value from thread (do not return -1)
//	return: -1 if invalid param or thread is still running
int a3appGetThreadReturnValue(HANDLE *handle);


//-----------------------------------------------------------------------------


#endif	// !__ANIMAL3D_APP_WINDOW_H
#endif	// (defined _WINDOWS || defined _WIN32)