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
	
	main_win.c
	Win32 application entry point (main).
*/


#if (defined _WINDOWS || defined _WIN32)


#include "a3_app_utils/Win32/a3_app_application.h"
#include "a3_app_utils/Win32/a3_app_window.h"

#include "../animal3D-LaunchApp.rc.h"


//-----------------------------------------------------------------------------
// main

int APIENTRY wWinMain(
	_In_ HINSTANCE hInstance,
	_In_opt_ HINSTANCE hPrevInstance,
	_In_ LPWSTR    lpCmdLine,
	_In_ int       nCmdShow)
{
	a3_WindowInterface wnd = { 0 };
	a3_WindowClass wndClass = { 0 };
	a3_RenderingContext_OpenGL renderContext = 0;
	const char *wndClassName = "A3_DEMO_CLASS";
	const char *wndName = "animal3D Demo";
	const unsigned int winWidth = 1024, winHeight = 640;
	int status;

	// initialize app
	status = a3appStartSingleInstanceSwitchExisting(wndClassName, wndName);

	// register window class
	status = a3appCreateDefaultRenderWindowClass(&wndClass, hInstance, wndClassName, sizeof(void *), IDI_ICON1, IDI_ICON1);
	if (status > 0)
	{
		// create rendering context
		status = a3appCreateDefaultRenderingContext(&renderContext, &wndClass);
		if (status > 0)
		{
			// create window
			status = a3appCreateWindow(&wnd, &wndClass, &renderContext, wndName, winWidth, winHeight, 0);
			if (status > 0)
			{
				// main loop
				status = a3appBeginWindowMainLoop(&wnd);
			}

			// kill rendering context
			status = a3appReleaseRenderingContext(&renderContext);
		}

		// kill window class
		status = a3appReleaseWindowClass(&wndClass, hInstance);
	}

	// the end
	return 0;
}


#endif	// (defined _WINDOWS || defined _WIN32)