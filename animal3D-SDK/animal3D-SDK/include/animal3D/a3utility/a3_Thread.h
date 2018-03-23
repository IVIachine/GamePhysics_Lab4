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
	
	a3_Thread.h
	Basic thread starter interface.

	**DO NOT MODIFY THIS FILE**
*/

#ifndef __ANIMAL3D_THREAD_H
#define __ANIMAL3D_THREAD_H


#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_Thread a3_Thread;
#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// A3: Thread function alias; any function declared with this format may 
	//	be used to launch a thread: 
	//	-> returns value of type "long" (i.e. an integer)
	//	-> single parameter of type void pointer (strict)
	typedef long(*a3_threadfunc)(void *);


//-----------------------------------------------------------------------------

	// A3: Simple structure for thread descriptor
	//	member threadID: thread ID number that appears in debugger; non-zero 
	//		when initialized; resets to zero when thread terminates
	//	member result: value returned when thread exits
	//	member handle: internal thread handle(s)
	//	member threadName: optional name of thread, renamed on launch
	//	member threadArgs: pointer to pass to threadFunc on launch
	//	member threadFunc: function to launch when thread starts
	struct a3_Thread
	{
		unsigned long threadID;
		long result;
		void *handle[2];
		char *threadName;
		void *threadArgs;
		a3_threadfunc threadFunc;
	};


//-----------------------------------------------------------------------------

	// A3: Launch thread.
	//	param thread_out: pointer to unused thread descriptor
	//	param func: non-null pointer to function to launch when thread starts
	//	param args_opt: optional pointer to function arguments
	//	param name_opt: optional cstring to rename thread; default name if null
	//	return: 1 if success (thread started)
	//	return: 0 if fail (already initialized)
	//	return: -1 if invalid params (null pointers)
	int a3threadLaunch(a3_Thread *thread_out, a3_threadfunc func, void *args_opt, char *name_opt);

	// A3: Terminate a running thread immediately.
	//	param thread: container of thread to terminate
	//	return: 1 if success (thread terminated)
	//	return: 0 if fail
	//	return: -1 if invalid params (null pointers)
	int a3threadTerminate(a3_Thread *thread);

	// A3: Wait indefinitely for a running thread to finish.
	//	param thread: container of thread to wait for
	//	return: 1 if success (thread ended)
	//	return: 0 if fail
	//	return: -1 if invalid params (null pointers)
	int a3threadWait(a3_Thread *thread);

	// A3: Check if a thread is still running.
	//	param thread: container of thread to check
	//	return: thread status (1 for running, 0 for not running)
	//	return: -1 if invalid params (null pointers)
	int a3threadIsRunning(const a3_Thread *thread);


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// __ANIMAL3D_THREAD_H