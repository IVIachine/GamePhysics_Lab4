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
	
	a3_DemoState.h
	Demo state interface and programmer function declarations.

	********************************************
	*** THIS IS YOUR DEMO'S MAIN HEADER FILE ***
	********************************************
*/

#ifndef __ANIMAL3D_DEMOSTATE_H
#define __ANIMAL3D_DEMOSTATE_H


//-----------------------------------------------------------------------------
// animal3D framework includes

#include "animal3D/animal3D.h"


//-----------------------------------------------------------------------------
// other demo includes

#include "_utilities/a3_DemoSceneObject.h"
#include "_utilities/a3_DemoShaderProgram.h"


//-----------------------------------------------------------------------------
// physics includes

#include "physics/a3_PhysicsWorld.h"


//-----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_DemoState					a3_DemoState;
#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// object maximum counts for easy array storage
	// good idea to make these numbers greater than what you actually need 
	//	and if you end up needing more just increase the count... there's 
	//	more than enough memory to hold extra objects
	enum a3_DemoStateObjectMaxCounts
	{
		demoStateMaxCount_sceneObject = 32,
		demoStateMaxCount_camera = 1,
		demoStateMaxCount_timer = 1,
		demoStateMaxCount_thread = 1,
		demoStateMaxCount_texture = 2,
		demoStateMaxCount_drawDataBuffer = 1,
		demoStateMaxCount_vertexArray = 4,
		demoStateMaxCount_drawable = 16,
		demoStateMaxCount_shaderProgram = 8,
	};


//-----------------------------------------------------------------------------

	// persistent demo state data structure
	struct a3_DemoState
	{
		//---------------------------------------------------------------------
		// general variables pertinent to the state

		// terminate key pressed
		int exitFlag;

		// global vertical axis: Z = 0, Y = 1
		int verticalAxis;

		// asset streaming between loads enabled (careful!)
		int streaming;

		// window and full-frame dimensions
		unsigned int windowWidth, windowHeight;
		unsigned int frameWidth, frameHeight;
		int frameBorder;


		//---------------------------------------------------------------------
		// objects that have known or fixed instance count in the whole demo

		// text renderer
		int textInit, showText;
		a3_TextRenderer text[1];

		// input
		a3_MouseInput mouse[1];
		a3_KeyboardInput keyboard[1];
		a3_XboxControllerInput xcontrol[4];

		// pointer to fast trig table
		float trigTable[4096 * 4];


		//---------------------------------------------------------------------
		// scene variables and objects

		unsigned int demoMode, demoModeCount;


		// physics world embedded in demo
		a3_PhysicsWorld physicsWorld[1];

		// object counts
		unsigned int rigidbodyObjects, particleObjects;


		// extra display options
		int displayGrid, displayAxes, displaySkybox;
		int displayPhysicsText;


		// ray and prevailing hit
		a3_Ray ray[1];
		a3_RayHit rayHit[1];
		int hitIndex;

		// collision tracking
		int colliding[physicsMaxCount_rigidbody];


		// dummy drawable for rays
		a3_VertexDrawable dummyDrawable[1];

		// pointer to drawable attached to each physics object
		const a3_VertexDrawable *rbDrawable[physicsMaxCount_rigidbody];


		//---------------------------------------------------------------------
		// object arrays: organized as anonymous unions for two reasons: 
		//	1. easy to manage entire sets of the same type of object using the 
		//		array component
		//	2. at the same time, variables are named pointers

		// scene objects
		union {
			a3_DemoSceneObject sceneObject[demoStateMaxCount_sceneObject];
			struct {
				a3_DemoSceneObject
					cameraObject[1];					// transform for camera
			};
			struct {
				a3_DemoSceneObject
					graphicsObjects[1];					// strictly graphics
				a3_DemoSceneObject
					physicsRigidbodies[20], 
					physicsParticles[1];				// physics aligned rigid bodies & particles
			};
			struct {
				a3_DemoSceneObject graphicsCam[1];

				// rigid body objects
				a3_DemoSceneObject groundObject[1];
				a3_DemoSceneObject sphereObject[4];
				a3_DemoSceneObject cylinderObject[4];
				a3_DemoSceneObject boxObject[4];

				// no named particles here
			};
		};

		// cameras
		union {
			a3_DemoCamera camera[demoStateMaxCount_camera];
			struct {
				a3_DemoCamera
					sceneCamera[1];						// scene viewing camera
			};
		};

		// timers
		union {
			a3_Timer timer[demoStateMaxCount_timer];
			struct {
				a3_Timer
					renderTimer[1];						// render FPS timer
			};
		};

		// threads
		union {
			a3_Thread thread[demoStateMaxCount_thread];
			struct {
				a3_Thread
					physicsThread[1];					// independent physics simulation thread
			};
		};


		// textures
		union {
			a3_Texture texture[demoStateMaxCount_texture];
			struct {
				a3_Texture
					tex_skybox_clouds[1],				// clouds skybox
					tex_checker[1];						// checkered texture
			};
		};


		// draw data buffers
		union {
			a3_VertexBuffer drawDataBuffer[demoStateMaxCount_drawDataBuffer];
			struct {
				a3_VertexBuffer
					vbo_staticSceneObjectDrawBuffer[1];			// buffer to hold all data for static scene objects (e.g. grid)
			};
		};

		// vertex array objects
		union {
			a3_VertexArrayDescriptor vertexArray[demoStateMaxCount_vertexArray];
			struct {
				a3_VertexArrayDescriptor
					vao_position[1],							// VAO for vertex format with only position
					vao_position_color[1],						// VAO for vertex format with position and color
					vao_position_texcoord[1],					// VAO for vertex format with position and UVs
					vao_tangent_basis[1];						// VAO for vertex format with full tangent basis
			};
		};

		// drawables
		union {
			a3_VertexDrawable drawable[demoStateMaxCount_drawable];
			struct {
				a3_VertexDrawable
					draw_grid[1],								// wireframe ground plane to emphasize scaling
					draw_axes[1],								// coordinate axes at the center of the world
					draw_skybox[1],								// skybox cube mesh
					draw_arrow[1],								// unit arrow mesh (1 unit, small base)
					draw_plane[1],								// unit plane mesh (1x1)
					draw_box[1],								// unit box mesh (1x1x1)
					draw_sphere[1],								// unit sphere mesh (r=1)
					draw_cylinder[1],							// unit cylinder mesh (r=1,h=1)
					draw_torus[1],								// unit torus mesh (R=1,r=1... basically a very plump donut)
					draw_teapot[1];								// can't not have a Utah teapot
			};
		};


		// shader programs and uniforms
		union {
			a3_DemoStateShaderProgram shaderProgram[demoStateMaxCount_shaderProgram];
			struct {
				a3_DemoStateShaderProgram
					prog_drawRay[1],					// draw a ray
					prog_drawLambert[1],				// draw Lambert shading
					prog_drawColorUnifTexture[1],		// draw texture mixed with solid color
					prog_drawColor[1],					// draw color attribute
					prog_drawColorUnif[1];				// draw uniform color
			};
		};


		//---------------------------------------------------------------------
	};

	
//-----------------------------------------------------------------------------

	// demo-related functions

	// other utilities
	void a3demo_setDefaultGraphicsState();

	// loading and unloading
	void a3demo_loadTextures(a3_DemoState *demoState);
	void a3demo_loadGeometry(a3_DemoState *demoState);
	void a3demo_loadShaders(a3_DemoState *demoState);

	void a3demo_unloadTextures(a3_DemoState *demoState);
	void a3demo_unloadGeometry(a3_DemoState *demoState);
	void a3demo_unloadShaders(a3_DemoState *demoState);

	void a3demo_initScene(a3_DemoState *demoState);

	void a3demo_refresh(a3_DemoState *demoState);

	void a3demo_validateUnload(const a3_DemoState *demoState);

	// main loop
	void a3demo_input(a3_DemoState *demoState, double dt);
	void a3demo_update(a3_DemoState *demoState, double dt);
	void a3demo_render(const a3_DemoState *demoState);


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_DEMOSTATE_H