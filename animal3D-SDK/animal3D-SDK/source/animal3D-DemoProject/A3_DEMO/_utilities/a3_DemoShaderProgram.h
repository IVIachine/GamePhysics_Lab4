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
	
	a3_DemoShaderProgram.h
	Wrapper for shaders used in demo state.
*/

#ifndef __ANIMAL3D_DEMOSHADERPROGRAM_H
#define __ANIMAL3D_DEMOSHADERPROGRAM_H


//-----------------------------------------------------------------------------
// animal3D framework includes

#include "animal3D/a3graphics/a3_ShaderProgram.h"


//-----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_DemoStateShaderProgram	a3_DemoStateShaderProgram;
#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// maximum number of uniforms in a program
	enum a3_DemoShaderProgramMaxCounts
	{
		demoStateMaxCount_shaderProgramUniform = 16,
	};


	// structure to help with shader program and uniform management
	struct a3_DemoStateShaderProgram
	{
		a3_ShaderProgram program[1];
		union {
			int uniformLocation[demoStateMaxCount_shaderProgramUniform];
			struct {
				int
					// common vertex shader uniforms
					uMVP,						// model-view-projection transform
					uVP,						// view-projection transform (camera matrix)
					uM,							// model matrix
					uMV,						// model-view matrix
					uP,							// projection matrix
					uNrm,						// normal matrix (inverse transpose model-view or model to cancel non-uniform scale)

					// common geometry shader uniforms
					uRayOrigin_world,
					uRayDirection_world,
					uRayLength,

					// common fragment shader uniforms
					uLightPos_eye,				// light position in eye space
					uTex_dm,					// diffuse map/surface texture
					uColor;						// uniform color
			};
		};
	};


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_DEMOSHADERPROGRAM_H