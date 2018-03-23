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
	
	a3_ShaderProgram.h
	GLSL shader and program creation and management.

	**DO NOT MODIFY THIS FILE**
*/

#ifndef __ANIMAL3D_SHADERPROGRAM_H
#define __ANIMAL3D_SHADERPROGRAM_H


#include "a3_GraphicsObjectHandle.h"


#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_Shader				a3_Shader;
	typedef struct a3_ShaderProgram			a3_ShaderProgram;
	typedef enum a3_ShaderType				a3_ShaderType;
	typedef enum a3_UniformType				a3_UniformType;
	typedef enum a3_UniformMatType			a3_UniformMatType;
#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// A3: Enumeration of shader types.
	enum a3_ShaderType
	{
		a3shader_vertex,
		a3shader_tessellationControl,
		a3shader_tessellationEvaluation,
		a3shader_geometry,
		a3shader_fragment,
		a3shader_compute,
	};


	// A3: GLSL shader object descriptor.
	//	member handle: graphics handle container
	//	member type: which type of shader is this
	//	member compiled: compile status of shader (should be 1 after creation)
	//	member attached: counter for the number of programs using this shader 
	//		(can be used more than once; in fact, it's better to load a shared 
	//		shader and have it be used multiple times instead of loading it 
	//		one time per each program that uses it)
	struct a3_Shader
	{
		a3_GraphicsObjectHandle handle[1];
		a3_ShaderType type;
		short compiled, attached;
	};


//-----------------------------------------------------------------------------
	
	// A3: GLSL shader program object descriptor.
	//	member handle: graphics handle container
	//	member linked: link status of program (can modify if not linked)
	//	member validated: validation status of program (optional)
	//	member shadersAttached: handles of each shader type attached to the 
	//		program; only one shader of each type can be attached; compute 
	//		shaders must have a program of their own
	struct a3_ShaderProgram
	{
		a3_GraphicsObjectHandle handle[1];
		short linked, validated;
		short shadersAttached[6];
	};

	
//-----------------------------------------------------------------------------

	// A3: Multi-purpose uniform types.
	enum a3_UniformType
	{
		a3unif_single,
		a3unif_vec2,
		a3unif_vec3,
		a3unif_vec4,
	};

	// A3: Multi-purpose uniform matrix types.
	enum a3_UniformMatType
	{
		a3unif_mat2,
		a3unif_mat3,
		a3unif_mat4,
		a3unif_mat2x3,
		a3unif_mat2x4,
		a3unif_mat3x4,
	};


//-----------------------------------------------------------------------------
	// shader functions

	// A3: Create GLSL shader from raw source.
	//	param shader_out: non-null pointer to uninitialized shader descriptor
	//	param type: shader type
	//	param source: non-null cstring of GLSL shader source
	//	return: 1 if success
	//	return: 0 if compile failed or invalid string
	//	return: -1 if invalid params or shader is already initialized
	int a3shaderCreateFromSource(a3_Shader *shader_out, const a3_ShaderType type, const char *source);

	// A3: Create GLSL shader from multiple raw sources.
	//	param shader_out: non-null pointer to uninitialized shader descriptor
	//	param type: shader type
	//	param sourceList: non-null array of cstrings of GLSL shader sources; 
	//		null or empty strings will be ignored
	//	param count: non-zero source string count
	//	return: number of sources compiled into shader if success
	//	return: 0 if compile failed or no valid strings
	//	return: -1 if invalid params or shader is already initialized
	int a3shaderCreateFromSourceList(a3_Shader *shader_out, const a3_ShaderType type, const char **sourceList, const unsigned int count);

	// A3: Create GLSL shader from text file.
	//	param shader_out: non-null pointer to uninitialized shader descriptor
	//	param type: shader type
	//	param filePath: non-null cstring of GLSL shader file path
	//	return: 1 if success
	//	return: 0 if compile failed or invalid string
	//	return: -1 if invalid params or shader is already initialized
	int a3shaderCreateFromFile(a3_Shader *shader_out, const a3_ShaderType type, const char *filePath);

	// A3: Create GLSL shader from multiple text files.
	//	param shader_out: non-null pointer to uninitialized shader descriptor
	//	param type: shader type
	//	param filePathList: non-null array of cstrings of shader file paths; 
	//		null or empty strings will be ignored
	//	param count: non-zero file path string count
	//	return: number of sources compiled into shader if success
	//	return: 0 if compile failed or no valid strings
	//	return: -1 if invalid params or shader is already initialized
	int a3shaderCreateFromFileList(a3_Shader *shader_out, const a3_ShaderType type, const char **filePathList, const unsigned int count);

	// A3: Reference shader.
	//	param shader: non-null pointer to shader descriptor to be referenced
	//	return: reference count if success
	//	return: -1 if invalid param or internal handle is zero (uninitialized)
	int a3shaderReference(a3_Shader *shader);

	// A3: Release shader.
	//	param shader: non-null pointer to shader descriptor to be released
	//	return: positive if removed reference
	//	return: 0 if object was permanently released; descriptor is reset
	//	return: -1 if invalid param or internal handle is zero (uninitialized)
	int a3shaderRelease(a3_Shader *shader);

	// A3: Update shader release function (hotload quick-fix).
	//	param shader: non-null pointer to shader to update
	//	return: 1 if success
	//	return: -1 if invalid params
	int a3shaderHandleUpdateReleaseCallback(a3_Shader *shader);


//-----------------------------------------------------------------------------
	// program functions

	// A3: Create GLSL program object.
	//	param program_out: non-null pointer to uninitialized program descriptor
	//	return: 1 if success
	//	return: 0 if creation failed
	//	return: -1 if invalid param or program is already initialized
	int a3shaderProgramCreate(a3_ShaderProgram *program_out);

	// A3: Attach shader object to a program object before linking program.
	//	param program: non-null pointer to program descriptor
	//	param shader: non-null pointer to shader descriptor
	//	return: 1 if success
	//	return: 0 if a shader of that type is already attached or cannot attach
	//	return: -1 if invalid params or internal handles are zero
	int a3shaderProgramAttachShader(a3_ShaderProgram *program, a3_Shader *shader);

	// A3: Detach shader object from program object before linking program.
	//	param program: non-null pointer to program descriptor
	//	param type: type of shader to be detached
	//	return: 1 if success
	//	return: 0 if a shader of the specified type is not attached
	//	return: -1 if invalid params or internal handle is zero
	int a3shaderProgramDetachShaderType(a3_ShaderProgram *program, const a3_ShaderType type);

	// A3: Link program; required before using program.
	//	param program: non-null pointer to program descriptor
	//	return: 1 if success; program object cannot be changed
	//	return: 0 if failed or already linked
	//	return: -1 if invalid param or internal handle is zero
	int a3shaderProgramLink(a3_ShaderProgram *program);

	// A3: Validate program; optional step after linking. This is used to 
	//		identify hardware-specific problems that may occur with the 
	//		attached shader objects.
	//	param program: non-null pointer to program descriptor
	//	return: 1 if success
	//	return: 0 if failed or already validated
	//	return: -1 if invalid param or internal handle is zero
	int a3shaderProgramValidate(a3_ShaderProgram *program);

	// A3: Activate/deactivate program; primitives drawn after calling this 
	//		function will flow through the pipeline defined by this program.
	//	param program: pointer to program descriptor; pass uninitialized 
	//		program or null pointer to deactivate the current program (i.e. 
	//		revert to fixed-function pipeline)
	//	return: 1 if program was activated (using programmable pipeline)
	//	return: 0 if active program was deactivated (using fixed-function)
	int a3shaderProgramActivate(const a3_ShaderProgram *program);

	// A3: Deactivate any active program (i.e. use fixed-function pipeline).
	//	return: 0
	int a3shaderProgramDeactivate();

	// A3: Reference program.
	//	param program: non-null pointer to program descriptor to be released
	//	return: new reference count if success
	//	return: -1 if invalid param or internal handle is zero (uninitialized)
	int a3shaderProgramReference(a3_ShaderProgram *program);

	// A3: Release program.
	//	param program: non-null pointer to program descriptor to be released
	//	return: positive if removed reference
	//	return: 0 if object was permanently released; descriptor is reset
	//	return: -1 if invalid param or internal handle is zero (uninitialized)
	int a3shaderProgramRelease(a3_ShaderProgram *program);

	// A3: Update shader program release function (hotload quick-fix).
	//	param program: non-null pointer to shader program to update
	//	return: 1 if success
	//	return: -1 if invalid params
	int a3shaderProgramHandleUpdateReleaseCallback(a3_ShaderProgram *program);


//-----------------------------------------------------------------------------
	// uniform functions

	// A3: Get location of a uniform in a program.
	//	param program: non-null pointer to program descriptor
	//	param uniformName: non-null, non-empty cstring of uniform variable name
	//	return: location of uniform in program (non-negative) if success
	//	return: -1 if invalid params or uniform was not found
	int a3shaderUniformGetLocation(const a3_ShaderProgram *program, const char *uniformName);

	// A3: Send integer or integer vector uniforms to active program.
	//	param uniformType: integer type of uniform (see enums above)
	//	param uniformLocation: non-negative uniform location in active program
	//	param count: non-zero quantity of this uniform type being sent
	//	param values: non-null array of integer-type values to be sent
	//	return: uniforms sent
	//	return: -1 if invalid params
	int a3shaderUniformSendInt(const a3_UniformType uniformType, const int uniformLocation, const unsigned int count, const int *values);

	// A3: Send float or float vector uniforms to active program.
	//	param uniformType: float type of uniform (see enums above)
	//	param uniformLocation: non-negative uniform location in active program
	//	param count: non-zero quantity of this uniform type being sent
	//	param values: non-null array of float-type values to be sent
	//	return: uniforms sent
	//	return: -1 if invalid params
	int a3shaderUniformSendFloat(const a3_UniformType uniformType, const int uniformLocation, const unsigned int count, const float *values);

	// A3: Send double or double vector uniforms to active program.
	//	param uniformType: integer type of uniform (see enums above)
	//	param uniformLocation: non-negative uniform location in active program
	//	param count: non-zero quantity of this uniform type being sent
	//	param values: non-null array of double-type values to be sent
	//	return: uniforms sent
	//	return: -1 if invalid params
	int a3shaderUniformSendDouble(const a3_UniformType uniformType, const int uniformLocation, const unsigned int count, const double *values);

	// A3: Send float matrix uniforms to active program.
	//	param uniformType: float matrix type of uniform (see enums above)
	//	param transpose: whether or not the data should be transposed
	//	param uniformLocation: non-negative uniform location in active program
	//	param count: non-zero quantity of this uniform type being sent
	//	param values: non-null array of float matrix-type values to be sent
	//	return: uniforms sent
	//	return: -1 if invalid params
	int a3shaderUniformSendFloatMat(const a3_UniformMatType uniformType, const int transpose, const int uniformLocation, const unsigned int count, const float *values);

	// A3: Send double matrix uniforms to active program.
	//	param uniformType: double matrix type of uniform (see enums above)
	//	param transpose: whether or not the data should be transposed
	//	param uniformLocation: non-negative uniform location in active program
	//	param count: non-zero quantity of this uniform type being sent
	//	param values: non-null array of double matrix-type values to be sent
	//	return: uniforms sent
	//	return: -1 if invalid params
	int a3shaderUniformSendDoubleMat(const a3_UniformMatType uniformType, const int transpose, const int uniformLocation, const unsigned int count, const double *values);

	// A3: Get maximum number of uniform locations.
	//	return: max uniform locations
	int a3shaderUniformMaxLocations();


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#include "_inl/a3_ShaderProgram.inl"


#endif	// !__ANIMAL3D_SHADERPROGRAM_H