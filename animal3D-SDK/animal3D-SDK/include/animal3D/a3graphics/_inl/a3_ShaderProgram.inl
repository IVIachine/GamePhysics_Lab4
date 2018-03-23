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
	
	a3_ShaderProgram.inl
	GLSL shader and program inline definitions.

	**DO NOT MODIFY THIS FILE**
*/

#ifdef __ANIMAL3D_SHADERPROGRAM_H
#ifndef __ANIMAL3D_SHADERPROGRAM_INL
#define __ANIMAL3D_SHADERPROGRAM_INL


//-----------------------------------------------------------------------------

inline int a3shaderCreateFromSource(a3_Shader *shader_out, const a3_ShaderType type, const char *source)
{
	const char *sourceList[] = { source };
	return a3shaderCreateFromSourceList(shader_out, type, sourceList, 1);
}

inline int a3shaderCreateFromFile(a3_Shader *shader_out, const a3_ShaderType type, const char *filePath)
{
	const char *filePathList[] = { filePath };
	return a3shaderCreateFromFileList(shader_out, type, filePathList, 1);
}

inline int a3shaderReference(a3_Shader *shader)
{
	if (shader)
		return a3handleIncrementCount(shader->handle);
	return -1;
}

inline int a3shaderRelease(a3_Shader *shader)
{
	if (shader)
	{
		int ret = a3handleDecrementCount(shader->handle);
		if (ret == 0)
			shader->compiled = shader->attached = 0;
		return ret;
	}
	return -1;
}


//-----------------------------------------------------------------------------

inline int a3shaderProgramReference(a3_ShaderProgram *program)
{
	if (program)
		return a3handleIncrementCount(program->handle);
	return -1;
}

inline int a3shaderProgramRelease(a3_ShaderProgram *program)
{
	if (program)
	{
		int ret = a3handleDecrementCount(program->handle);
		if (ret == 0)
			program->linked = program->validated = 0;
		return ret;
	}
	return -1;
}


//-----------------------------------------------------------------------------


#endif	// !__ANIMAL3D_SHADERPROGRAM_INL
#endif	// __ANIMAL3D_SHADERPROGRAM_H