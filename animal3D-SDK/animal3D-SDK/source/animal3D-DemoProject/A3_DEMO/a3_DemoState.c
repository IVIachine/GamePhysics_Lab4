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

	a3_DemoState.c/.cpp
	Demo state function implementations.

	********************************************
	*** THIS IS YOUR DEMO'S MAIN SOURCE FILE ***
	*** Implement your demo logic here.      ***
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

//Tyler merged raycast and hull collisions code with help from Charlie

#include "a3_DemoState.h"


//-----------------------------------------------------------------------------

// OpenGL
#ifdef _WIN32
#include <Windows.h>
#include <GL/GL.h>
#else	// !_WIN32
#include <OpenGL/gl3.h>
#endif	// _WIN32


#include <stdio.h>
#include <stdlib.h>


//-----------------------------------------------------------------------------
// GENERAL UTILITIES

inline void a3demo_applyScale_internal(a3_DemoSceneObject *sceneObject, a3real4x4p s, const int overrideInverse)
{
	// override inverse without scale first
	if (overrideInverse)
		a3real4x4TransformInverseIgnoreScale(sceneObject->modelMatInv.m, sceneObject->modelMat.m);

	if (sceneObject->scaleMode)
	{
		if (sceneObject->scaleMode == 1)
		{
			s[0][0] = s[1][1] = s[2][2] = sceneObject->scale.x;
		}
		else
		{
			s[0][0] = sceneObject->scale.x;
			s[1][1] = sceneObject->scale.y;
			s[2][2] = sceneObject->scale.z;
		}
		a3real4x4ConcatL(sceneObject->modelMat.m, s);
	}
}


//-----------------------------------------------------------------------------
// SETUP AND TERMINATION UTILITIES

// set default GL state
void a3demo_setDefaultGraphicsState()
{
	const float lineWidth = 2.0f;
	const float pointSize = 4.0f;

	// lines and points
	glLineWidth(lineWidth);
	glPointSize(pointSize);

	// backface culling
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	// alpha blending
	// result = ( new*[new alpha] ) + ( old*[1 - new alpha] )
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// textures
	glEnable(GL_TEXTURE_2D);

	// background
	glClearColor(0.0f, 0.0f, 0.0, 0.0f);
}


//-----------------------------------------------------------------------------
// LOADING AND UNLOADING

// utility to load textures
void a3demo_loadTextures(a3_DemoState *demoState)
{
	// pointer to texture
	a3_Texture *tex;
	unsigned int i;

	// list of texture files to load
	const char *texFiles[] = {
		"../../../../resource/tex/bg/sky_clouds.png",
		"../../../../resource/tex/sprite/checker.png",
	};
	const unsigned int numTextures = sizeof(texFiles) / sizeof(const char *);

	for (i = 0; i < numTextures; ++i)
	{
		tex = demoState->texture + i;
		a3textureCreateFromFile(tex, texFiles[i]);
	}

	// change settings on a per-texture basis
	a3textureActivate(demoState->tex_skybox_clouds, a3tex_unit00);
	a3textureDefaultSettings();	// nearest filtering, repeat on both axes

	a3textureActivate(demoState->tex_checker, a3tex_unit00);
	a3textureDefaultSettings();	// nearest filtering, repeat on both axes


	// done
	a3textureDeactivate(a3tex_unit00);
}

// utility to load geometry
void a3demo_loadGeometry(a3_DemoState *demoState)
{
	// static model transformations
	static const a3mat4 downscale20x = {
		+0.05f, 0.0f, 0.0f, 0.0f,
		0.0f, +0.05f, 0.0f, 0.0f,
		0.0f, 0.0f, +0.05f, 0.0f,
		0.0f, 0.0f, 0.0f, +1.0f,
	};

	// pointer to shared vbo/ibo
	a3_VertexBuffer *vbo_ibo;
	a3_VertexArrayDescriptor *vao;
	a3_VertexDrawable *currentDrawable;
	unsigned int sharedVertexStorage = 0, sharedIndexStorage = 0;
	unsigned int numVerts = 0;
	unsigned int i;


	// file streaming (if requested)
	a3_FileStream fileStream[1] = { 0 };
	const char *const geometryStream = "./data/geom_data.dat";

	// geometry data
	a3_GeometryData sceneShapesData[3] = { 0 };
	a3_GeometryData proceduralShapesData[6] = { 0 };
	a3_GeometryData loadedModelsData[1] = { 0 };
	const unsigned int sceneShapesCount = sizeof(sceneShapesData) / sizeof(a3_GeometryData);
	const unsigned int proceduralShapesCount = sizeof(proceduralShapesData) / sizeof(a3_GeometryData);
	const unsigned int loadedModelsCount = sizeof(loadedModelsData) / sizeof(a3_GeometryData);

	// common index format
	a3_IndexFormatDescriptor sceneCommonIndexFormat[1] = { 0 };


	// procedural scene objects
	// attempt to load stream if requested
	if (demoState->streaming && a3fileStreamOpenRead(fileStream, geometryStream))
	{
		// read from stream

		// static scene objects
		for (i = 0; i < sceneShapesCount; ++i)
			a3fileStreamReadObject(fileStream, sceneShapesData + i, (a3_FileStreamReadFunc)a3geometryLoadDataBinary);

		// procedurally-generated objects
		for (i = 0; i < proceduralShapesCount; ++i)
			a3fileStreamReadObject(fileStream, proceduralShapesData + i, (a3_FileStreamReadFunc)a3geometryLoadDataBinary);

		// loaded model objects
		for (i = 0; i < loadedModelsCount; ++i)
			a3fileStreamReadObject(fileStream, loadedModelsData + i, (a3_FileStreamReadFunc)a3geometryLoadDataBinary);

		// done
		a3fileStreamClose(fileStream);
	}
	// not streaming or stream doesn't exist
	else if (!demoState->streaming || a3fileStreamOpenWrite(fileStream, geometryStream))
	{
		// create new data
		a3_ProceduralGeometryDescriptor sceneShapes[3] = { a3geomShape_none };
		a3_ProceduralGeometryDescriptor proceduralShapes[6] = { a3geomShape_none };
		a3_ProceduralGeometryDescriptor loadedModelShapes[1] = { a3geomShape_none };

		// static scene procedural objects
		//	(axes, grid, skybox)
		a3proceduralCreateDescriptorAxes(sceneShapes + 0, a3geomFlag_wireframe, 0.0f, 1);
		a3proceduralCreateDescriptorPlane(sceneShapes + 1, a3geomFlag_wireframe, a3geomAxis_default, 32.0f, 32.0f, 32, 32);
		a3proceduralCreateDescriptorBox(sceneShapes + 2, a3geomFlag_texcoords, 100.0f, 100.0f, 100.0f, 1, 1, 1);
		for (i = 0; i < sceneShapesCount; ++i)
		{
			a3proceduralGenerateGeometryData(sceneShapesData + i, sceneShapes + i);
			a3fileStreamWriteObject(fileStream, sceneShapesData + i, (a3_FileStreamWriteFunc)a3geometrySaveDataBinary);
		}

		// other procedurally-generated objects
		a3proceduralCreateDescriptorPyramid(proceduralShapes + 0, a3geomFlag_vanilla, a3geomAxis_default, 0.1f, 1.0f);
		a3proceduralCreateDescriptorPlane(proceduralShapes + 1, a3geomFlag_tangents, a3geomAxis_default, 1.0f, 1.0f, 1, 1);
		a3proceduralCreateDescriptorBox(proceduralShapes + 2, a3geomFlag_tangents, 1.0f, 1.0f, 1.0f, 1, 1, 1);
		a3proceduralCreateDescriptorSphere(proceduralShapes + 3, a3geomFlag_tangents, a3geomAxis_default, 1.0f, 32, 24);
		a3proceduralCreateDescriptorCylinder(proceduralShapes + 4, a3geomFlag_tangents, a3geomAxis_default, 1.0f, 1.0f, 32, 1, 1);
		a3proceduralCreateDescriptorTorus(proceduralShapes + 5, a3geomFlag_tangents, a3geomAxis_default, 1.0f, 1.0f, 32, 24);
		for (i = 0; i < proceduralShapesCount; ++i)
		{
			a3proceduralGenerateGeometryData(proceduralShapesData + i, proceduralShapes + i);
			a3fileStreamWriteObject(fileStream, proceduralShapesData + i, (a3_FileStreamWriteFunc)a3geometrySaveDataBinary);
		}

		// objects loaded from mesh files
		a3modelLoadOBJ(loadedModelsData + 0, "../../../../resource/obj/teapot/teapot.obj", a3model_calculateVertexTangents, downscale20x.mm);
		for (i = 0; i < loadedModelsCount; ++i)
			a3fileStreamWriteObject(fileStream, loadedModelsData + i, (a3_FileStreamWriteFunc)a3geometrySaveDataBinary);

		// done
		a3fileStreamClose(fileStream);
	}


	// GPU data upload process: 
	//	- determine storage requirements
	//	- allocate buffer
	//	- create vertex arrays using unique formats
	//	- create drawable and upload data

	// get storage size
	sharedVertexStorage = numVerts = 0;
	for (i = 0; i < sceneShapesCount; ++i)
	{
		sharedVertexStorage += a3geometryGetVertexBufferSize(sceneShapesData + i);
		numVerts += sceneShapesData[i].numVertices;
	}
	for (i = 0; i < proceduralShapesCount; ++i)
	{
		sharedVertexStorage += a3geometryGetVertexBufferSize(proceduralShapesData + i);
		numVerts += proceduralShapesData[i].numVertices;
	}
	for (i = 0; i < loadedModelsCount; ++i)
	{
		sharedVertexStorage += a3geometryGetVertexBufferSize(loadedModelsData + i);
		numVerts += loadedModelsData[i].numVertices;
	}

	// common index format required for shapes that share vertex formats
	a3geometryCreateIndexFormat(sceneCommonIndexFormat, numVerts);
	sharedIndexStorage = 0;
	for (i = 0; i < sceneShapesCount; ++i)
		sharedIndexStorage += a3indexStorageSpaceRequired(sceneCommonIndexFormat, sceneShapesData[i].numIndices);
	for (i = 0; i < proceduralShapesCount; ++i)
		sharedIndexStorage += a3indexStorageSpaceRequired(sceneCommonIndexFormat, proceduralShapesData[i].numIndices);
	for (i = 0; i < loadedModelsCount; ++i)
		sharedIndexStorage += a3indexStorageSpaceRequired(sceneCommonIndexFormat, loadedModelsData[i].numIndices);


	// create shared buffer
	vbo_ibo = demoState->vbo_staticSceneObjectDrawBuffer;
	a3bufferCreateSplit(vbo_ibo, a3buffer_vertex, sharedVertexStorage, sharedIndexStorage, 0, 0);
	sharedVertexStorage = 0;


	// create vertex formats and drawables
	// axes
	vao = demoState->vao_position_color;
	a3geometryGenerateVertexArray(vao, sceneShapesData + 0, vbo_ibo, sharedVertexStorage);
	currentDrawable = demoState->draw_axes;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, sceneShapesData + 0, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);

	// grid: position attribute only
	// overlay objects are also just position
	vao = demoState->vao_position;
	a3geometryGenerateVertexArray(vao, sceneShapesData + 1, vbo_ibo, sharedVertexStorage);
	currentDrawable = demoState->draw_grid;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, sceneShapesData + 1, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);

	// skybox: position and texture coordinates
	vao = demoState->vao_position_texcoord;
	a3geometryGenerateVertexArray(vao, sceneShapesData + 2, vbo_ibo, sharedVertexStorage);
	currentDrawable = demoState->draw_skybox;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, sceneShapesData + 2, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);

	// arrowhead thingy: just position
	vao = demoState->vao_position;
	currentDrawable = demoState->draw_arrow;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, proceduralShapesData + 0, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);

	// other scene objects: full tangent basis (in case we want them to be fancy)
	vao = demoState->vao_tangent_basis;
	a3geometryGenerateVertexArray(vao, loadedModelsData + 0, vbo_ibo, sharedVertexStorage);
	currentDrawable = demoState->draw_plane;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, proceduralShapesData + 1, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);
	currentDrawable = demoState->draw_box;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, proceduralShapesData + 2, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);
	currentDrawable = demoState->draw_sphere;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, proceduralShapesData + 3, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);
	currentDrawable = demoState->draw_cylinder;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, proceduralShapesData + 4, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);
	currentDrawable = demoState->draw_torus;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, proceduralShapesData + 5, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);

	// use fancy format for teapot
	currentDrawable = demoState->draw_teapot;
	sharedVertexStorage += a3geometryGenerateDrawable(currentDrawable, loadedModelsData + 0, vao, vbo_ibo, sceneCommonIndexFormat, 0, 0);

	// release data when done
	for (i = 0; i < sceneShapesCount; ++i)
		a3geometryReleaseData(sceneShapesData + i);
	for (i = 0; i < proceduralShapesCount; ++i)
		a3geometryReleaseData(proceduralShapesData + i);
	for (i = 0; i < loadedModelsCount; ++i)
		a3geometryReleaseData(loadedModelsData + i);


	// DUMMY DRAWABLE: DO NOT RELEASE
	// single point duplicate of any other mesh with 
	//	the sole purpose of triggering vertex shader
	*demoState->dummyDrawable = *demoState->draw_axes;
	demoState->dummyDrawable->count = 1;
	demoState->dummyDrawable->primitive = GL_POINTS;
}


// utility to load shaders
void a3demo_loadShaders(a3_DemoState *demoState)
{
	// direct to demo programs
	a3_DemoStateShaderProgram *currentDemoProg;
	int *currentUnif, uLocation;
	unsigned int i, j;

	// list of uniform names: align with uniform list in demo struct!
	const char *uniformNames[demoStateMaxCount_shaderProgramUniform] = {
		// common vertex
		"uMVP",
		"uVP",
		"uM",
		"uMV",
		"uP",
		"uNrm",

		// common geometry
		"uRayOrigin_world",
		"uRayDirection_world",
		"uRayLength",

		// common fragment
		"uLightPos_eye",
		"uTex_dm",
		"uColor",
	};


	// some default uniform values
	const float defaultColor[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
	const int defaultTexUnits[] = { 0, 1, 2, 3, 4, 5, 6, 7 };
	const int defaultInt[] = { 0 };
	const float defaultFloat[] = { 0.0f };


	// list of all unique shaders
	// this is a good idea to avoid multi-loading 
	//	those that are shared between programs
	union {
		struct {
			// vertex shaders
			a3_Shader passLambertComponents_eye_transform_vs[1];
			a3_Shader passTexcoord_transform_vs[1];
			a3_Shader passthru_transform_vs[1];
			a3_Shader passColor_transform_vs[1];
			a3_Shader dummy_vs[1];

			// geometry shaders
			a3_Shader drawRay_gs[1];

			// fragment shaders
			a3_Shader drawLambert_eye_fs[1];
			a3_Shader drawColorUnifTexture_fs[1];
			a3_Shader drawColorUnif_fs[1];
			a3_Shader drawColorAttrib_fs[1];
		};
	} shaderList = { 0 };
	a3_Shader *const shaderListPtr = (a3_Shader *)(&shaderList);
	const unsigned int numUniqueShaders = sizeof(shaderList) / sizeof(a3_Shader);

	// descriptors to help load shaders; aligned with above list
	struct {
		a3_ShaderType shaderType;
		unsigned int srcCount;
		const char *filePath[8];	// max number of source files per shader
	} shaderDescriptor[] = {
		{ a3shader_vertex,		1, { "../../../../resource/glsl/4x/vs/e/passLambertComponents_eye_transform_vs4x.glsl" } },
		{ a3shader_vertex,		1, { "../../../../resource/glsl/4x/vs/e/passTexcoord_transform_vs4x.glsl" } },
		{ a3shader_vertex,		1, { "../../../../resource/glsl/4x/vs/e/passthru_transform_vs4x.glsl" } },
		{ a3shader_vertex,		1, { "../../../../resource/glsl/4x/vs/e/passColor_transform_vs4x.glsl" } },
		{ a3shader_vertex,		1, { "../../../../resource/glsl/4x/vs/e/dummy_vs4x.glsl" } },

		{ a3shader_geometry,	1, { "../../../../resource/glsl/4x/gs/e/drawRay_gs4x.glsl" } },

		{ a3shader_fragment,	1, { "../../../../resource/glsl/4x/fs/e/drawLambert_eye_fs4x.glsl" } },
		{ a3shader_fragment,	1, { "../../../../resource/glsl/4x/fs/e/drawColorUnifTexture_fs4x.glsl" } },
		{ a3shader_fragment,	1, { "../../../../resource/glsl/4x/fs/e/drawColorUnif_fs4x.glsl" } },
		{ a3shader_fragment,	1, { "../../../../resource/glsl/4x/fs/e/drawColorAttrib_fs4x.glsl" } },
	};

	// load unique shaders: 
	//	- load file contents
	//	- create and compile shader object
	//	- release file contents
	for (i = 0; i < numUniqueShaders; ++i)
	{
		a3shaderCreateFromFileList(shaderListPtr + i, shaderDescriptor[i].shaderType,
			shaderDescriptor[i].filePath, shaderDescriptor[i].srcCount);
	}


	// setup programs: 
	//	- create program object
	//	- attach shader objects

	// draw ray program
	currentDemoProg = demoState->prog_drawRay;
	a3shaderProgramCreate(currentDemoProg->program);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.dummy_vs);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.drawRay_gs);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.drawColorUnif_fs);

	// Lambert shading program
	currentDemoProg = demoState->prog_drawLambert;
	a3shaderProgramCreate(currentDemoProg->program);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.passLambertComponents_eye_transform_vs);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.drawLambert_eye_fs);

	// color with texture program
	currentDemoProg = demoState->prog_drawColorUnifTexture;
	a3shaderProgramCreate(currentDemoProg->program);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.passTexcoord_transform_vs);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.drawColorUnifTexture_fs);

	// color attrib program
	currentDemoProg = demoState->prog_drawColor;
	a3shaderProgramCreate(currentDemoProg->program);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.passColor_transform_vs);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.drawColorAttrib_fs);

	// uniform color program
	currentDemoProg = demoState->prog_drawColorUnif;
	a3shaderProgramCreate(currentDemoProg->program);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.passthru_transform_vs);
	a3shaderProgramAttachShader(currentDemoProg->program, shaderList.drawColorUnif_fs);


	// activate a primitive for validation
	// makes sure the specified geometry can draw using programs
	// good idea to activate the drawable with the most attributes
	a3vertexActivateDrawable(demoState->draw_axes);

	// link and validate all programs
	for (i = 0; i < demoStateMaxCount_shaderProgram; ++i)
	{
		currentDemoProg = demoState->shaderProgram + i;
		a3shaderProgramLink(currentDemoProg->program);
		a3shaderProgramValidate(currentDemoProg->program);
	}

	// if linking fails, contingency plan goes here
	// otherwise, release shaders
	for (i = 0; i < numUniqueShaders; ++i)
		a3shaderRelease(shaderListPtr + i);


	// prepare uniforms algorithmically instead of manually for all programs
	for (i = 0; i < demoStateMaxCount_shaderProgram; ++i)
	{
		// activate program
		currentDemoProg = demoState->shaderProgram + i;
		a3shaderProgramActivate(currentDemoProg->program);

		// get uniform locations
		currentUnif = currentDemoProg->uniformLocation;
		for (j = 0; j < demoStateMaxCount_shaderProgramUniform; ++j)
			currentUnif[j] = a3shaderUniformGetLocation(currentDemoProg->program, uniformNames[j]);

		// set default values for all programs that have a uniform that will 
		//	either never change or is consistent for all programs
		if ((uLocation = currentDemoProg->uMVP) >= 0)
			a3shaderUniformSendFloatMat(a3unif_mat4, 0, uLocation, 1, a3identityMat4.mm);
		if ((uLocation = currentDemoProg->uVP) >= 0)
			a3shaderUniformSendFloatMat(a3unif_mat4, 0, uLocation, 1, a3identityMat4.mm);
		if ((uLocation = currentDemoProg->uM) >= 0)
			a3shaderUniformSendFloatMat(a3unif_mat4, 0, uLocation, 1, a3identityMat4.mm);
		if ((uLocation = currentDemoProg->uMV) >= 0)
			a3shaderUniformSendFloatMat(a3unif_mat4, 0, uLocation, 1, a3identityMat4.mm);
		if ((uLocation = currentDemoProg->uP) >= 0)
			a3shaderUniformSendFloatMat(a3unif_mat4, 0, uLocation, 1, a3identityMat4.mm);
		if ((uLocation = currentDemoProg->uNrm) >= 0)
			a3shaderUniformSendFloatMat(a3unif_mat4, 0, uLocation, 1, a3identityMat4.mm);

		if ((uLocation = currentDemoProg->uRayOrigin_world) >= 0)
			a3shaderUniformSendFloat(a3unif_vec4, uLocation, 1, a3wVec4.v);
		if ((uLocation = currentDemoProg->uRayDirection_world) >= 0)
			a3shaderUniformSendFloat(a3unif_vec4, uLocation, 1, a3zeroVec4.v);
		if ((uLocation = currentDemoProg->uRayLength) >= 0)
			a3shaderUniformSendFloat(a3unif_single, uLocation, 1, defaultFloat);

		if ((uLocation = currentDemoProg->uLightPos_eye) >= 0)
			a3shaderUniformSendFloat(a3unif_vec4, uLocation, 1, a3wVec4.v);
		if ((uLocation = currentDemoProg->uTex_dm) >= 0)
			a3shaderUniformSendInt(a3unif_single, uLocation, 1, defaultTexUnits + 0);
		if ((uLocation = currentDemoProg->uColor) >= 0)
			a3shaderUniformSendFloat(a3unif_vec4, uLocation, 1, defaultColor);
	}

	//done
	a3shaderProgramDeactivate();
	a3vertexDeactivateDrawable();
}


//-----------------------------------------------------------------------------
// release objects
// this is where the union array style comes in handy; don't need a single 
//	release statement for each and every object... just iterate and release!

// utility to unload textures
void a3demo_unloadTextures(a3_DemoState *demoState)
{
	a3_Texture *currentTex = demoState->texture,
		*const endTex = currentTex + demoStateMaxCount_texture;

	while (currentTex < endTex)
		a3textureRelease(currentTex++);
}

// utility to unload geometry
void a3demo_unloadGeometry(a3_DemoState *demoState)
{
	a3_BufferObject *currentBuff = demoState->drawDataBuffer,
		*const endBuff = currentBuff + demoStateMaxCount_drawDataBuffer;
	a3_VertexArrayDescriptor *currentVAO = demoState->vertexArray,
		*const endVAO = currentVAO + demoStateMaxCount_vertexArray;
	a3_VertexDrawable *currentDraw = demoState->drawable,
		*const endDraw = currentDraw + demoStateMaxCount_drawable;

	while (currentBuff < endBuff)
		a3bufferRelease(currentBuff++);
	while (currentVAO < endVAO)
		a3vertexArrayReleaseDescriptor(currentVAO++);
	while (currentDraw < endDraw)
		a3vertexReleaseDrawable(currentDraw++);
}


// utility to unload shaders
void a3demo_unloadShaders(a3_DemoState *demoState)
{
	a3_DemoStateShaderProgram *currentProg = demoState->shaderProgram,
		*const endProg = currentProg + demoStateMaxCount_shaderProgram;

	while (currentProg < endProg)
		a3shaderProgramRelease((currentProg++)->program);
}


//-----------------------------------------------------------------------------

// initialize non-asset objects
void a3demo_initScene(a3_DemoState *demoState)
{
	unsigned int i;
	const float cameraAxisPos = 15.0f;

	// all objects
	for (i = 0; i < demoStateMaxCount_sceneObject; ++i)
		a3demo_initSceneObject(demoState->sceneObject + i);

	a3demo_setCameraSceneObject(demoState->sceneCamera, demoState->cameraObject);
	for (i = 0; i < demoStateMaxCount_camera; ++i)
		a3demo_initCamera(demoState->sceneCamera + i);

	// cameras
	demoState->sceneCamera->znear = 1.00f;
	demoState->sceneCamera->zfar = 100.0f;
	demoState->sceneCamera->ctrlMoveSpeed = 10.0f;
	demoState->sceneCamera->ctrlRotateSpeed = 5.0f;
	demoState->sceneCamera->ctrlZoomSpeed = 5.0f;


	// camera's starting orientation depends on "vertical" axis
	// we want the exact same view in either case
	if (demoState->verticalAxis)
	{
		// vertical axis is Y
		demoState->sceneCamera->sceneObject->position.x = +cameraAxisPos;
		demoState->sceneCamera->sceneObject->position.y = +cameraAxisPos;
		demoState->sceneCamera->sceneObject->position.z = +cameraAxisPos;
		demoState->sceneCamera->sceneObject->euler.x = -30.0f;
		demoState->sceneCamera->sceneObject->euler.y = 45.0f;
		demoState->sceneCamera->sceneObject->euler.z = 0.0f;
	}
	else
	{
		// vertical axis is Z
		demoState->sceneCamera->sceneObject->position.x = +cameraAxisPos;
		demoState->sceneCamera->sceneObject->position.y = -cameraAxisPos;
		demoState->sceneCamera->sceneObject->position.z = +cameraAxisPos;
		demoState->sceneCamera->sceneObject->euler.x = 60.0f;
		demoState->sceneCamera->sceneObject->euler.y = 0.0f;
		demoState->sceneCamera->sceneObject->euler.z = 45.0f;
	}

	// same fovy to start
	demoState->sceneCamera->fovy = a3realSixty;


	// display physics info text
	demoState->displayAxes = 1;
	demoState->displayGrid = 1;
	demoState->displaySkybox = 0;
	demoState->displayPhysicsText = 0;


	// demo modes
	demoState->demoMode = 0;
	demoState->demoModeCount = 1;


	// drawable pointers for these objects
	demoState->rbDrawable[0] = demoState->draw_plane;

	demoState->rbDrawable[1] = demoState->draw_sphere;
	demoState->rbDrawable[2] = demoState->draw_sphere;
	demoState->rbDrawable[3] = demoState->draw_sphere;
	demoState->rbDrawable[4] = demoState->draw_sphere;

	demoState->rbDrawable[5] = demoState->draw_cylinder;
	demoState->rbDrawable[6] = demoState->draw_cylinder;
	demoState->rbDrawable[7] = demoState->draw_cylinder;
	demoState->rbDrawable[8] = demoState->draw_cylinder;

	demoState->rbDrawable[9] = demoState->draw_box;
	demoState->rbDrawable[10] = demoState->draw_box;
	demoState->rbDrawable[11] = demoState->draw_box;
	demoState->rbDrawable[12] = demoState->draw_box;


	// reset
	demoState->hitIndex = -1;
}


//-----------------------------------------------------------------------------

// the handle release callbacks are no longer valid; since the library was 
//	reloaded, old function pointers are out of scope!
// could reload everything, but that would mean rebuilding GPU data...
//	...or just set new function pointers!
void a3demo_refresh(a3_DemoState *demoState)
{
	a3_Texture *currentTex = demoState->texture,
		*const endTex = currentTex + demoStateMaxCount_texture;
	a3_BufferObject *currentBuff = demoState->drawDataBuffer,
		*const endBuff = currentBuff + demoStateMaxCount_drawDataBuffer;
	a3_VertexArrayDescriptor *currentVAO = demoState->vertexArray,
		*const endVAO = currentVAO + demoStateMaxCount_vertexArray;
	a3_DemoStateShaderProgram *currentProg = demoState->shaderProgram,
		*const endProg = currentProg + demoStateMaxCount_shaderProgram;

	while (currentTex < endTex)
		a3textureHandleUpdateReleaseCallback(currentTex++);
	while (currentBuff < endBuff)
		a3bufferHandleUpdateReleaseCallback(currentBuff++);
	while (currentVAO < endVAO)
		a3vertexArrayHandleUpdateReleaseCallback(currentVAO++);
	while (currentProg < endProg)
		a3shaderProgramHandleUpdateReleaseCallback((currentProg++)->program);
}


// confirm that all graphics objects were unloaded
void a3demo_validateUnload(const a3_DemoState *demoState)
{
	unsigned int handle;
	const a3_Texture *currentTex = demoState->texture,
		*const endTex = currentTex + demoStateMaxCount_texture;
	const a3_BufferObject *currentBuff = demoState->drawDataBuffer,
		*const endBuff = currentBuff + demoStateMaxCount_drawDataBuffer;
	const a3_VertexArrayDescriptor *currentVAO = demoState->vertexArray,
		*const endVAO = currentVAO + demoStateMaxCount_vertexArray;
	const a3_DemoStateShaderProgram *currentProg = demoState->shaderProgram,
		*const endProg = currentProg + demoStateMaxCount_shaderProgram;

	handle = 0;
	currentTex = demoState->texture;
	while (currentTex < endTex)
		handle += (currentTex++)->handle->handle;
	if (handle)
		printf("\n A3 Warning: One or more textures not released.");

	handle = 0;
	currentBuff = demoState->drawDataBuffer;
	while (currentBuff < endBuff)
		handle += (currentBuff++)->handle->handle;
	if (handle)
		printf("\n A3 Warning: One or more draw data buffers not released.");

	handle = 0;
	currentVAO = demoState->vertexArray;
	while (currentVAO < endVAO)
		handle += (currentVAO++)->handle->handle;
	if (handle)
		printf("\n A3 Warning: One or more vertex arrays not released.");

	handle = 0;
	currentProg = demoState->shaderProgram;
	while (currentProg < endProg)
		handle += (currentProg++)->program->handle->handle;
	if (handle)
		printf("\n A3 Warning: One or more shader programs not released.");
}


//-----------------------------------------------------------------------------
// MAIN LOOP

void a3demo_input(a3_DemoState *demoState, double dt)
{
	a3real ctrlRotateSpeed = 1.0f;
	a3real azimuth = 0.0f;
	a3real elevation = 0.0f;
	int rotatingCamera = 0, movingCamera = 0, changingParam = 0;

	// using Xbox controller
	if (a3XboxControlIsConnected(demoState->xcontrol))
	{
		// move and rotate camera using joysticks
		double lJoystick[2], rJoystick[2], lTrigger[1], rTrigger[1];
		a3XboxControlGetJoysticks(demoState->xcontrol, lJoystick, rJoystick);
		a3XboxControlGetTriggers(demoState->xcontrol, lTrigger, rTrigger);

		movingCamera = a3demo_moveSceneObject(demoState->camera->sceneObject, (float)dt * demoState->camera->ctrlMoveSpeed,
			(a3real)(rJoystick[0]),
			(a3real)(*rTrigger - *lTrigger),
			(a3real)(-rJoystick[1])
		);
		// rotate
		{
			ctrlRotateSpeed = 10.0f;
			azimuth = (a3real)(-lJoystick[0]);
			elevation = (a3real)(lJoystick[1]);

			// this really defines which way is "up"
			// mouse's Y motion controls pitch, but X can control yaw or roll
			// controlling yaw makes Y axis seem "up", roll makes Z seem "up"
			rotatingCamera = a3demo_rotateSceneObject(demoState->camera->sceneObject,
				ctrlRotateSpeed * (float)dt * demoState->camera->ctrlRotateSpeed,
				// pitch: vertical tilt
				elevation,
				// yaw/roll depends on "vertical" axis: if y, yaw; if z, roll
				demoState->verticalAxis ? azimuth : a3realZero,
				demoState->verticalAxis ? a3realZero : azimuth);
		}
	}

	// using mouse and keyboard
	else
	{
		// move using WASDEQ
		movingCamera = a3demo_moveSceneObject(demoState->camera->sceneObject, (float)dt * demoState->camera->ctrlMoveSpeed,
			(a3real)a3keyboardGetDifference(demoState->keyboard, a3key_D, a3key_A),
			(a3real)a3keyboardGetDifference(demoState->keyboard, a3key_E, a3key_Q),
			(a3real)a3keyboardGetDifference(demoState->keyboard, a3key_S, a3key_W)
		);
		if (a3mouseIsHeld(demoState->mouse, a3mouse_left))
		{
			azimuth = -(a3real)a3mouseGetDeltaX(demoState->mouse);
			elevation = -(a3real)a3mouseGetDeltaY(demoState->mouse);

			// this really defines which way is "up"
			// mouse's Y motion controls pitch, but X can control yaw or roll
			// controlling yaw makes Y axis seem "up", roll makes Z seem "up"
			rotatingCamera = a3demo_rotateSceneObject(demoState->camera->sceneObject,
				ctrlRotateSpeed * (float)dt * demoState->camera->ctrlRotateSpeed,
				// pitch: vertical tilt
				elevation,
				// yaw/roll depends on "vertical" axis: if y, yaw; if z, roll
				demoState->verticalAxis ? azimuth : a3realZero,
				demoState->verticalAxis ? a3realZero : azimuth);
		}
	}
}

void a3demo_update(a3_DemoState *demoState, double dt)
{
	const unsigned int graphicsObjectCount = sizeof(demoState->graphicsObjects) / sizeof(a3_DemoSceneObject);
	const unsigned int rigidBodyObjectCount = sizeof(demoState->physicsRigidbodies) / sizeof(a3_DemoSceneObject);
	const unsigned int particleObjectCount = sizeof(demoState->physicsParticles) / sizeof(a3_DemoSceneObject);
	unsigned int i;

	a3_DemoSceneObject *tmpObject;
	a3mat4 tmpScale = a3identityMat4;

	//a3_ConvexHullCollision testCollision[1];


	// update scene objects
	for (i = 0; i < graphicsObjectCount; ++i)
		a3demo_updateSceneObject(demoState->graphicsObjects + i);

	// update cameras
	for (i = 0; i < demoStateMaxCount_camera; ++i)
		a3demo_updateCameraViewProjection(demoState->camera + i);


	// grab stuff from physics
	// lock world
	if (a3physicsLockWorld(demoState->physicsWorld) > 0)
	{
		// copy state
		const a3_PhysicsWorldState worldState[1] = { *(demoState->physicsWorld->state) };

		// copy other

		// unlock
		a3physicsUnlockWorld(demoState->physicsWorld);

		// copy from physics world to demo state here
		for (i = 0; i < worldState->count_rb; ++i)
		{
			//	- convert rigid body state to graphics object state
			//	- apply graphics scale later (see below)

			demoState->physicsRigidbodies[i].modelMat = worldState->transform_rb[i];
		}
		demoState->rigidbodyObjects = a3minimum(i, rigidBodyObjectCount);

		// ditto
		for (i = 0; i < worldState->count_p; ++i)
		{
			//	- convert particle state to graphics object state
			//	- apply graphics scale later (see below)

			demoState->physicsParticles[i].modelMat = a3identityMat4;
			demoState->physicsParticles[i].modelMat.v3 = worldState->position_p[i];
		}
		demoState->particleObjects = a3minimum(i, particleObjectCount);
	}

	// apply scales
	for (i = 0, tmpObject = demoState->physicsRigidbodies + i; i < demoState->rigidbodyObjects; ++i, ++tmpObject)
		a3demo_applyScale_internal(demoState->physicsRigidbodies + i, tmpScale.m, 1);
	for (i = 0, tmpObject = demoState->physicsParticles + i; i < demoState->particleObjects; ++i, ++tmpObject)
		a3demo_applyScale_internal(demoState->physicsParticles + i, tmpScale.m, 1);
}



void a3demo_render(const a3_DemoState *demoState)
{
	const a3_VertexDrawable *currentDrawable;
	const a3_DemoStateShaderProgram *currentDemoProgram;

	const int useVerticalY = demoState->verticalAxis;


	// grid lines highlight
	// if Y axis is up, give it a greenish hue
	// if Z axis is up, a bit of blue
	const float gridColor[] = {
		0.15f,
		useVerticalY ? 0.25f : 0.20f,
		useVerticalY ? 0.20f : 0.25f,
		1.0f
	};


	// RGB
	const float rgba4[] = {
		1.0f, 0.0f, 0.0f, 1.0f,	// red
		0.0f, 1.0f, 0.0f, 1.0f,	// green
		0.0f, 0.0f, 1.0f, 1.0f,	// blue
		0.0f, 1.0f, 1.0f, 1.0f,	// cyan
		1.0f, 0.0f, 1.0f, 1.0f,	// magenta
		1.0f, 1.0f, 0.0f, 1.0f,	// yellow
		1.0f, 0.5f, 0.0f, 1.0f,	// orange
		0.0f, 0.5f, 1.0f, 1.0f,	// sky blue
		0.5f, 0.5f, 0.5f, 1.0f,	// solid grey
		0.5f, 0.5f, 0.5f, 0.5f,	// translucent grey
	},
	*const red = rgba4 + 0, *const green = rgba4 + 4, *const blue = rgba4 + 8,
	*const cyan = rgba4 + 12, *const magenta = rgba4 + 16, *const yellow = rgba4 + 20,
	*const orange = rgba4 + 24, *const skyblue = rgba4 + 28,
	*const grey = rgba4 + 32, *const grey_t = rgba4 + 36;


	// model transformations (if needed)
	const a3mat4 convertY2Z = {
		+1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, +1.0f, 0.0f,
		0.0f, -1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, +1.0f,
	};
	const a3mat4 convertZ2Y = {
		+1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, -1.0f, 0.0f,
		0.0f, +1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, +1.0f,
	};
	const a3mat4 convertZ2X = {
		0.0f, 0.0f, -1.0f, 0.0f,
		0.0f, +1.0f, 0.0f, 0.0f,
		+1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, +1.0f,
	};

	// final model matrix and full matrix stack
	a3mat4 modelMat = a3identityMat4, //modelMatInv = a3identityMat4, modelMatOrig = a3identityMat4,
		modelViewProjectionMat = a3identityMat4, modelViewMat = a3identityMat4, normalMat = a3identityMat4;

	// current scene object being rendered, for convenience
	const a3_DemoSceneObject *currentSceneObject;

	unsigned int i;


	// particle colors
	const float *rigidbodyColor[] = {
		red, green, blue, orange,
	};
	const float *particleColor[] = {
		cyan, magenta, yellow, skyblue,
	};


	// other default values
	const float defaultLength[1] = { 100.0f };


	// reset viewport and clear buffers
	a3framebufferDeactivateSetViewport(a3fbo_depth24, -demoState->frameBorder, -demoState->frameBorder, demoState->frameWidth, demoState->frameHeight);

	if (demoState->displaySkybox)
	{
		a3textureActivate(demoState->tex_skybox_clouds, a3tex_unit00);
		currentDemoProgram = demoState->prog_drawColorUnifTexture;
		a3shaderProgramActivate(currentDemoProgram->program);
		currentDrawable = demoState->draw_skybox;
		modelMat = demoState->verticalAxis ? a3identityMat4 : convertY2Z;
		modelMat.v3 = demoState->cameraObject->modelMat.v3;
		a3real4x4Product(modelViewProjectionMat.m, demoState->camera->viewProjectionMat.m, modelMat.m);
		a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uMVP, 1, modelViewProjectionMat.mm);
		a3shaderUniformSendFloat(a3unif_vec4, currentDemoProgram->uColor, 1, a3oneVec4.v);
		glCullFace(GL_FRONT);
		glDepthFunc(GL_ALWAYS);
		a3vertexActivateAndRenderDrawable(currentDrawable);
		glDepthFunc(GL_LEQUAL);
		glCullFace(GL_BACK);
	}
	else
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}


	// draw grid aligned to world
	if (demoState->displayGrid)
	{
		currentDemoProgram = demoState->prog_drawColorUnif;
		a3shaderProgramActivate(currentDemoProgram->program);
		currentDrawable = demoState->draw_grid;
		modelViewProjectionMat = demoState->camera->viewProjectionMat;
		if (useVerticalY)
			a3real4x4ConcatL(modelViewProjectionMat.m, convertZ2Y.m);
		a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uMVP, 1, modelViewProjectionMat.mm);
		a3shaderUniformSendFloat(a3unif_vec4, currentDemoProgram->uColor, 1, gridColor);
		a3vertexActivateAndRenderDrawable(currentDrawable);
	}


	// draw objects: 
	//	- correct "up" axis if needed
	//	- calculate full MVP matrix
	//	- move lighting objects' positions into object space
	//	- send uniforms
	//	- draw

	// activate common texture
	a3textureActivate(demoState->tex_checker, a3tex_unit00);

	// draw models
	currentDemoProgram = demoState->prog_drawLambert;
	a3shaderProgramActivate(currentDemoProgram->program);

	// send common uniforms
	a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uP, 1, demoState->sceneCamera->projectionMat.mm);

	// rigid body shapes
	for (i = 0; i < demoState->rigidbodyObjects; ++i)
	{
		currentSceneObject = demoState->physicsRigidbodies + i;
		modelMat = currentSceneObject->modelMat;
		currentDrawable = demoState->rbDrawable[i];
		a3vertexActivateDrawable(currentDrawable);
		a3real4x4Product(modelViewMat.m, demoState->cameraObject->modelMatInv.m, modelMat.m);
		a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uMV, 1, modelViewMat.mm);
		a3real4x4TransformInverse(normalMat.m, modelViewMat.m);
		a3real4x4Transpose(normalMat.m);
		a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uNrm, 1, normalMat.mm);
		a3shaderUniformSendFloat(a3unif_vec4, currentDemoProgram->uColor, 1, demoState->hitIndex == i ||
			demoState->colliding[i] ? red : a3oneVec4.v);
		a3vertexRenderActiveDrawable();
	}


	// visible points
	currentDemoProgram = demoState->prog_drawColorUnif;
	a3shaderProgramActivate(currentDemoProgram->program);
	a3shaderUniformSendFloat(a3unif_vec4, currentDemoProgram->uColor, 1, red);

	currentDrawable = demoState->draw_sphere;
	a3vertexActivateDrawable(currentDrawable);
	for (i = 0; i < demoState->particleObjects; ++i)
	{
		currentSceneObject = demoState->physicsParticles + i;
		modelMat = currentSceneObject->modelMat;
		a3real4x4Product(modelViewProjectionMat.m, demoState->camera->viewProjectionMat.m, modelMat.m);
		a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uMVP, 1, modelViewProjectionMat.mm);
		a3vertexRenderActiveDrawable();
	}


	// draw ray
	currentDemoProgram = demoState->prog_drawRay;
	a3shaderProgramActivate(currentDemoProgram->program);
	a3shaderUniformSendFloat(a3unif_vec4, currentDemoProgram->uColor, 1, demoState->hitIndex >= 0 ? green : red);
	a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uVP, 1, demoState->sceneCamera->viewProjectionMat.mm);

	currentDrawable = demoState->dummyDrawable;
	a3vertexActivateDrawable(currentDrawable);
	a3shaderUniformSendFloat(a3unif_vec4, currentDemoProgram->uRayOrigin_world, 1, demoState->ray->origin.v);
	a3shaderUniformSendFloat(a3unif_vec4, currentDemoProgram->uRayDirection_world, 1, demoState->ray->direction.v);
	a3shaderUniformSendFloat(a3unif_single, currentDemoProgram->uRayLength, 1, defaultLength);
	a3vertexRenderActiveDrawable();

	// draw hit
	if (demoState->hitIndex >= 0)
	{
		currentDemoProgram = demoState->prog_drawColorUnif;
		a3shaderProgramActivate(currentDemoProgram->program);

		currentDrawable = demoState->draw_sphere;
		a3vertexActivateDrawable(currentDrawable);

		a3shaderUniformSendFloat(a3unif_vec4, currentDemoProgram->uColor, 1, orange);
		a3real4x4SetScale(modelMat.m, 0.05f);

		modelMat.v3 = demoState->rayHit->hit0;
		a3real4x4Product(modelViewProjectionMat.m, demoState->camera->viewProjectionMat.m, modelMat.m);
		a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uMVP, 1, modelViewProjectionMat.mm);
		a3vertexRenderActiveDrawable();

		modelMat.v3 = demoState->rayHit->hit1;
		a3real4x4Product(modelViewProjectionMat.m, demoState->camera->viewProjectionMat.m, modelMat.m);
		a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uMVP, 1, modelViewProjectionMat.mm);
		a3vertexRenderActiveDrawable();
	}


	/*
		// teapot
		i = 0;
		currentDrawable = demoState->draw_teapot;
		currentSceneObject = demoState->planetObject + i;	// the sun is now a teapot

		modelMatOrig = currentSceneObject->modelMat;
		if (!useVerticalY)	// teapot's axis is Y
			a3real4x4Product(modelMat.m, modelMatOrig.m, convertY2Z.m);
		else
			modelMat = modelMatOrig;
		a3real4x4TransformInverseIgnoreScale(modelMatInv.m, modelMat.m);
		a3real4x4Product(modelViewProjectionMat.m, demoState->camera->viewProjectionMat.m, modelMat.m);

		a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uMVP, 1, modelViewProjectionMat.mm);
		a3shaderUniformSendFloat(a3unif_vec4, currentDemoProgram->uColor, 1, planetColor[i]);
		a3vertexActivateAndRenderDrawable(currentDrawable);
	*/

	if (demoState->displayAxes)
	{
		glDisable(GL_DEPTH_TEST);

		// draw coordinate axes in front of everything
		currentDemoProgram = demoState->prog_drawColor;
		a3shaderProgramActivate(currentDemoProgram->program);
		currentDrawable = demoState->draw_axes;
		a3vertexActivateDrawable(currentDrawable);

		// center of world
		modelViewProjectionMat = demoState->camera->viewProjectionMat;
		a3shaderUniformSendFloatMat(a3unif_mat4, 0, currentDemoProgram->uMVP, 1, modelViewProjectionMat.mm);
		a3vertexRenderActiveDrawable();

		glEnable(GL_DEPTH_TEST);
	}


	// deactivate things
	a3vertexDeactivateDrawable();
	a3shaderProgramDeactivate();


	// display planet names
	if (demoState->showText && demoState->textInit && demoState->displayPhysicsText)
	{
		const char *rigidbodyDescriptions[] = {
			"rb ground",
			"rb sphere A",
			"rb sphere B",
			"rb sphere C",
			"rb sphere D",
			"rb cylinder A",
			"rb cylinder B",
			"rb cylinder C",
			"rb cylinder D",
			"rb box A",
			"rb box B",
			"rb box C",
			"rb box D",
		};
		const char *particleDescriptions[] = {
			"particle",
		};
		a3vec4 pos_model = a3wVec4, pos_ndc;

		glDisable(GL_DEPTH_TEST);
		for (i = 0; i < demoState->rigidbodyObjects; ++i)
		{
			// transform world position to clip space
			pos_model = demoState->physicsRigidbodies[i].modelMat.v3;
			a3real4Real4x4Product(pos_ndc.v, demoState->camera->viewProjectionMat.m, pos_model.v);

			// perspective divide for NDC
			a3real3DivS(pos_ndc.v, pos_ndc.w);

			// displace
			pos_ndc.x += 0.05f;

			// display text
			a3textDraw(demoState->text, pos_ndc.x, pos_ndc.y, pos_ndc.z, magenta[0], magenta[1], magenta[2], 1.0f,
				rigidbodyDescriptions[i]);
		}
		for (i = 0; i < demoState->particleObjects; ++i)
		{
			pos_model = demoState->physicsParticles[i].modelMat.v3;
			a3real4Real4x4Product(pos_ndc.v, demoState->camera->viewProjectionMat.m, pos_model.v);
			a3real3DivS(pos_ndc.v, pos_ndc.w);
			pos_ndc.x += 0.05f;
			a3textDraw(demoState->text, pos_ndc.x, pos_ndc.y, pos_ndc.z, magenta[0], magenta[1], magenta[2], 1.0f,
				particleDescriptions[i]);
		}
		glEnable(GL_DEPTH_TEST);
	}


	// HUD
	if (demoState->textInit && demoState->showText)
	{
		// display mode info
		const char *demoModeText[] = {
			"Rigid body dynamics: ray picking",
		};

		const float col[3] = { 0.1f, 0.3f, 0.5f };

		glDisable(GL_DEPTH_TEST);

		a3textDraw(demoState->text, -0.98f, +0.90f, -1.0f, col[0], col[1], col[2], 1.0f,
			"Demo mode (%u / %u): ", demoState->demoMode + 1, demoState->demoModeCount);
		a3textDraw(demoState->text, -0.98f, +0.80f, -1.0f, col[0], col[1], col[2], 1.0f,
			"    %s", demoModeText[demoState->demoMode]);


		// physics reload
		a3textDraw(demoState->text, -0.98f, +0.70f, -1.0f, col[0], col[1], col[2], 1.0f,
			"RELOAD PHYSICS: 'R'");

		a3textDraw(demoState->text, -0.98f, +0.60f, -1.0f, col[0], col[1], col[2], 1.0f,
			"Toggle GRID:     'g' | Toggle AXES:     'x'");
		a3textDraw(demoState->text, -0.98f, +0.50f, -1.0f, col[0], col[1], col[2], 1.0f,
			"Toggle CAPTIONS: 'c' | Toggle SKYBOX:   'b'");


		// display controls
		if (a3XboxControlIsConnected(demoState->xcontrol))
		{
			a3textDraw(demoState->text, -0.98f, -0.50f, -1.0f, col[0], col[1], col[2], 1.0f,
				"Xbox controller camera control: ");
			a3textDraw(demoState->text, -0.98f, -0.60f, -1.0f, col[0], col[1], col[2], 1.0f,
				"    Left joystick = rotate | Right joystick, triggers = move");
		}
		else
		{
			a3textDraw(demoState->text, -0.98f, -0.50f, -1.0f, col[0], col[1], col[2], 1.0f,
				"Keyboard/mouse camera control: ");
			a3textDraw(demoState->text, -0.98f, -0.60f, -1.0f, col[0], col[1], col[2], 1.0f,
				"    Left click and drag = rotate | WASDEQ = move | wheel = zoom");
		}

		a3textDraw(demoState->text, -0.98f, -0.70f, -1.0f, col[0], col[1], col[2], 1.0f,
			"    Toggle demo mode:           ',' prev | next '.' ");
		a3textDraw(demoState->text, -0.98f, -0.80f, -1.0f, col[0], col[1], col[2], 1.0f,
			"    Toggle text display:        't' (toggle) | 'T' (alloc/dealloc) ");
		a3textDraw(demoState->text, -0.98f, -0.90f, -1.0f, col[0], col[1], col[2], 1.0f,
			"    Reload all shader programs: 'P' ****CHECK CONSOLE FOR ERRORS!**** ");

		glEnable(GL_DEPTH_TEST);
	}
}


//-----------------------------------------------------------------------------
