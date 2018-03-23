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
	
	a3_GeometryData.inl
	Inline definitions for geometry data and vertex descriptors.

	**DO NOT MODIFY THIS FILE**
*/

#ifdef __ANIMAL3D_GEOMETRYDATA_H
#ifndef __ANIMAL3D_GEOMETRYDATA_INL
#define __ANIMAL3D_GEOMETRYDATA_INL


//-----------------------------------------------------------------------------

inline int a3geometryGetVertexCount(const a3_GeometryData *geom)
{
	if (geom && geom->data)
		return geom->numVertices;
	return -1;
}

inline int a3geometryGetIndexCount(const a3_GeometryData *geom)
{
	if (geom && geom->data)
		return geom->numIndices;
	return -1;
}

inline int a3geometryGetVertexBufferSize(const a3_GeometryData *geom)
{
	if (geom && geom->data)
		return a3vertexStorageSpaceRequired(geom->vertexFormat, geom->numVertices);
	return 0;
}

inline int a3geometryGetIndexBufferSize(const a3_GeometryData *geom)
{
	if (geom && geom->data)
		return a3indexStorageSpaceRequired(geom->indexFormat, geom->numIndices);
	return 0;
}


//-----------------------------------------------------------------------------

inline int a3geometryGenerateVertexArray(a3_VertexArrayDescriptor *vertexArray_out, const a3_GeometryData *geom, a3_VertexBuffer *vertexBuffer, const unsigned int vertexBufferOffset)
{
	if (geom && geom->data)
		return a3vertexArrayCreateDescriptor(vertexArray_out, vertexBuffer, geom->vertexFormat, vertexBufferOffset);
	return -1;
}


//-----------------------------------------------------------------------------


#endif	// !__ANIMAL3D_GEOMETRYDATA_INL
#endif	// __ANIMAL3D_GEOMETRYDATA_H