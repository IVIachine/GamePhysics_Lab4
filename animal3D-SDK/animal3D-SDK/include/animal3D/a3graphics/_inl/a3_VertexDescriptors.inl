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
	
	a3_VertexDescriptors.inl
	Inline definitions for vertex descriptors.

	**DO NOT MODIFY THIS FILE**
*/

#ifdef __ANIMAL3D_VERTEXDESCRIPTORS_H
#ifndef __ANIMAL3D_VERTEXDESCRIPTORS_INL
#define __ANIMAL3D_VERTEXDESCRIPTORS_INL


//-----------------------------------------------------------------------------
// internal utility to help pad

#define a3pad(n, p)	( (n) + ( (p) - (n) % (p) ) % (p) )
#define a3pad4(n)	a3pad(n, 4)


//-----------------------------------------------------------------------------

inline int a3vertexAttribGetElementsPerAttrib(const a3_VertexAttributeType attribType)
{
	static const char elementsPerAttrib[] = { 0, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4 };
	return elementsPerAttrib[attribType];
}

inline int a3vertexAttribGetBytesPerElement(const a3_VertexAttributeType attribType)
{
	static const char bytesPerElement[] = { 0, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8 };
	return bytesPerElement[attribType];
}

inline int a3indexGetBytesPerIndex(const a3_IndexType indexType)
{
	static const char bytesPerIndex[] = { 0, 1, 2, 4 };
	return bytesPerIndex[indexType];
}

inline int a3vertexStorageSpaceRequired(const a3_VertexFormatDescriptor *vertexFormat, const unsigned int count)
{
	if (vertexFormat)
	{
		// padding may be required
		const unsigned int bytes = (vertexFormat->vertexSize * count);
		return a3pad4(bytes);
	}
	return -1;
}

inline int a3indexStorageSpaceRequired(const a3_IndexFormatDescriptor *indexFormat, const unsigned int count)
{
	if (indexFormat)
	{
		// padding may be required
		const unsigned int bytes = (indexFormat->indexSize * count);
		return a3pad4(bytes);
	}
	return -1;
}


//-----------------------------------------------------------------------------

inline int a3vertexAttribCreateDescriptor(a3_VertexAttributeDescriptor *attrib_out, const a3_VertexAttributeName attribName, const a3_VertexAttributeType attribType)
{
	if (attrib_out)
	{
		attrib_out->name = attribName;
		attrib_out->type = attribType;
		return attribType;
	}
	return -1;
}

inline int a3vertexAttribDataCreateDescriptor(a3_VertexAttributeDataDescriptor *attrib_out, const a3_VertexAttributeName attribName, const void *attribData)
{
	if (attrib_out)
	{
		attrib_out->name = attribName;
		attrib_out->data = attribData;
		return 1;
	}
	return -1;
}

inline int a3vertexFormatHasAttribute(const a3_VertexFormatDescriptor *vertexFormat, const a3_VertexAttributeName attribName)
{
	if (vertexFormat)
		return (vertexFormat->attribType[attribName]);
	return -1;
}

inline int a3indexFormatIsActive(const a3_IndexFormatDescriptor *indexFormat)
{
	if (indexFormat)
		return (indexFormat->indexType);
	return -1;
}


//-----------------------------------------------------------------------------


#endif	// !__ANIMAL3D_VERTEXDESCRIPTORS_INL
#endif	// __ANIMAL3D_VERTEXDESCRIPTORS_H