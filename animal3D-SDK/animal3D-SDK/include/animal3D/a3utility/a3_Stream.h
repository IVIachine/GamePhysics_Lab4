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
	
	a3_Stream.h
	Basic stream reading and writing interface; can be used to load and save 
		files as well as copy data to and from strings.

	**DO NOT MODIFY THIS FILE**
*/

#ifndef __ANIMAL3D_STREAM_H
#define __ANIMAL3D_STREAM_H



#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_Stream		a3_Stream;
	typedef struct a3_FileStream	a3_FileStream;
#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// A3: Basic string stream info structure.
	//	member contents: cstring of stream contents
	//	member ptr: iterator for streaming
	//	member length: length of contents in bytes
	//	member offset: offset of ptr from contents start
	struct a3_Stream
	{
		const char *contents;
		char *ptr;
		unsigned int length;
		unsigned int offset;
	};


	// A3: Basic wrapper for file streaming.
	//	member stream: alias for stdio's file pointer
	//	member bytes: total number of bytes read/written since opening stream
	//	member lastStream: number of bytes read/written in the last streaming call
	//	member mode: stream usage: read/write, for internal use
	struct a3_FileStream
	{
		void *stream;
		unsigned int bytes;
		unsigned int lastStream;
		int mode;
	};


	// A3: General form of external string streaming functions: 
	//	Read: takes a pointer to an object and a constant cstring
	//	Write: takes a pointer to a constant object and a cstring
	//	return: number of bytes read/written
	typedef int(*a3_StreamReadFunc)(void *, const char *);
	typedef int(*a3_StreamWriteFunc)(const void *, char *);


	// A3: General form of external file streaming functions: 
	//	Read: takes a pointer to an object and a file stream pointer
	//	Write: takes a pointer to a constant object and a file stream pointer
	//	return: number of bytes read/written
	typedef int(*a3_FileStreamReadFunc)(void *, const a3_FileStream *);
	typedef int(*a3_FileStreamWriteFunc)(const void *, const a3_FileStream *);


//-----------------------------------------------------------------------------

	// A3: Allocate an empty stream.
	//	param stream_out: non-null pointer to unused stream info structure
	//	param size: non-zero number of bytes
	//	return: number of bytes allocated if success
	//	return: 0 if failed to allocate
	//	return: -1 if invalid params or already in-use
	int a3streamAllocate(a3_Stream *stream_out, const unsigned int size);

	// A3: Load contents of file.
	//	param stream_out: non-null pointer to unused stream info structure
	//	param filePath: non-null, non-empty cstring of relative / absolute path
	//	return: file length if success
	//	return: 0 if file load failed
	//	return: -1 if invalid params or already in-use
	int a3streamLoadContents(a3_Stream *stream_out, const char *filePath);

	// A3: Write contents to file.
	//	param stream: non-null pointer to stream descriptor to write
	//	param filePath: non-null, non-empty cstring of relative / absolute path
	//	return: file length if success
	//	return: -1 if invalid params or stream does not have contents
	int a3streamSaveContents(const a3_Stream *stream, const char *filePath);

	// A3: Use external stream function to read data.
	//	param stream: non-null pointer to stream descriptor
	//	param object: non-null pointer to object that is streaming
	//	param streamReadFunc: non-null function to use for streaming
	//	return: number of bytes streamed
	//	return: -1 if invalid params or function returns invalid
	int a3streamObjectRead(a3_Stream *stream, void *object, const a3_StreamReadFunc streamReadFunc);

	// A3: Use external stream function to store data.
	//		NOTE: Enough space must already be allocated!
	//	param stream: non-null pointer to stream descriptor
	//	param object: non-null pointer to object that is streaming
	//	param streamWriteFunc: non-null function to use for streaming
	//	return: number of bytes streamed
	//	return: -1 if invalid params or function returns invalid
	int a3streamObjectWrite(a3_Stream *stream, const void *object, const a3_StreamWriteFunc streamWriteFunc);

	// A3: Reset stream's pointer without releasing contents.
	//	param stream: non-null pointer to stream descriptor
	//	return: old offset
	//	return: -1 if invalid param
	int a3streamReset(a3_Stream *stream);

	// A3: Release stream contents; only use this on streams loaded using the 
	//		above functions to avoid heap corruption.
	//	param stream: non-null pointer to file descriptor to release
	//	return: file length if success
	//	return: -1 if invalid params or stream does not have contents
	int a3streamReleaseContents(a3_Stream *stream);


//-----------------------------------------------------------------------------

	// A3: Open file stream for reading; closes any existing stream.
	//	param fileStream: non-null pointer to file stream wrapper
	//	param filePath: non-null, non-empty cstring of file path to read from
	//	return: 1 if success
	//	return: 0 if failed
	//	return: -1 if invalid params
	int a3fileStreamOpenRead(a3_FileStream *fileStream, const char *filePath);

	// A3: Open file stream for writing; closes any existing stream.
	//	param fileStream: non-null pointer to file stream wrapper
	//	param filePath: non-null, non-empty cstring of file path to write to
	//	return: 1 if success
	//	return: 0 if failed
	//	return: -1 if invalid params
	int a3fileStreamOpenWrite(a3_FileStream *fileStream, const char *filePath);

	// A3: Read object from file using file stream.
	//	param fileStream: non-null pointer to opened file stream in read mode
	//	param object: non-null pointer to object that is streaming
	//	param streamReadFunc: non-null function to use for streaming
	//	return: number of bytes read
	//	return: -1 if invalid params or streaming function returns invalid
	int a3fileStreamReadObject(a3_FileStream *fileStream, void *object, const a3_FileStreamReadFunc streamReadFunc);

	// A3: Write object to file using file stream.
	//	param fileStream: non-null pointer to opened file stream in write mode
	//	param object: non-null pointer to constant object that is streaming
	//	param streamWriteFunc: non-null function to use for streaming
	//	return: number of bytes written
	//	return: -1 if invalid params or streaming function returns invalid
	int a3fileStreamWriteObject(a3_FileStream *fileStream, const void *object, const a3_FileStreamWriteFunc streamWriteFunc);

	// A3: Close file stream; counters are retained.
	//	param fileStream: non-null pointer to file stream wrapper
	//	return: total number of bytes read/written if success
	//	return: -1 if invalid params
	int a3fileStreamClose(a3_FileStream *fileStream);


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#endif	// !__ANIMAL3D_STREAM_H