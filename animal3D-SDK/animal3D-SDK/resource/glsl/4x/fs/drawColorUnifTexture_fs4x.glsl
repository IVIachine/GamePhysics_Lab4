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
	
	drawColorUnifTexture_fs4x.glsl
	Draw uniform solid color mixed with texture sample.

	**DO NOT MODIFY THIS FILE**
*/

#version 410

// ****TO-DO: 
//	1) declare varying to receive texture coordinate from vertex shader
//	2) declare uniform variable for color; see demo code for hints
//	3) declare uniform texture sampler; see demo code for hints
//	4) sample texture to get surface color
//	5) multiply texture sample by uniform color
//	6) copy resulting product to output

out vec4 rtFragColor;

void main()
{
	// DUMMY OUTPUT: all fragments are OPAQUE WHITE
	rtFragColor = vec4(1.0, 1.0, 1.0, 1.0);
}
