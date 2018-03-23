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
	
	a3_KeyboardInput.inl
	Basic keyboard input state tracking implementation.

	**DO NOT MODIFY THIS FILE**
*/


#ifdef __ANIMAL3D_KEYBOARDINPUT_H
#ifndef __ANIMAL3D_KEYBOARDINPUT_INL
#define __ANIMAL3D_KEYBOARDINPUT_INL


//-----------------------------------------------------------------------------

inline int a3keyboardGetState(const a3_KeyboardInput *keyboard, const a3_KeyboardKey key)
{
	if (keyboard)
		return (keyboard->key.key[key]);
	return a3input_invalid;
}

inline int a3keyboardIsHeld(const a3_KeyboardInput *keyboard, const a3_KeyboardKey key)
{
	if (keyboard)
		return (keyboard->key.key[key] && keyboard->key0.key[key]);
	return a3input_invalid;
}

inline int a3keyboardIsChanged(const a3_KeyboardInput *keyboard, const a3_KeyboardKey key)
{
	if (keyboard)
		return (keyboard->key.key[key] - keyboard->key0.key[key]);
	return a3input_invalid;
}

inline int a3keyboardGetStateASCII(const a3_KeyboardInput *keyboard, const char key)
{
	if (keyboard)
		return (keyboard->keyASCII.key[key]);
	return a3input_invalid;
}

inline int a3keyboardIsHeldASCII(const a3_KeyboardInput *keyboard, const char key)
{
	if (keyboard)
		return (keyboard->keyASCII.key[key] && keyboard->keyASCII0.key[key]);
	return a3input_invalid;
}

inline int a3keyboardIsChangedASCII(const a3_KeyboardInput *keyboard, const char key)
{
	if (keyboard)
		return (keyboard->keyASCII.key[key] - keyboard->keyASCII0.key[key]);
	return a3input_invalid;
}

inline int a3keyboardGetDifference(const a3_KeyboardInput *keyboard, const a3_KeyboardKey keyPositive, const a3_KeyboardKey keyNegative)
{
	if (keyboard)
		return (keyboard->key.key[keyPositive] - keyboard->key.key[keyNegative]);
	return a3input_invalid;
}

inline int a3keyboardGetDifferenceASCII(const a3_KeyboardInput *keyboard, const char keyPositive, const char keyNegative)
{
	if (keyboard)
		return (keyboard->keyASCII.key[keyPositive] - keyboard->keyASCII.key[keyNegative]);
	return a3input_invalid;
}

inline int a3keyboardSetState(a3_KeyboardInput *keyboard, const a3_KeyboardKey key, const a3_InputState state)
{
	if (keyboard)
		return (keyboard->key.key[key] = state);
	return a3input_invalid;
}

inline int a3keyboardSetStateASCII(a3_KeyboardInput *keyboard, const char key)
{
	if (keyboard)
		return (keyboard->keyASCII.key[key] = a3input_down);
	return a3input_invalid;
}

inline int a3keyboardUpdate(a3_KeyboardInput *keyboard)
{
	static const struct a3_KeyboardInput reset = { 0 };

	if (keyboard)
	{
		keyboard->key0 = keyboard->key;
		keyboard->keyASCII0 = keyboard->keyASCII;
		keyboard->keyASCII = reset.keyASCII;
		return 1;
	}
	return a3input_invalid;
}

inline int a3keyboardReset(a3_KeyboardInput *keyboard)
{
	static const struct a3_KeyboardInput reset = { 0 };

	if (keyboard)
	{
		keyboard->key = reset.key;
		keyboard->keyASCII = reset.keyASCII;
		return 1;
	}
	return a3input_invalid;
}


//-----------------------------------------------------------------------------


#endif	// !__ANIMAL3D_KEYBOARDINPUT_INL
#endif	// __ANIMAL3D_KEYBOARDINPUT_H