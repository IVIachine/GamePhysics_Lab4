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
	
	a3_XboxControllerInput.h
	Basic Xbox controller input state tracking interface.

	**DO NOT MODIFY THIS FILE**
*/

#ifndef __ANIMAL3D_XBOXCONTROLLERINPUT_H
#define __ANIMAL3D_XBOXCONTROLLERINPUT_H


#include "a3_InputState.h"


#ifdef __cplusplus
extern "C"
{
#else	// !__cplusplus
	typedef struct a3_XboxControllerInput	a3_XboxControllerInput;
	typedef enum a3_XboxControllerButton	a3_XboxControllerButton;
	typedef enum a3_XboxControllerID		a3_XboxControllerID;
#endif	// __cplusplus


//-----------------------------------------------------------------------------

	// A3: Xbox controller state tracker; a wrapper for XInput's structure.
	//	members ctrl, ctrl0: wrapper for Xinput state with button, trigger and 
	//		joystick state, along with normalized versions of the triggers and 
	//		joysticks (for external flexibility)
	//	member rumble: vibration/rumble motor states
	//	member id: controller/player number (0-3)
	//	members connected, connected0: current and previous connection state
	struct a3_XboxControllerInput
	{
		struct {
			// internally represented as Xinput state
			unsigned long packet;
			unsigned short buttons;
			unsigned char lTrigger;
			unsigned char rTrigger;
			short lThumbX;
			short lThumbY;
			short rThumbX;
			short rThumbY;

			// normalized axis controls
			double lTrigger_unit;
			double rTrigger_unit;
			double lThumbX_unit;
			double lThumbY_unit;
			double rThumbX_unit;
			double rThumbY_unit;
		} ctrl, ctrl0;

		unsigned short rumble[2];
		short id;
		char connected, connected0;
	};


	// A3: Enumeration of controller buttons.
	enum a3_XboxControllerButton
	{
		a3xbox_null,
		a3xbox_DPAD_up = 0x0001,
		a3xbox_DPAD_down = 0x0002,
		a3xbox_DPAD_left = 0x0004,
		a3xbox_DPAD_right = 0x0008,
		a3xbox_start = 0x0010,
		a3xbox_back = 0x0020,
		a3xbox_leftThumb = 0x0040,
		a3xbox_rightThumb = 0x0080,
		a3xbox_leftBumper = 0x0100,
		a3xbox_rightBumper = 0x0200,
		a3xbox_A = 0x1000,
		a3xbox_B = 0x2000,
		a3xbox_X = 0x4000,
		a3xbox_Y = 0x8000,
	};


	// A3: Enumeration of controller IDs by player number.
	enum a3_XboxControllerID
	{
		a3xbox_player1,
		a3xbox_player2,
		a3xbox_player3,
		a3xbox_player4,
	};


//-----------------------------------------------------------------------------

	// A3: Initialize controller given ID and update that controller.
	//	param control: non-null pointer to Xbox controller state
	//	param id: player number index, determines which controller this is
	//	return: 1 if success
	//	return: a3input_invalid if invalid params
	int a3XboxControlSetID(a3_XboxControllerInput *control, const a3_XboxControllerID id);

	// A3: Get the state of a button.
	//	param control: non-null pointer to Xbox controller state
	//	param button: controller button to check the state of
	//	return:	1 if button is pressed
	//	return: 0 if button is not pressed
	//	return: a3input_invalid if invalid params
	int a3XboxControlGetState(const a3_XboxControllerInput *control, const a3_XboxControllerButton button);

	// A3: Check if button is held down from last update.
	//	param control: non-null pointer to Xbox controller state
	//	param button: controller button to check
	//	return:	1 if button is down and was also down last update
	//	return: 0 if button is not down or has changed
	//	return: a3input_invalid if invalid params
	int a3XboxControlIsHeld(const a3_XboxControllerInput *control, const a3_XboxControllerButton button);

	// A3: Check if button was pressed or released from last update.
	//	param control: non-null pointer to Xbox controller state
	//	param button: controller button to check
	//	return:	+1 if button changed from released to pressed
	//	return:	-1 if button changed from pressed to released
	//	return: 0 if key did not change
	//	return: a3input_invalid if invalid params
	int a3XboxControlIsChanged(const a3_XboxControllerInput *control, const a3_XboxControllerButton button);

	// A3: Get left trigger state as normalized value.
	//	param control: non-null pointer to Xbox controller state
	//	param trigger_out: non-null pointer to double to receive trigger state
	//	return: 1 if success, trigger_out holds value of trigger
	//	return: a3input_invalid if invalid params
	int a3XboxControlGetLeftTrigger(const a3_XboxControllerInput *control, double *trigger_out);

	// A3: Get right trigger state as normalized value.
	//	param control: non-null pointer to Xbox controller state
	//	param trigger_out: non-null pointer to double to receive trigger state
	//	return: 1 if success, trigger_out holds value of trigger
	//	return: a3input_invalid if invalid params
	int a3XboxControlGetRightTrigger(const a3_XboxControllerInput *control, double *trigger_out);

	// A3: Get both trigger states as normalized values.
	//	param control: non-null pointer to Xbox controller state
	//	param leftTrigger_out: non-null pointer to receive left trigger state
	//	param rightTrigger_out: non-null pointer to receive right trigger state
	//	return: 1 if success, pointers hold the respective trigger values
	//	return: a3input_invalid if invalid params
	int a3XboxControlGetTriggers(const a3_XboxControllerInput *control, double *leftTrigger_out, double *rightTrigger_out);

	// A3: Get left joystick state as normalized values.
	//	param control: non-null pointer to Xbox controller state
	//	param joystick_out: non-null array of 2 values to store joystick state
	//	return: 1 if success, trigger_out holds joystick X in position [0] 
	//		and holds joystick Y in position [1]
	//	return: a3input_invalid if invalid params
	int a3XboxControlGetLeftJoystick(const a3_XboxControllerInput *control, double *joystick_out);

	// A3: Get right joystick state as normalized values.
	//	param control: non-null pointer to Xbox controller state
	//	param joystick_out: non-null array of 2 values to store joystick state
	//	return: 1 if success, trigger_out holds joystick X in position [0] 
	//		and holds joystick Y in position [1]
	//	return: a3input_invalid if invalid params
	int a3XboxControlGetRightJoystick(const a3_XboxControllerInput *control, double *joystick_out);

	// A3: Get both joystick states as normalized values.
	//	param control: non-null pointer to Xbox controller state
	//	param leftJoystick_out: non-null array to store left joystick state
	//	param rightJoystick_out: non-null array to store right joystick state
	//	return: 1 if success, pointers hold the respective joystick values
	//	return: a3input_invalid if invalid params
	int a3XboxControlGetJoysticks(const a3_XboxControllerInput *control, double *leftJoystick_out, double *rightJoystick_out);

	// A3: Get connection status.
	//	param control: non-null pointer to Xbox controller state
	//	return: 1 if connected
	//	return: 0 if not connected
	//	return: a3input_invalid if invalid params
	int a3XboxControlIsConnected(const a3_XboxControllerInput *control);

	// A3: Check if connection was dropped or reconnected from last update.
	//	param control: non-null pointer to Xbox controller state
	//	return: +1 if reconnected
	//	return: -1 if disconnected
	//	return: 0 if no change
	//	return: a3input_invalid if invalid params
	int a3XboxControlConnectionChanged(const a3_XboxControllerInput *control);

	// A3: Set controller vibration motor speeds. Note: left motor is for 
	//		low-frequency rumble, right motor is for high-frequency rumble.
	//	param control: non-null pointer to Xbox controller state
	//	param leftMotor: the left motor speed
	//	param rightMotor: the right motor speed
	//	return: 1 if success
	//	return: 0 if failed (not connected)
	//	return: a3input_invalid if invalid params
	int a3XboxControlSetRumble(a3_XboxControllerInput *control, const unsigned short leftMotor, const unsigned short rightMotor);

	// A3: Update controller; checks connection state and sets button state.
	//	param control: non-null pointer to Xbox controller state
	//	return: 1 if connected and status updated
	//	return: 0 if not connected
	//	return: a3input_invalid if invalid param
	int a3XboxControlUpdate(a3_XboxControllerInput *control);

	// A3: Reset state; sets all control states to zero.
	//	param control: non-null pointer to Xbox controller state
	//	return: 1 if success
	//	return: a3input_invalid if invalid param
	int a3XboxControlReset(a3_XboxControllerInput *control);


//-----------------------------------------------------------------------------


#ifdef __cplusplus
}
#endif	// __cplusplus


#include "_inl/a3_XboxControllerInput.inl"


#endif	// !__ANIMAL3D_XBOXCONTROLLERINPUT_H