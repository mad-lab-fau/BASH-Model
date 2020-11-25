#include "input.h"

//
// Input
//
bool Input::MouseButtonPressed(sf::Mouse::Button button) {
	if (sf::Mouse::isButtonPressed(button)) {
		if (mouseButtons[button] == false) {
			mouseButtons[button] = true;
			return true;
		}
	} else {
		mouseButtons[button] = false;
	}
	return false;
}

bool Input::KeyPressed(sf::Keyboard::Key key) {
	if (sf::Keyboard::isKeyPressed(key)) {
		if (keys[key] == false) {
			keys[key] = true;
			return true;
		}
	} else {
		keys[key] = false;
	}
	return false;
}

bool Input::KeyPressedContinously(sf::Keyboard::Key key) {
	return sf::Keyboard::isKeyPressed(key);
}

