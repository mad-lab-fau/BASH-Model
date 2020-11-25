#pragma once

#include "global.h"

#include <SFML/Window.hpp>

enum MouseButton {
	None,
	LeftPressed,
	LeftReleased,
	RightPressed,
	RightReleased,
	MiddlePressed,
	MiddleReleased
};

//
// Input
//
class Input {
private:
	std::unordered_map<int, bool> mouseButtons;
	std::unordered_map<int, bool> keys;

	Input() {};
	Input(const Input&) = delete;
	Input& operator=(const Input&) = delete;
public:
	static Input& GetInstance() {
		static Input instance;
		return instance;
	}

	bool MouseButtonPressed(sf::Mouse::Button button);
	bool KeyPressed(sf::Keyboard::Key key);
	bool KeyPressedContinously(sf::Keyboard::Key key);
};
