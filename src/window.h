#pragma once

#include "util.h"
#include "world.h"
#include "settings.h"

#include <render/game.h>
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

//
// Window
//
class Window : public Game {
private:
	int frameLimit;
	bool VSync;

	int windowMode;
	sf::RenderWindow* renderWindow;
	int windowWidth, windowHeight;

	World* world;
public:
	Window(int windowMode);

	void ProcessInput();
	void Update(sf::Time delta);
	void Render(float alpha, sf::Time delta);
	void Cleanup();
};


