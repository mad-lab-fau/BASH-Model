#pragma once

#include "global.h"

#include <SFML/Window.hpp>

struct Renderer {
	virtual void Render(float alpha, sf::Time delta) = 0;
};

//
// Game
//
class Game {
private:
	const sf::Time delta;
	int fps, ups;
	bool running;

public:
	Game(sf::Time delta);
	~Game() {}

	virtual void Run();
	virtual void Stop();

	int GetFPS();
	int GetUPS();

	virtual void ProcessInput() = 0;
	virtual void Update(sf::Time delta) = 0;
	virtual void Render(float alpha, sf::Time delta) = 0;
	virtual void Cleanup() = 0;
};
