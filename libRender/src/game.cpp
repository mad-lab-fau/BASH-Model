#include "game.h"


//
// Game
//
Game::Game(sf::Time delta) :
	delta(delta), fps(0), ups(0), running(true) {
}

void Game::Run() {
	// http://gafferongames.com/game-physics/fix-your-timestep/
	sf::Clock clock;

	sf::Time currentTime = clock.getElapsedTime();
	sf::Time frameTime;
	sf::Time newTime;
	sf::Time accumulator;

	// http://noobtuts.com/cpp/frames-per-second
	sf::Clock fpsClock;
	int fpsCounter = 0;
	int upsCounter = 0;

	running = true;
	while (running) {

		newTime = clock.getElapsedTime();
		frameTime = newTime - currentTime;
		if (frameTime.asMilliseconds() > 250) {
			frameTime = sf::milliseconds(250);
		}
		currentTime = newTime;
		accumulator += frameTime;

		while (accumulator >= delta) {
			ProcessInput();
			Update(delta);
			upsCounter++;
			accumulator -= delta;
		}

		const float alpha = accumulator.asMicroseconds() / (delta.asMicroseconds() * 1.0f);

		Render(alpha, delta);
		fpsCounter++;
		if (fpsClock.getElapsedTime() > sf::milliseconds(1000)) {
			ups = upsCounter;
			fps = fpsCounter;
			fpsCounter = 0;
			upsCounter = 0;
			fpsClock.restart();
		}
	}
	Cleanup();
}

void Game::Stop() {
	running = false;
}

int Game::GetFPS() {
	return fps;
}

int Game::GetUPS() {
	return ups;
}

