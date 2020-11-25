#include "window.h"

#include <render/input.h>

//
// Window
//
Window::Window(int wm) :
	Game(sf::milliseconds(15)), windowMode(wm), frameLimit(60), VSync(false) {

	sf::ContextSettings params;
	params.depthBits = 24;
	params.stencilBits = 8;
	params.antialiasingLevel = 0;
	params.majorVersion = 4;
	params.minorVersion = 4;

	int style;
	sf::VideoMode videoMode;
	if (windowMode == WindowMode::desktop) {
		style = sf::Style::Default;
		sf::VideoMode tmpVideoMode = sf::VideoMode::getDesktopMode();
		videoMode = sf::VideoMode(tmpVideoMode.width - 16, tmpVideoMode.height - 78);
	} else if (windowMode == WindowMode::fullscreen) {
		style = sf::Style::Fullscreen;
		videoMode = sf::VideoMode::getFullscreenModes()[0];
	} else { // WindowMode::window
		style = sf::Style::Default;
		videoMode = sf::VideoMode(DEFAULT_WINDOW_WIDTH, DEFAULT_WINDOW_HEIGHT);
	}
	renderWindow = new sf::RenderWindow(videoMode, "RenderWindow", style, params);

	windowWidth = renderWindow->getSize().x;
	windowHeight = renderWindow->getSize().y;

	if (windowMode == WindowMode::desktop) {
		renderWindow->setPosition({ 0, 0 });
	}

	if (frameLimit != -1) {
		renderWindow->setFramerateLimit(frameLimit);
	}
	PRINT_VAR(frameLimit);

	renderWindow->setVerticalSyncEnabled(VSync);
	PRINT_VAR(VSync);

	renderWindow->setActive(true);
	renderWindow->requestFocus();

	glewExperimental = GL_TRUE;
	glewInit();
	PRINT("GL_Version: " << glGetString(GL_VERSION));

	glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
	glEnable(GL_POINT_SMOOTH);

	sf::ContextSettings settings = renderWindow->getSettings();
	PRINT_VAR(settings.depthBits);
	PRINT_VAR(settings.stencilBits);
	PRINT_VAR(settings.antialiasingLevel);
	PRINT_VAR(settings.majorVersion);

	// World
	world = new World(windowWidth, windowHeight);
}

void Window::Cleanup() {
	SAFE_DELETE(world);
	SAFE_DELETE(renderWindow);
}

//
// Input
//
void Window::ProcessInput() {
	sf::Vector2i mousePos = sf::Mouse::getPosition(*renderWindow);
	int mouseWheel = 0;
	int mouseButton = MouseButton::None;

	sf::Event event;
	while (renderWindow->pollEvent(event)) {
		switch (event.type) {
		case sf::Event::Closed:
			PRINT("--- Window closed ---");
			Stop();
			break;
		case sf::Event::Resized:
			windowWidth = event.size.width;
			windowHeight = event.size.height;
			world->Resize(windowWidth, windowHeight);

			renderWindow->setView(sf::View(sf::FloatRect(0, 0, windowWidth, windowHeight)));

			PRINT_VAR(windowWidth);
			PRINT_VAR(windowHeight);
			break;
		case sf::Event::MouseWheelMoved:
			mouseWheel = event.mouseWheel.delta;
			break;

		case sf::Event::MouseButtonPressed:
			if (event.mouseButton.button == sf::Mouse::Button::Left) {
				mouseButton = MouseButton::LeftPressed;
			}
			if (event.mouseButton.button == sf::Mouse::Button::Right) {
				mouseButton = MouseButton::RightPressed;
			}
			if (event.mouseButton.button == sf::Mouse::Button::Middle) {
				mouseButton = MouseButton::MiddlePressed;
			}
			break;

		case sf::Event::MouseButtonReleased:
			if (event.mouseButton.button == sf::Mouse::Button::Left) {
				mouseButton = MouseButton::LeftReleased;
			}
			if (event.mouseButton.button == sf::Mouse::Button::Right) {
				mouseButton = MouseButton::RightReleased;
			}
			if (event.mouseButton.button == sf::Mouse::Button::Middle) {
				mouseButton = MouseButton::MiddleReleased;
			}
			break;

		default:
			break;
		}
	}


	// Process only focused window
	if (!renderWindow->hasFocus() || mousePos.x < 0 || mousePos.y < 0 || mousePos.x > windowWidth || mousePos.y > windowHeight) {
		return;
	}

	// Close Window
	if (Input::GetInstance().KeyPressed(KEY_CLOSE_WINDOW)) {
		PRINT("--- Window closed ---");
		Stop();
		return;
	}

	// Toggle VSync
	if (Input::GetInstance().KeyPressed(KEY_TOGGLE_VSYNC)) {
		VSync = !VSync;
		renderWindow->setVerticalSyncEnabled(VSync);
		PRINT_VAR(VSync);
	}

	world->ProcessInput();
	world->ProcessMouseInput(mousePos, mouseButton, mouseWheel);
}

//
// Update
//
void Window::Update(sf::Time delta) {
	world->Update(delta);
}

//
// Render
//
void Window::Render(float alpha, sf::Time delta) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	renderWindow->setTitle("RenderWindow: " + std::to_string(GetFPS()) + " fps / " + std::to_string(GetUPS()) + " ups");

	// Render stuff
	world->Render(alpha, delta);

	// SFML's OpenGL state
	glUseProgram(0);
	glEnable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	// renderWindow->resetGLStates();

	renderWindow->display();
}
