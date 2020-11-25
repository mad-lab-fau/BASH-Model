
#include "util.h"
#include "io.h"
#include "osim.h"
#include "scapeWrapper.h"
#include "window.h"

//
// ArgumentParser
//
class ArgumentParser {
public:
	ArgumentParser(int argc, char** argv) {
		for (uint i = 1; i < argc; i++) {
			args.push_back(std::string(argv[i]));
		}
	}
	const std::string& Get(const std::string& option) const {
		std::vector<std::string>::const_iterator itr = std::find(args.begin(), args.end(), option);
		if (itr != args.end() && ++itr != args.end()) {
			return *itr;
		}
		return "";
	}
	bool Exists(const std::string& option) const {
		return std::find(args.begin(), args.end(), option) != args.end();
	}
private:
	std::vector<std::string> args;
};


//
// The main application
//
int main(int argc, char* argv[]) {

	PRINT("---------------------------------------------------------------");
	PRINT("Animation of 3D Human Surface Models for Biomechanical Analysis");
	PRINT("---------------------------------------------------------------");

	// parse arguments
	ArgumentParser args(argc, argv);

	// help
	if (args.Exists("-h") || args.Exists("--help")) {
		PRINT("Help");
		PRINT("Usage: ...");

		// TODO: explain hotkeys

		exit(EXIT_SUCCESS);
	}

	// OSIM
	if (args.Exists("--osim")) {
		const std::string& filenameOSIM = args.Get("--osim");
		if (!filenameOSIM.empty()) {
			Settings::GetInstance().inputModelOSIM = filenameOSIM;
			PRINT("OSIM file: " << Settings::GetInstance().inputModelOSIM);
		} else {
			PRINT_ERR("Argument .osim was empty.");
		}
	} else {
		PRINT_ERR("No .osim file was specified.");
	}

	// Scaling
	if (args.Exists("--scale")) {
		const std::string& filenameScale = args.Get("--scale");
		if (!filenameScale.empty()) {
			Settings::GetInstance().inputModelScale = filenameScale;
			PRINT("Scale file: " << Settings::GetInstance().inputModelScale);
		} else {
			PRINT_ERR("Argument ScaleFile was empty.");
		}
	} else {
		PRINT("No ScaleFile file was specified. Using generic model.");
	}

	// MOT
	if (args.Exists("--mot")) {
		const std::string& filenameMOT = args.Get("--mot");
		if (!filenameMOT.empty()) {
			Settings::GetInstance().inputModelMOT = filenameMOT;
			PRINT("MOT file: " << Settings::GetInstance().inputModelMOT);
		} else {
			PRINT_ERR("Argument .mot was empty.");
		}
	} else {
		PRINT_ERR("No .mot file was specified.");
	}

	// STO
	if (args.Exists("--sto")) {
		const std::string& filenameSTO = args.Get("--sto");
		if (!filenameSTO.empty()) {
			Settings::GetInstance().inputModelSTO = filenameSTO;
			Settings::GetInstance().visualizeMuscleActivity = true;
			PRINT("STO file: " << Settings::GetInstance().inputModelSTO);
		} else {
			Settings::GetInstance().visualizeMuscleActivity = false;
			Settings::GetInstance().showMuscles = false;
			PRINT("No .sto file was specified. Skipping visualization of muscle activity...");
		}
	} else {
		Settings::GetInstance().visualizeMuscleActivity = false;
		Settings::GetInstance().showMuscles = false;
		PRINT("No .sto file was specified. Skipping visualization of muscle activity...");
	}

	// Baseline Model
	if (args.Exists("--model")) {
		const std::string& filenameModel = args.Get("--model");
		if (!filenameModel.empty()) {
			Settings::GetInstance().baselineModelDir = filenameModel;
		} else {
			PRINT("No Baseline-Model directory was specified. Using default...");
		}
	} else {
		PRINT("No Baseline-Model directory was specified. Using default...");
	}
	PRINT("Baseline-Model: " << Settings::GetInstance().baselineModelDir);
	
	// Frames
	if (args.Exists("--frames")) {
		const std::string& limitFrames = args.Get("--frames");
		if (!limitFrames.empty()) {
			Settings::GetInstance().limitFrames = std::stoi(limitFrames);
			PRINT("Limited animation to " << Settings::GetInstance().limitFrames << " frames.");
		}
	}


	PRINT("---------------------------------------------------------------");

	// create necessary folders
	if (!CreateFolder(FILEPATH_CACHE_MAPPING)) {
		PRINT_ERR("Unable to create folder: " + std::string(FILEPATH_CACHE_MAPPING));
	}
	if (CACHE_ALL_MESHES) {
		std::string osimfile = GetFileFromPath(Settings::GetInstance().inputModelOSIM);
		std::string scaleFile = GetFileFromPath(Settings::GetInstance().inputModelScale);
		std::string motFile = GetFileFromPath(Settings::GetInstance().inputModelMOT);
		Settings::GetInstance().filepathModelCache = std::string(FILEPATH_CACHE_MESH) + osimfile + "/" + scaleFile + "/" + motFile + "/";
		if (!CreateFolder(Settings::GetInstance().filepathModelCache)) {
			PRINT_ERR("Unable to create folder: " + Settings::GetInstance().filepathModelCache);
		}
	}

	// initialize the global model mapper instance
	Model::GetInstance().InitModel();

	// create and start the render window for the application
	Window window(WindowMode::desktop);
	window.Run();
}