#include "osim.h"
#include "settings.h"
#include "io.h"

//
// OSIM
//
OSIM::OSIM(const std::string& filepath_osim, const std::string& filepath_scale, const std::string& filepath_mot, const std::string& filepath_sto) {
	try {
		// load the osim file
		osimModel = OpenSim::Model(filepath_osim);

		// Initialize the system and get the state representing the state system
		state = osimModel.initSystem();

		// load the scale factors from the .xml configuration file
		if (filepath_scale != "") {
			LoadScaleFactors(filepath_scale);
		}

		// load the motion file (.mot)
		if (filepath_mot != "") {
			motStorage = LoadStorageFile(filepath_mot);

			// load the muscle activation file (.sto)
			if (filepath_sto != "") {
				LoadMuscleActivationFile(filepath_sto);
			}
		}

	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}

int OSIM::GetNumFrames() {
	if (Settings::GetInstance().limitFrames != -1) {
		return Settings::GetInstance().limitFrames;
	} else {
		return motStorage.getSize();
	}
}

// loads a scaleSet from a .xml file
void OSIM::LoadScaleFactors(const std::string& filepath) {
	try {
		OpenSim::ScaleSet scaleSet;
		SimTK::State& currentState = osimModel.initSystem();
		osimModel.getMultibodySystem().realize(currentState, SimTK::Stage::Position);

		// load scaleSet from file
		scaleSet = OpenSim::ScaleSet(filepath);

		// fill scaleFactors vector
		for (uint i = 0; i < scaleSet.getSize(); i++) {
			scaleFactors[scaleSet[i].getSegmentName()] = toGLM(scaleSet[i].getScaleFactors());
		}

		// apply scaling to input model?
		if (APPLY_SCALEFACTORS_INPUTMODEL) {
			osimModel.scale(currentState, scaleSet, false);
			PRINT("Scaled the input model on given scale factors.");
		}
	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}

// get scale factors
std::map<std::string, glm::vec3> OSIM::GetScaleFactors() {
	return scaleFactors;
}


// TODO: use InverseKinematicsSolver instead of InverseKinematicsTool which only handles filepaths
void OSIM::InverseKinematics(const std::string& filepath_input_trc, const std::string& filepath_output_mot) {
	try {
		OpenSim::InverseKinematicsTool ik_tool;
		ik_tool.setModel(osimModel);
		ik_tool.setMarkerDataFileName(filepath_input_trc); // temporarily store the marker file
		ik_tool.setOutputMotionFileName(filepath_output_mot);
		ik_tool.setResultsDir(FILEPATH_CACHE_MAPPING);
		ik_tool.run();
	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}


// looad the muscle activation data from a .sto file
void OSIM::LoadMuscleActivationFile(const std::string& filepath) {
	try {
		// load .sto file
		OpenSim::TimeSeriesTable stoTable = OpenSim::STOFileAdapter().read(filepath);

		// number of time states has to match the loaded .mot file
		if (stoTable.getNumRows() != motStorage.getSize()) {
			PRINT_ERR("Number of time states has to match in file: " + filepath);
		}

		muscleActivations.resize(stoTable.getNumRows());
		for (int frameID = 0; frameID < stoTable.getNumRows(); frameID++) {
			double time = stoTable.getIndependentColumn().at(frameID);

			// match time
			if (time != motStorage.getStateVector(frameID)->getTime()) {
				PRINT_ERR("Time state values do not match in file: " + filepath);
			}

			// iterate through muscles
			OpenSim::Array<std::string> muscleNames;
			osimModel.getMuscles().getNames(muscleNames);

			for (int muscleID = 0; muscleID < muscleNames.size(); muscleID++) {
				std::string muscleName = muscleNames[muscleID];
				const OpenSim::Muscle& m = osimModel.getMuscles().get(muscleName);

				std::string muscleName_sto = stoTable.getColumnLabel(muscleID);
				if (muscleName != muscleName_sto) {
					PRINT_ERR("Muscle name '" + muscleName + "' does not match '" + muscleName_sto + "' in " + filepath);
				}

				double muscleActivation = stoTable.getRowAtIndex(frameID)[muscleID];

				muscleActivations.at(frameID)[muscleName] = muscleActivation;
			}
		}
	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}

// load OSIM storage format
OpenSim::Storage OSIM::LoadStorageFile(const std::string& filepath) {
	OpenSim::Storage storage;

	// load the storage file
	OpenSim::Storage storageFromFile = OpenSim::Storage(filepath);

	// if data is in degrees -> convert to radians
	if (storageFromFile.isInDegrees()) {
		osimModel.getSimbodyEngine().convertDegreesToRadians(storageFromFile);
		PRINT("Converted motion data from degrees to radians.");
	}

	// link motion data to the model and receive a new reorderd storage that matches the model internal states
	osimModel.formStateStorage(storageFromFile, storage, false);

	// per time states
	PRINT("Number of time states: " << storage.getSize());

	return storage;
}

void OSIM::SetStateToDefault() {
	try {
		osimModel.assemble(state);
	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}

// update the internal state of the model to a certain time frameID
void OSIM::SetStateToFrame(int frameID) {
	try {
		OpenSim::StateVector* stateVector = motStorage.getStateVector(frameID);
		state.setTime(stateVector->getTime());

		// convert string arrary -> vector
		OpenSim::Array<double> stateVectorData = stateVector->getData();
		SimTK::Vector stateVariableValues(stateVectorData.size());
		for (int i = 0; i < stateVectorData.size(); i++) {
			stateVariableValues.set(i, stateVectorData[i]);
		}

		// update osim model with state
		osimModel.setStateVariableValues(state, stateVariableValues);
		osimModel.assemble(state);
		//osimModel.realizeAcceleration(state);
		//osimModel.realizeDynamics(state);
		//osimModel.realizePosition(state);
		//osimModel.realizeReport(state);
		//osimModel.realizeTime(state);
		//osimModel.realizeVelocity(state);
		osimModel.equilibrateMuscles(state);
	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}


// get markers - in default pose
std::map<std::string, Data::Marker> OSIM::GetMarkers() {
	SetStateToDefault();
	return LoadMarkerData();
}
// get markers - for a certain frame
std::map<std::string, Data::Marker> OSIM::GetMarkers(int frameID) {
	SetStateToFrame(frameID);
	return LoadMarkerData();
}
// get markers
std::map<std::string, Data::Marker> OSIM::LoadMarkerData() {
	std::map<std::string, Data::Marker> markerData;
	try {
		const OpenSim::MarkerSet& markerSet = osimModel.getMarkerSet();
		OpenSim::Array<std::string> markerNames;
		markerSet.getNames(markerNames);

		// iterate through markers
		for (int markerID = 0; markerID < markerNames.size(); markerID++) {
			std::string markerName = markerNames[markerID];
			if (markerName.rfind("CP_", 0) == 0) { // drop markers named "CP_*"
				continue;
			}
			const OpenSim::Marker& m = markerSet.get(markerName);
			const OpenSim::Frame& f = m.getParentFrame();

			Data::Marker marker;
			marker.name = markerName;
			marker.globalPosition = toGLM(m.getLocationInGround(state));
			marker.attachedBody = f.getName();
			marker.localPosition = toGLM(m.findLocationInFrame(state, f));

			markerData[markerName] = marker;
		}
		return markerData;
	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}

// get bodies - in default pose
std::map<std::string, Data::Body> OSIM::GetBodies() {
	SetStateToDefault();
	return LoadBodyData();
}
// get bodies - for a certain frame
std::map<std::string, Data::Body> OSIM::GetBodies(int frameID) {
	SetStateToFrame(frameID);
	return LoadBodyData();
}
// get body bodies
std::map<std::string, Data::Body> OSIM::LoadBodyData() {
	std::map<std::string, Data::Body> bodyData;
	try {
		OpenSim::Array<std::string> bodyNames;
		osimModel.getBodySet().getNames(bodyNames);

		// iterate through all bodies
		for (int bodyID = 0; bodyID < bodyNames.size(); bodyID++) {
			std::string bodyName = bodyNames[bodyID];
			const OpenSim::Body& b = osimModel.getBodySet().get(bodyName);
			const SimTK::Transform& globalTransform = b.getTransformInGround(state);
			const SimTK::Transform& localTransform = b.findTransformInBaseFrame();

			// geometry
			const OpenSim::Geometry& frameGeometry = b.get_frame_geometry();
			const SimTK::Vec3& fg_scaleFactors = frameGeometry.get_scale_factors();

			int numAttachedGeometry = b.getProperty_attached_geometry().size();
			if (numAttachedGeometry > 0) {
				for (int geometryID = 0; geometryID < numAttachedGeometry; geometryID++) {
					const OpenSim::Geometry& attachedGeometry = b.get_attached_geometry(geometryID);

					const SimTK::Vec3& ag_scaleFactors = attachedGeometry.get_scale_factors();

					const OpenSim::AbstractProperty& a = attachedGeometry.getPropertyByName("mesh_file");
					std::string mesh_file = a.getValue<std::string>();

					// TODO: parse meshes with VTK
				}
			}

			Data::Body body;
			body.name = bodyName;
			body.globalTransform = toGLM(globalTransform);
			body.localTransform = toGLM(localTransform);

			bodyData[bodyName] = body;
		}
		return bodyData;

	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}

// get joints - in default pose
std::map<std::string, Data::Joint> OSIM::GetJoints() {
	SetStateToDefault();
	return LoadJointData();
}
// get joints - for a certain frame
std::map<std::string, Data::Joint> OSIM::GetJoints(int frameID) {
	SetStateToFrame(frameID);
	return LoadJointData();
}
// get joints
std::map<std::string, Data::Joint> OSIM::LoadJointData() {
	std::map<std::string, Data::Joint> jointData;

	try {
		const OpenSim::JointSet& jointSet = osimModel.getJointSet();
		OpenSim::Array<std::string> jointNames;
		jointSet.getNames(jointNames);

		// iterate through all joints
		for (int jointID = 0; jointID < jointNames.size(); jointID++) {
			std::string jointName = jointNames[jointID];
			const OpenSim::Joint& j = jointSet.get(jointName);

			const OpenSim::Frame& childFrame = j.getChildFrame();
			const OpenSim::Frame& parentFrame = j.getParentFrame().findBaseFrame();

			const SimTK::Transform& globalTransformChild = childFrame.getTransformInGround(state);
			const SimTK::Transform& globalTransformParent = parentFrame.getTransformInGround(state);

			const SimTK::Transform& localTransform = childFrame.findTransformBetween(state, parentFrame);

			// coordinates from joints (= raw joint angles like in the .mot file)
			for (int coordID = 0; coordID < j.numCoordinates(); coordID++) {
				const OpenSim::Coordinate& coordinate = j.get_coordinates(coordID);
				double coordinateValue = coordinate.getValue(state);
				//PRINT(coordinate.getName() << ": " << coordinateValue);
			}

			Data::Joint joint;
			joint.name = jointName;
			joint.childName = childFrame.getName();
			joint.parentName = parentFrame.getName();
			joint.globalTransformChild = toGLM(globalTransformChild);
			joint.globalTransformParent = toGLM(globalTransformParent);
			joint.localTransform = toGLM(localTransform);

			jointData[jointName] = joint;
		}
		return jointData;

	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}

// get muscles - in default pose
std::map<std::string, Data::Muscle> OSIM::GetMuscles() {
	SetStateToDefault();
	return LoadMuscleData();
}
// get muscles - for a certain frame
std::map<std::string, Data::Muscle> OSIM::GetMuscles(int frameID) {
	SetStateToFrame(frameID);

	// store the manually loaded muscle activations
	std::map<std::string, Data::Muscle> muscleData = LoadMuscleData();
	for (auto& m : muscleData) {
		Data::Muscle& muscle = m.second;

		if (Settings::GetInstance().visualizeMuscleActivity) {
			muscle.activation = muscleActivations.at(frameID).at(muscle.name);
		}
	}
	return muscleData;
}
// get muscles
std::map<std::string, Data::Muscle> OSIM::LoadMuscleData() {
	std::map<std::string, Data::Muscle> muscleData;

	try {
		OpenSim::Array<std::string> muscleNames;
		osimModel.getMuscles().getNames(muscleNames);

		// iterate through all muscles
		for (int muscleID = 0; muscleID < muscleNames.size(); muscleID++) {
			std::string muscleName = muscleNames[muscleID];
			const OpenSim::Muscle& m = osimModel.getMuscles().get(muscleName);

			Data::Muscle muscle;
			muscle.name = muscleName;
			muscle.activation = 0.0;
			muscle.maxForce = m.getMaxIsometricForce();

			// muscle line segments
			for (uint i = 0; i < m.getGeometryPath().getPathPointSet().getSize(); i++) {
				const SimTK::Vec3& p = m.getGeometryPath().getPathPointSet().get(i).getLocationInGround(state);
				muscle.lineSet.push_back(toGLM(p));
			}

			muscleData[muscleName] = muscle;
		}
		return muscleData;

	} catch (const std::exception & ex) {
		PRINT_ERR(ex.what());
	}
}

