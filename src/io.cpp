#include "io.h"


//
// load vertices and indices from a .ply file
//
Data::Mesh IO::ReadPLY(const std::string& filepath) {
	Data::Mesh mesh;
	std::ifstream file(filepath);
	if (!file.good()) {
		PRINT_ERR("Unable to read file " + filepath);
	}
	int numVertices = 0;
	int numFaces = 0;

	// read header
	std::string inputString;
	file >> inputString;
	while (inputString.compare("end_header") != 0) {
		file >> inputString;
		if (inputString.compare("element") == 0) {
			file >> inputString;
			if (inputString.compare("vertex") == 0) {
				file >> numVertices;
			} else if (inputString.compare("face") == 0) {
				file >> numFaces;
			}
		}
	}

	// read vertices
	double x, y, z;
	mesh.vertices.clear();
	for (long int i = 0; i < numVertices; ++i) {
		file >> x >> y >> z;
		mesh.vertices.push_back(glm::vec3(x, y, z));
	}

	// read faces
	long int faceSize;
	mesh.indices.clear();
	for (long int i = 0; i < numFaces; ++i) {
		file >> faceSize;
		file >> x >> y >> z;
		if (faceSize != 3) {
			PRINT_ERR("Only triangles are supported");
		}
		mesh.indices.push_back(x);
		mesh.indices.push_back(y);
		mesh.indices.push_back(z);
	}
	file.close();

	return mesh;
}


//
// write vertices and indices to a .ply file
//
void IO::WritePLY(const std::string& filepath, const Data::Mesh& mesh) {
	std::ofstream file(filepath, std::ios::trunc);
	if (file.fail()) {
		PRINT_ERR("Unable to open file " + filepath);
	}
	if (mesh.indices.size() % 3 != 0) {
		PRINT_ERR("Only triangles are supported");
	}
	int numFaces = mesh.indices.size() / 3;
	int numVertices = mesh.vertices.size();

	// write header
	file << "ply" << std::endl
		<< "format ascii 1.0" << std::endl
		<< "element vertex " << numVertices << std::endl
		<< "property float x" << std::endl
		<< "property float y" << std::endl
		<< "property float z" << std::endl
		<< "element face " << numFaces << std::endl
		<< "property list uchar int vertex_indices" << std::endl
		<< "end_header" << std::endl;

	// write vertices
	for (int i = 0; i < numVertices; ++i) {
		file << mesh.vertices.at(i).x << " " << mesh.vertices.at(i).y << " " << mesh.vertices.at(i).z << std::endl;
	}
	// write faces
	for (int i = 0; i < mesh.indices.size(); i += 3) {
		file << "3 " << mesh.indices.at(i) << " " << mesh.indices.at(i + 1) << " " << mesh.indices.at(i + 2) << std::endl;
	}
	file.close();
}


//
// write vertices, indices and optionally (vertex-)normals to a .obj file
//
Data::Mesh IO::ReadOBJ(const std::string& filepath) {
	Data::Mesh mesh;
	mesh.vertices.clear();
	mesh.normals.clear();
	mesh.indices.clear();

	std::ifstream file(filepath);
	if (!file.good()) {
		PRINT_ERR("Unable to read file " + filepath);
	}

	std::string line;
	while (file.good() && !file.eof() && std::getline(file, line)) {
		std::string key = "";
		std::stringstream lineStream(line);
		lineStream >> key;

		if (key == "v") { // vertices
			float x, y, z;
			lineStream >> x >> y >> z;
			mesh.vertices.push_back(glm::vec3(x, y, z));
		} else if (key == "vn") { // normals
			float x, y, z;
			lineStream >> x >> y >> z;
			mesh.normals.push_back(glm::vec3(x, y, z));
		} else if (key == "f") { // faces
			uint i0, i1, i2;
			lineStream >> i0 >> i1 >> i2;
			// indices (.obj starts counting at 1 instead of 0)
			mesh.indices.push_back(i0 - 1);
			mesh.indices.push_back(i1 - 1);
			mesh.indices.push_back(i2 - 1);
		}
	}

	file.close();
	return mesh;
};


//
// write vertices, indices and optionally (vertex-)normals to a .obj file
//
void IO::WriteOBJ(const std::string& filepath, const Data::Mesh& mesh) {
	std::ofstream file(filepath, std::ios::trunc);
	if (file.fail()) {
		PRINT_ERR("Unable to open file " + filepath);
	}
	if (mesh.indices.size() % 3 != 0) {
		PRINT_ERR("Only triangles are supported");
	}
	file << "# numVertices=" << mesh.vertices.size() << "\n";
	file << "# numFaces=" << mesh.indices.size() / 3 << "\n";

	// write vertices (and normals if given)
	for (int vertexID = 0; vertexID < mesh.vertices.size(); vertexID++) {
		file << "v " << mesh.vertices.at(vertexID).x << " " << mesh.vertices.at(vertexID).y << " " << mesh.vertices.at(vertexID).z << "\n";

		if (mesh.normals.size() > 0) {
			if (mesh.normals.size() != mesh.vertices.size()) {
				PRINT_ERR("Number of normals has to match the number of vertices");
			}
			file << "vn " << mesh.normals.at(vertexID).x << " " << mesh.normals.at(vertexID).y << " " << mesh.normals.at(vertexID).z << "\n";
		}
	}

	// write indices (.obj starts counting at 1 instead of 0)
	for (int i = 0; i < mesh.indices.size(); i += 3) {
		file << "f " << mesh.indices.at(i) + 1 << " " << mesh.indices.at(i + 1) + 1 << " " << mesh.indices.at(i + 2) + 1 << "\n";
	}

	file.close();
};


//
// read an .obj file to receive stored markers
//
std::map<std::string, Data::Marker> IO::ReadOBJ_Markers(const std::string& filepath) {
	std::map<std::string, Data::Marker> markers;

	std::ifstream file(filepath);
	if (!file.good()) {
		PRINT_ERR("Unable to read file " + filepath);
	}
	std::string line;
	while (file.good() && !file.eof() && std::getline(file, line)) {
		std::string key = "";
		std::stringstream lineStream(line);
		lineStream >> key;

		if (key == "o") {
			std::string markerName;
			lineStream >> markerName;

			// parse line after "o"
			std::getline(file, line);
			lineStream = std::stringstream(line);
			lineStream >> key;

			if (key == "v") {
				double x, y, z;
				lineStream >> x >> y >> z; // marker coordinates
				Data::Marker marker;
				marker.name = markerName;
				marker.globalPosition = glm::vec3(x, y, z);
				markers[markerName] = marker;
			}
		}
	}
	return markers;
}


//
// write markers to a .trc file (standard extension to store MoCap-Markers)
//
void IO::WriteTRC_Markers(const std::string& filepath, const std::map<std::string, Data::Marker>& markers) {
	std::ofstream file(filepath, std::ios::trunc);
	if (file.fail()) {
		PRINT_ERR("Unable to open file " + filepath);
	}
	std::string line1 = "Frame#\tTime\t";
	std::string line2 = "\t\t";
	std::string line3 = "1\t0\t";
	int counter = 1;
	for (const auto& m : markers) {
		Data::Marker marker = m.second;
		line1 += marker.name + "\t" + "\t" + "\t";
		line2 += "X" + std::to_string(counter) + "\t" + "Y" + std::to_string(counter) + "\t" + "Z" + std::to_string(counter) + "\t";
		line3 += std::to_string(marker.globalPosition.x) + "\t" + std::to_string(marker.globalPosition.y) + "\t" + std::to_string(marker.globalPosition.z) + "\t";
		counter++;
	}
	file << "PathFileType\t4\t(X / Y / Z)\tmarkers.trc\nDataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n200\t200\t1\t" << markers.size() << "\tm\t200\t1\t1" << std::endl;
	file << line1 << std::endl << line2 << std::endl << line3 << std::endl;

	file.close();
}



//
// writes a tupel (string, vec3) to an .obj file
//
void IO::WriteOBJ_Points(const std::string& filepath, const std::map<std::string, glm::vec3>& points) {
	std::ofstream file(filepath, std::ios::trunc);
	if (file.fail()) {
		PRINT_ERR("Unable to open file " + filepath);
	}
	file << "# in Blender: Set Origin to Geometry & Visualize names" << std::endl;
	int counter = 1;
	for (auto& p : points) {
		file << "o " << p.first << std::endl;
		file << "v " << p.second.x << " " << p.second.y << " " << p.second.z << std::endl;
		file << "f " << counter++ << std::endl; // to recognize the points as geometry
	}
	file.close();
}

//
// writes a tupel (string, line) to an .obj file
//
void IO::WriteOBJ_Lines(const std::string& filepath, const std::map<std::string, Line>& lines) {
	std::ofstream file(filepath, std::ios::trunc);
	if (file.fail()) {
		PRINT_ERR("Unable to open file " + filepath);
	}
	file << "# in Blender: Set Origin to Geometry & Visualize names" << std::endl;
	int counter = 1;
	for (auto& l : lines) {
		file << "o " << l.first << std::endl;
		file << "v " << l.second.A.x << " " << l.second.A.y << " " << l.second.A.z << std::endl;
		file << "v " << l.second.B.x << " " << l.second.B.y << " " << l.second.B.z << std::endl;
		file << "l " << counter++ << " ";
		file << counter++ << std::endl;
	}
	file.close();
}


