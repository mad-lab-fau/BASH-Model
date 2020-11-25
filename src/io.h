#pragma once

#include "util.h"


// helpful functions
inline std::string replaceAll(const std::string& input, const std::string& toSearch, const std::string& replaceStr) {
	std::string output = input;
	size_t pos = input.find(toSearch);
	while (pos != std::string::npos) {
		output.replace(pos, toSearch.size(), replaceStr);
		pos = output.find(toSearch, pos + replaceStr.size());
	}
	return output;
}
inline bool FileExists(const std::string& filepath) {
	std::ifstream file(filepath);
	if (file.good()) {
		file.close();
		return true;
	} else {
		file.close();
		return false;
	}
}
inline std::string GetFileFromPath(const std::string& path) {
	return path.substr(path.find_last_of("/\\") + 1);
}
inline bool CreateFolder(const std::string& path) {
	std::string mkdirCommand = "mkdir " + replaceAll(path, "/", "\\");
	system(mkdirCommand.c_str());
	return true;
}


//
// IO
//
class IO {
private:
	bool init = false;

	IO() {}
	IO(const IO&) = delete;
	IO& operator=(const IO&) = delete;

public:
	static IO& GetInstance() {
		static IO instance;
		return instance;
	}

	// .ply
	Data::Mesh ReadPLY(const std::string& filepath);
	void WritePLY(const std::string& filepath, const Data::Mesh& mesh);

	// .obj
	Data::Mesh ReadOBJ(const std::string& filepath);
	void WriteOBJ(const std::string& filepath, const Data::Mesh& mesh);

	// markers
	std::map<std::string, Data::Marker> ReadOBJ_Markers(const std::string& filepath);
	void WriteTRC_Markers(const std::string& filepath, const std::map<std::string, Data::Marker>& markers);

	// helpfull .obj
	void WriteOBJ_Lines(const std::string& filepath, const std::map<std::string, Line>& lines);
	void WriteOBJ_Points(const std::string& filepath, const std::map<std::string, glm::vec3>& points);

	// File content format: "numEntries\nEntry1 Entry2 Entry3..."
	template <typename T>
	inline std::vector<T> ReadVectorDataFromFile(const std::string& filepath) {
		std::ifstream file(filepath);
		if (!file.good()) {
			PRINT_ERR("Unable to read file " + filepath);
		}
		uint numEntries = 0;
		std::vector<T> entries;
		file >> numEntries;
		for (uint i = 0; i < numEntries; i++) {
			T entry;
			file >> entry;
			entries.push_back(entry);
		}
		if (numEntries != entries.size()) {
			PRINT_ERR("Number of entries does not match the read value.")
		}
		return entries;
	}
};

