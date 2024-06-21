#include "modules/scene_file.h"

namespace modules {
  using namespace std;

  void splitString(const std::string& str, std::vector<string> &cont, char delim) {
      std::stringstream ss(str);
      std::string token;
      while (std::getline(ss, token, delim)) {
        cont.push_back(token);
      }
  }

  std::string getDirectoryFromPath(std::string str) {
    using namespace std;
    vector<string> parts;
    splitString(str, parts, '/');

    int nParts = parts.size();
    if (nParts == 1) return "./";
    
    string path = "";

    for (int i = 0; i < nParts - 1; i++) {
        path = path + parts[i] + "/";
    }

    return path;
  }

  void processLine(SceneObject &scene, std::string directory, std::vector<std::string> &parts) {
    string key = parts[0];

    if (key == "#" || key == "//") {
      return;
    }

    if (key == "curve") {
      scene.curveFileName = directory + parts[1];
    } else if (key == "mesh") {
      scene.meshFileName = directory + parts[1];
    } else if (key == "curve_mesh") {
      scene.meshFileName = directory + parts[1];
    } else if (key == "dmat") {
      scene.dmatFilename = directory + parts[1];
    } else if (key == "scalar") {
      scene.scalarFileName = directory + parts[1];
    } else if (key == "radius") {
      scene.radius = stod(parts[1]);
      scene.h = scene.radius * igl::PI / 20;
      scene.rmax = scene.radius * 5;
    } else if (key == "timestep") {
      scene.timestep = stod(parts[1]);
    } else if (key == "h") {
      scene.h = stod(parts[1]);
    } else if (key == "p") {
      scene.p = stod(parts[1]);
    } else if (key == "q") {
      scene.q = stod(parts[1]);
    } else if (key == "rmax") {
      scene.rmax = stod(parts[1]);
    } else if (key == "field_aligned") {
      scene.w_fieldAlignedness = stod(parts[1]);
    } else if (key == "curxvature_aligned") {
      scene.w_curvatureAlignedness = stod(parts[1]);
    } else if (key == "bilaplacian") {
      scene.w_bilaplacian = stod(parts[1]);
    } else if (key == "varying_alpha") {
      scene.varyingAlpha = true;
    } else if (key == "geodesic_medial_axis") {
      scene.useGeodesicMedialAxis = true;
    } else if (key == "excecute_only") {
      scene.excecuteOnly = true;
    }
  }

  SceneObject read_scene(std::string filename) {
    string directory = getDirectoryFromPath(filename);

    ifstream inFile;
    inFile.open(filename);

    if (!inFile) {
      cerr << "Could not open file " << filename << endl;
      exit(1);
    }

    SceneObject scene;

    std::vector<std::string> parts;
    for (std::string line; std::getline(inFile, line ); ) {
      if (line == "" || line == "\n") continue;
      parts.clear();
      splitString(line, parts, ' ');
      processLine(scene, directory, parts);
    }

    inFile.close();
    return scene;
  }
}
