#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <igl/PI.h>

namespace modules {
  enum CurveType {
    Closed,
    Open
  };

  struct SceneObject {
    std::string curveFileName = "";
    std::string fixedCurveFileName = "";
    std::string meshFileName = "";
    std::string dmatFilename = "";
    std::string scalarFileName = "";
    CurveType curveType = CurveType::Closed;
    std::vector<int> startIds;
    std::vector<int> endIds;
    double velocity = 1.0;
    double radius = 0.1;
    double alpha = 1;
    double timestep = 1.;
    double h = igl::PI * radius / 20;
    double p = 2;
    double q = 2;
    double rmax = radius * 5;
    double w_fieldAlignedness = 0;
    double w_curvatureAlignedness = 0;
    double w_bilaplacian = 0;
    bool varyingAlpha = false;
    bool useGeodesicMedialAxis = false;
    bool excecuteOnly = false;
  };

  void splitString(const std::string &str, std::vector<std::string> &cont, char delim = ' ');

  SceneObject read_scene(std::string filename);
}
