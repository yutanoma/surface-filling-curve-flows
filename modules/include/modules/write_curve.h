#pragma once

#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
void write_curve(
  const std::string &fname,
  std::vector<Vector3> &nodes,
  std::vector<std::array<int, 2>> &segments
);
}