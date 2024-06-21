#pragma once

#include <geometrycentral/utilities/vector2.h>

using namespace geometrycentral;

namespace modules {
// refine the curve s.t. all the lengths of the curve segments are between h and 2h
std::vector<std::vector<Vector2>> curve_remesh_2d(
  const std::vector<std::vector<Vector2>> &curves,
  // add node if distance between two nodes is larger than h
  const double &h
);
}
