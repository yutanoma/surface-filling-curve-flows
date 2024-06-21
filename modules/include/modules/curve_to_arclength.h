#pragma once

#include <vector>

#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// compute the arclength parametrization of a curve
// only works for relatively simple curves (e.g., single closed loop or single open loop)
namespace modules {
  std::vector<double> curve_to_arclength(
    const std::vector<Vector3>& nodes,
    const std::vector<std::array<int, 2>> &segments
  );
}
