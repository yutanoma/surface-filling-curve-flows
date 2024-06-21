#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
  void curve_evolution(
    std::vector<Vector2> &nodes,
    const std::vector<int> &globalNodes,
    const std::vector<int> &localNodes,
    const double timeStep
  );
}
