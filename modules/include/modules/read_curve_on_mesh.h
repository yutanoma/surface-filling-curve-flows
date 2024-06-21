#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
  std::tuple<
    std::vector<SurfacePoint>,
    std::vector<std::array<int, 2>>,
    std::vector<bool>
  >
  readCurveOnMesh(
      const std::string filename,
      ManifoldSurfaceMesh &mesh,
      VertexPositionGeometry &geometry);
}
