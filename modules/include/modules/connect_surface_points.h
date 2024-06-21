#pragma once

#include <geometrycentral/surface/surface_point.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>

#include <modules/exact_geodesics.h>

#include <vector>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
  std::tuple<
    std::vector<std::vector<SurfacePoint>>,
    std::vector<double>>
  connect_surface_points(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &meshGeom,
    const std::vector<SurfacePoint> &nodes,
    const std::vector<std::array<int, 2>> &segments);
}
