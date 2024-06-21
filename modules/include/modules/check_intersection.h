#pragma once

#include <geometrycentral/surface/surface_point.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include <modules/surface_point_to_cartesian.h>

using namespace geometrycentral;
using namespace surface;

namespace modules {
  bool check_intersection(
    const ManifoldSurfaceMesh& mesh,
    const VertexPositionGeometry& meshGeom,
    const std::vector<SurfacePoint> &nodes,
    const std::vector<std::array<int, 2>> &segments,
    const std::vector<std::vector<SurfacePoint>> &edgeSurfacePoints,
    const std::vector<double> &edgeLengths
  );
}
