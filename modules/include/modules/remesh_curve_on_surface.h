#pragma once

#include <geometrycentral/surface/surface_point.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>

#include <tuple>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
std::tuple<
  std::vector<SurfacePoint>,
  std::vector<std::array<int, 2>>,
  std::vector<std::vector<SurfacePoint>>,
  std::vector<double>,
  std::vector<bool>
> remesh_curve_on_surface(
  ManifoldSurfaceMesh &mesh,
  VertexPositionGeometry &geometry,
  const std::vector<SurfacePoint> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<std::vector<SurfacePoint>> &segmentSurfacePoints,
  const std::vector<double> &segmentLengths,
  const std::vector<bool> &isFixedNode,
  const double &h
);
}
