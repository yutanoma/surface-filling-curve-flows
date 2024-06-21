#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/trace_geodesic.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
std::tuple<
  std::vector<SurfacePoint> /* nodes */,
  std::vector<std::array<int, 2>> /* segments */,
  std::vector<std::vector<SurfacePoint>> /* segment paths */,
  std::vector<double> /* segment lengths */,
  std::vector<bool> /* isFixedNode */,
  std::vector<std::vector<SurfacePoint>> /* retraction paths */
> surface_path_evolution(
  ManifoldSurfaceMesh& mesh,
  VertexPositionGeometry& geometry,
  const std::vector<SurfacePoint> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<std::vector<SurfacePoint>> &segmentSurfacePoints,
  const std::vector<double> &segmentLengths,
  const std::vector<bool> &isFixedNode,
  const double h,
  const std::vector<Vector3> &direction
);
}
