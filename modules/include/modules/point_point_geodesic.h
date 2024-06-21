#pragma once

#include <geometrycentral/surface/surface_point.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/manifold_surface_mesh.h>

#include <modules/exact_geodesics.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
std::tuple<
  std::vector<SurfacePoint>
> point_point_geodesic(
  ManifoldSurfaceMesh &mesh,
  VertexPositionGeometry &geometry,
  GeodesicAlgorithmExact &mmp,
  const SurfacePoint &p0,
  const SurfacePoint &p1
);
}
