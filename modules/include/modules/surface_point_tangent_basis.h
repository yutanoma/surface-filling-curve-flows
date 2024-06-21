#pragma once

#include <geometrycentral/surface/surface_point.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include <iostream>

using namespace geometrycentral;
using namespace surface;

namespace modules {
  std::tuple<Vector3, Vector3> surface_point_tangent_basis(const VertexPositionGeometry &meshGeom, const SurfacePoint &surfacePoint);
}
