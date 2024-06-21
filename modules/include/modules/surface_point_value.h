#pragma once

#include <geometrycentral/surface/surface_point.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

using namespace geometrycentral;
using namespace surface;

namespace modules {
  double surface_point_value(
      const VertexPositionGeometry &meshGeom,
      const VertexData<double> &function,
      const SurfacePoint &surfacePoint);
}
