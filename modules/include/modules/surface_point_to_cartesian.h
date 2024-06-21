#pragma once

#include <geometrycentral/surface/surface_point.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

#include <TinyAD/ScalarFunction.hh>

using namespace geometrycentral;
using namespace surface;

namespace modules {
  // converts a SurfacePoint list into Vector3 list
  std::vector<Vector3> surface_point_to_cartesian(
    const ManifoldSurfaceMesh& mesh,
    const VertexPositionGeometry& meshGeom,
    const std::vector<SurfacePoint> &surfacePoints
  );
}
