#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
  std::vector<std::vector<Vector3>> medial_axis(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry,
    const std::vector<SurfacePoint> &nodes,
    const std::vector<std::array<int, 2>> &segments,
    const std::vector<Vector3> &cartesianCoords,
    const std::vector<Vector3> &nodeTangents,
    const std::vector<Vector3> &nodeNormals,
    const std::vector<Vector3> &nodeBitangents,
    const double maxRadius,
    const bool isGeodesic
  );
}
