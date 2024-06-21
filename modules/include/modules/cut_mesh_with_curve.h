#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
  // warn: this function is leaky and do not rely on it!!
  // it is only used for visualization purpose
  std::tuple<
    std::vector<Eigen::MatrixXd>,
    std::vector<Eigen::MatrixXi>,
    Eigen::MatrixXd,
    Eigen::MatrixXi
  > cut_mesh_with_curve(
    ManifoldSurfaceMesh &_mesh,
    VertexPositionGeometry &_geometry,
    const std::vector<SurfacePoint> &_nodes,
    const std::vector<std::array<int, 2>> &_segments,
    const std::vector<std::vector<SurfacePoint>> &_segmentSurfacePoints
  );
}
