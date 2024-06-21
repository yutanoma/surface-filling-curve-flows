#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"

namespace modules {
std::tuple<
  std::vector<geometrycentral::Vector3>, // node list on medial axis
  std::vector<std::array<int, 2>>, // segment list on medial axis
  std::vector<double>, // radius list on medial axis
  std::vector<std::vector<std::vector<int>>>, // node list to the nearest medial axis node
  Eigen::MatrixXd, // V of underlying mesh
  Eigen::MatrixXi, // F of underlying mesh
  Eigen::MatrixXi // T of underlying mesh
> get_medial_axis_curve_euclidian(
  geometrycentral::surface::ManifoldSurfaceMesh &mesh,
  geometrycentral::surface::VertexPositionGeometry &geometry,
  const Eigen::MatrixXd &meshV,
  const Eigen::MatrixXi &meshF,
  const std::vector<std::vector<geometrycentral::surface::SurfacePoint>> &curves
);
}
