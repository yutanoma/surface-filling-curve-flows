#pragma once

#include <geometrycentral/surface/surface_point.h>
#include <geometrycentral/surface/intrinsic_geometry_interface.h>
#include <geometrycentral/surface/vertex_position_geometry.h>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
  std::tuple<
    Eigen::MatrixXd, // V
    Eigen::MatrixXi, // F
    Eigen::MatrixXi, // EF
    Eigen::MatrixXi, // FE
    Eigen::MatrixXi // EV
  > gc_to_igl(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry
  );

  std::tuple<
    Eigen::MatrixXd, // V
    Eigen::MatrixXi // F
  > gc_to_igl_vf(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry
  );
}
