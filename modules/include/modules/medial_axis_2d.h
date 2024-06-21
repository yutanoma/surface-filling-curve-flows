#pragma once

#include <Eigen/Core>

#include <igl/predicates/delaunay_triangulation.h>
#include <igl/triangle/cdt.h>
#include <igl/edge_topology.h>

#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;

namespace modules {
  std::tuple<
    std::vector<Vector2>, // node list on medial axis
    std::vector<std::array<int, 2>>, // segment list on medial axis
    std::vector<double>, // radius list on medial axis
    std::vector<std::vector<std::vector<int>>>, // node list to the nearest medial axis node
    Eigen::MatrixXd, // V of underlying mesh
    Eigen::MatrixXi, // F of underlying mesh
    Eigen::MatrixXi // V2curveNode
  > get_medial_axis_2d(
    const std::vector<std::vector<Vector2>> &curves,
    const int &mode
  );
}
