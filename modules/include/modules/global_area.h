#pragma once

#include <Eigen/Core>
#include <vector>

namespace modules {
std::vector<int> extract_global_area_3d(
  const Eigen::MatrixXi &T,
  const Eigen::MatrixXi &T_FE,
  const std::vector<std::vector<int>> &T_EF,
  const std::vector<std::vector<int>> &T_ET,
  const std::vector<int> &innerEdges,
  const std::vector<int> &innerFaces,
  const std::vector<double> &radiusPerTet,
  const double _radius
);
}
