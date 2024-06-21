#pragma once

#include <Eigen/Core>
#include <igl/tet_tet_adjacency.h>

namespace modules {
  void tet_topology(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& T,
    // #E by 2 matrix of vertex indices
    Eigen::MatrixXi& EV,
    // #F by 3 matrix of edge indices
    Eigen::MatrixXi& FE,
    // #f by 2 matrix of tet indices
    Eigen::MatrixXi& FT,
    // #t by 6 matrix of edge indices
    Eigen::MatrixXi& TE,
    // #t by 4 matrix of face indices
    Eigen::MatrixXi& TF
  );
}
