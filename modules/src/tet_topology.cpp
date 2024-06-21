#include <modules/tet_topology.h>

#include <igl/tet_tet_adjacency.h>

#include <map>

#include <iostream>

namespace modules {
  std::tuple<
    Eigen::MatrixXi,
    Eigen::MatrixXi
  > get_e(
    const Eigen::MatrixXi &T
  ) {
    std::vector<std::array<int, 2>> T_EV;
    std::map<std::pair<int, int>, bool> edgeRegistered;
    std::map<std::pair<int, int>, int> edgeId;

    Eigen::MatrixXi TE(T.rows(), 6);
    TE.setConstant(-1);

    for (int i = 0; i < T.rows(); i++) {
      for (int j = 0; j < 4; j++) {
        for (int k = 0; k < 4; k++) {
          if (j == k) {
            continue;
          }
          
          int vid1 = T(i, j), vid2 = T(i, k);
          if (vid1 > vid2) {
            std::swap(vid1, vid2);
          }

          assert(vid1 < vid2);

          std::pair<int, int> edge = {vid1, vid2};

          if (!edgeRegistered[edge]) {
            T_EV.emplace_back(std::array<int, 2> {vid1, vid2});
            edgeRegistered[edge] = true;
            edgeId[edge] = T_EV.size() - 1;
          }
        }
      }

      // order is [0, 1] [1, 2] [2, 0] [0, 3] [1, 3] [2, 3]
      std::array<std::array<int, 2>, 6> edgeList = {
        std::array<int, 2> {0, 1},
        std::array<int, 2> {1, 2},
        std::array<int, 2> {2, 0},
        std::array<int, 2> {0, 3},
        std::array<int, 2> {1, 3},
        std::array<int, 2> {2, 3}
      };
      
      for (int j = 0; j < 6; j++) {
        int i0 = edgeList[j][0], i1 = edgeList[j][1];
        int v0 = T(i, i0), v1 = T(i, i1);

        if (v0 > v1) {
          std::swap(v0, v1);
        }

        std::pair<int, int> edge = {v0, v1};

        assert(edgeRegistered[edge]);

        TE(i, j) = edgeId[edge];
      }
    }

    Eigen::MatrixXi EV(T_EV.size(), 2);
    for (int i = 0; i < T_EV.size(); i++) {
      EV(i, 0) = T_EV[i][0];
      EV(i, 1) = T_EV[i][1];
    }

    return {EV, TE};
  }

  std::tuple<
    Eigen::MatrixXi,
    Eigen::MatrixXi
  > get_f(
    const Eigen::MatrixXi &T
  ) {
    Eigen::MatrixXi TT, TTi;
    igl::tet_tet_adjacency(T, TT, TTi);

    std::vector<std::array<int, 2>> ft, ftIndices;

    for (int i = 0; i < TT.rows(); i++) {
      for (int j = 0; j < TT.cols(); j++) {
        int tetId = i, otherTetId = TT(i, j);

        if (TT(i, j) == -1) {
          ft.emplace_back(std::array<int, 2> {tetId, -1});
          ftIndices.emplace_back(std::array<int, 2> {j, -1});
          continue;
        }

        assert(TT(TT(i, j), TTi(i, j)) == i);

        if (tetId < otherTetId) {
          ft.emplace_back(std::array<int, 2> {tetId, otherTetId});
          ftIndices.emplace_back(std::array<int, 2> {j, TTi(i, j)});
        }
      }
    }

    Eigen::MatrixXi FT(ft.size(), 2);
    for (int i = 0; i < ft.size(); i++) {
      FT(i, 0) = ft[i][0];
      FT(i, 1) = ft[i][1];
    }

    Eigen::MatrixXi TF(ft.size(), 4);
    TF.setConstant(-1);
    for (int i = 0; i < ft.size(); i++) {
      for (int j = 0; j < 2; j++) {
        int tetId = ft[i][j];
        int index = ftIndices[i][j];

        if (index == -1) {
          continue;
        }

        TF(tetId, index) = i;
      }
    }

    // std::cout << "facesnum: " << FT.rows() << std::endl;

    return {FT, TF};
  }

  Eigen::MatrixXi get_fe(
    const Eigen::MatrixXi &T,
    const Eigen::MatrixXi &TE,
    const Eigen::MatrixXi &TF,
    const int numFaces
  ) {
    Eigen::MatrixXi FE(numFaces, 3);
    FE.setConstant(-1);

    for (int i = 0; i < T.rows(); i++) {
      // in TE, order is 0: [0, 1], 1: [1, 2], 2: [2, 0], 3: [0, 3], 4: [1, 3], 5: [2, 3]
      // in TF, order is [0, 1, 2] [0, 1, 3] [1, 2, 3] [2, 0, 3]
      // thus, FE should be lik [0, 1, 2] [0, 4, 3] [1, 5, 4] [2, 3, 5]
      std::array<std::array<int, 3>, 4> feMapping = {
        std::array<int, 3> {0, 1, 2},
        std::array<int, 3> {0, 4, 3},
        std::array<int, 3> {1, 5, 4},
        std::array<int, 3> {2, 3, 5}
      };

      for (int j = 0; j < 4; j++) {
        int fid = TF(i, j);

        if (fid == -1) {
          continue;
        }

        assert(fid != -1);

        auto mapping = feMapping[j];

        for (int k = 0; k < 3; k++) {
          int eid = TE(i, mapping[k]);

          assert(eid != -1);

          FE(fid, k) = eid;
        }
      }
    }

    return FE;
  }

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
    // order is [0, 1] [1, 2] [2, 0] [0, 3] [1, 3] [2, 3]
    Eigen::MatrixXi& TE,
    // #t by 4 matrix of face indices
    // order is [0, 1, 2] [0, 1, 3] [1, 2, 3] [2, 0, 3]
    Eigen::MatrixXi& TF
  ) {
    // // 1. get list of EV
    std::tie(EV, TE) = get_e(T);

    // // 2. get list of F
    std::tie(FT, TF) = get_f(T);

    // // 3. get list of FE
    FE = get_fe(T, TE, TF, FT.rows());
  }
}
