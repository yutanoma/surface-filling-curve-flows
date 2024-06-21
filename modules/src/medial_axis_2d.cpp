#include <modules/medial_axis_2d.h>

#include <modules/circumcenter.h>
#include <modules/evolution_mode.h>

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
) {
  int numVertices = 0;
  std::vector<std::vector<int>> curveNodeList(curves.size(), std::vector<int> {});
  for (int i = 0; i < curves.size(); i++) {
    auto curve = curves[i];
    numVertices += curve.size();
    curveNodeList[i].resize(curve.size(), -1);
  }

  // 1. do delaunay triangulation
  Eigen::MatrixXi V2curveNode(numVertices, 2), _E(numVertices, 2);
  Eigen::MatrixXd V(numVertices, 2), _V;

  int _numVertices = 0;
  for (int i = 0; i < curves.size(); i++) {
    int firstIdx = _numVertices;

    for (int j = 0; j < curves[i].size(); j++) {
      V(_numVertices, 0) = curves[i][j].x;
      V(_numVertices, 1) = curves[i][j].y;
      V2curveNode(_numVertices, 0) = i;
      V2curveNode(_numVertices, 1) = j;
      curveNodeList[i][j] = _numVertices;

      _E(_numVertices, 0) = _numVertices;
      _E(_numVertices, 1) = j == curves[i].size() - 1
        ? firstIdx
        : _numVertices + 1;

      _numVertices++;
    }
  }

  // std::cout << _E << std::endl;

  Eigen::MatrixXi F, J;
  igl::triangle::triangulate(V, _E, Eigen::MatrixXi(0, 2), "cYS0Q", _V, F);

  assert(V.rows() == _V.rows());

  V = _V;

  // igl::predicates::delaunay_triangulation(V, F);
  // igl::copyleft::cgal::delaunay_triangulation(V, F);

  // 2. take the dual and get the voronoi diagram
  Eigen::MatrixXi EV, FE, EF;
  igl::edge_topology(V, F, EV, FE, EF);

  Eigen::VectorXi isFInside(F.rows());
  isFInside.setConstant(-1);

  std::vector<std::vector<int>> VE(V.rows(), std::vector<int> {});
  for (int i = 0; i < EV.rows(); i++) {
    int v0 = EV(i, 0), v1 = EV(i, 1);

    VE[v0].push_back(i);
    VE[v1].push_back(i);
  }

  for (int i = 0; i < curves.size(); i++) {
    for (int j = 0; j < curves[i].size(); j++) {
      int prev = (j - 1 + curves[i].size()) % curves[i].size();
      int curr = j;

      int v_prev = curveNodeList[i][prev];
      int v_curr = curveNodeList[i][curr];

      // get two faces
      int e = -1;
      for (int k = 0; k < VE[v_prev].size(); k++) {
        if (EV(VE[v_prev][k], 0) == v_curr || EV(VE[v_prev][k], 1) == v_curr) {
          e = VE[v_prev][k];
          break;
        }
      }

      if (e == -1) {
        std::cout << "edge flipped in curve[" << i << "][" << prev << " -> " << curr << "]" << std::endl;
        continue;
      }

      // assert(e != -1);

      for (int k = 0; k < 2; k++) {
        int f = EF(e, k);

        if (f == -1) {
          continue;
        }

        for (int l = 0; l < 3; l++) {
          if (F(f, l) == v_prev && F(f, (l+1)%3) == v_curr) {
            isFInside(f) = Mode::Outside;
            break;
          } else if (F(f, l) == v_curr && F(f, (l+1)%3) == v_prev) {
            isFInside(f) = Mode::Inside;
            break;
          }
        }
      }
    }
  }

  for (int i = 0; i < EF.rows(); i++) {
    if (EF(i, 0) == -1 || EF(i, 1) == -1) {
      continue;
    }

    if (isFInside(EF(i, 0)) == -1 && isFInside(EF(i, 1)) == -1) {
      continue;
    }

    assert(isFInside(EF(i, 0)) != -1 || isFInside(EF(i, 1)) != -1);

    if (isFInside(EF(i, 0)) == -1) {
      isFInside(EF(i, 0)) = isFInside(EF(i, 1));
    }

    if (isFInside(EF(i, 1)) == -1) {
      isFInside(EF(i, 1)) = isFInside(EF(i, 0));
    }
  }

  std::vector<Vector2> nodes = {};
  std::vector<std::array<int, 2>> segments = {};
  std::vector<double> ballRadii = {};

  Eigen::VectorXi F2Node(F.rows());
  F2Node.setConstant(-1);

  std::vector<std::vector<std::vector<int>>> nearestNodeList(curves.size(), std::vector<std::vector<int>> {});

  for (int i = 0; i < curves.size(); i++) {
    nearestNodeList[i].resize(curves[i].size(), std::vector<int> {});
  }

  // 3. get the medial axis
  for (int i = 0; i < F.rows(); i++) {
    if (isFInside(i) != mode && mode != Mode::Both) {
      continue;
    }

    int nodeId = nodes.size();
    F2Node(i) = nodeId;

    // compute the circumcenter of the triangle
    Vector2 circumCenter = modules::circumcenter(
      Vector2 {V(F(i, 0), 0), V(F(i, 0), 1)},
      Vector2 {V(F(i, 1), 0), V(F(i, 1), 1)},
      Vector2 {V(F(i, 2), 0), V(F(i, 2), 1)}
    );

    nodes.emplace_back(Vector2 {circumCenter.x, circumCenter.y} );

    // get the ball radii
    // ball radii is the circumCenter
    auto p = V.row(F(i, 0));
    double ballRadius = (circumCenter - Vector2 {p(0), p(1)}).norm();

    ballRadii.emplace_back(ballRadius);

    for (int j = 0; j < 3; j++) {
      int vid = F(i, j);

      int k = V2curveNode(vid, 0), l = V2curveNode(vid, 1);

      nearestNodeList[k][l].push_back(nodeId);
    }
  }

  for (int i = 0; i < EF.rows(); i++) {
    if (EF(i, 0) == -1 || EF(i, 1) == -1) {
      continue;
    }

    if (mode == Mode::Both && isFInside(EF(i, 0)) != isFInside(EF(i, 1)) || isFInside(EF(i, 0)) == -1) {
      continue;
    }

    if ((isFInside(EF(i, 0)) != mode || isFInside(EF(i, 1)) != mode) && mode != Mode::Both) {
      continue;
    }

    segments.emplace_back(std::array<int, 2> {
      F2Node(EF(i, 0)),
      F2Node(EF(i, 1))
    });
  }

  return {nodes, segments, ballRadii, nearestNodeList, V, F, V2curveNode};
}
}
