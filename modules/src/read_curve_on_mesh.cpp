#include "modules/read_curve_on_mesh.h"

#include "modules/scene_file.h"

#include "igl/point_mesh_squared_distance.h"
#include "igl/barycentric_coordinates.h"
#include "igl/AABB.h"

const double EPSILON = 1e-5;

namespace modules {
std::tuple<
  std::vector<SurfacePoint>,
  std::vector<std::array<int, 2>>,
  std::vector<bool>
> readCurveOnMesh(
  const std::string filename,
  ManifoldSurfaceMesh &mesh,
  VertexPositionGeometry &geometry
) {
  Eigen::MatrixXd V;
  V.resize(mesh.nVertices(), 3);
  for (int i = 0; i < mesh.nVertices(); i++) {
    auto p = geometry.vertexPositions[i];
    V.row(i) << p.x, p.y, p.z;
  }

  Eigen::MatrixXi F;
  F.resize(mesh.nFaces(), 3);
  for (int i = 0; i < mesh.nFaces(); i++) {
    auto f = mesh.face(i);
    auto vs = f.adjacentVertices();

    int j = 0;
    for (auto v : vs) {
      F(i, j) = v.getIndex();
      j++;
    }
  }

  igl::AABB<Eigen::MatrixXd, 3> tree;
  tree.init(V,F);

  std::vector<SurfacePoint> nodes = {};
  std::vector<std::array<int, 2>> segments = {};
  std::vector<bool> isFixedNode = {};

  std::ifstream inFile;
  inFile.open(filename);

  if (!inFile) {
    std::cerr << "Could not open file " << filename << std::endl;
    exit(1);
  }

  std::vector<Vector3> addedNodes = {};
  std::map<int, bool> isAddedNodeFixed = {};
  std::vector<std::array<int, 2>> addedSegments = {};

  for (std::string line; std::getline(inFile, line ); ) {
    std::vector<std::string> parts;
    modules::splitString(line, parts, ' ');

    if (parts.size() == 0) continue;

    if (parts[0] == "v" && parts.size() >= 4) {
      auto coord = Vector3 {
        std::stof(parts[1]),
        std::stof(parts[2]),
        std::stof(parts[3])
      };

      int vid = addedNodes.size();
      addedNodes.push_back(coord);
    } else if (parts[0] == "l") {
      for (int i = 2; i < parts.size(); i++) {
        int v0 = (std::stoi(parts[i-1]) - 1);
        int v1 = (std::stoi(parts[i]) - 1);
        std::array<int, 2> line = {v0, v1};
        addedSegments.push_back(line);
      }
    } else if (parts[0] == "fixed") {
      for (int i = 1; i < parts.size(); i++) {
        int id = (std::stoi(parts[i]) - 1);
        isAddedNodeFixed[id] = true;
      }
    }
  }

  std::map<int, bool> nodeAdded = {};
  std::map<int, int> addedNode2Index = {};
  std::vector<int> addedNodeIds = {};

  for (int i = 0; i < addedSegments.size(); i++) {
    auto seg = addedSegments[i];
    std::array<int, 2> newSeg = {-1, -1};
    for (int j = 0; j < 2; j++) {
      if (!nodeAdded[seg[j]]) {
        auto vid = addedNodeIds.size();
        addedNodeIds.push_back(seg[j]);
        isFixedNode.emplace_back(isAddedNodeFixed[seg[j]]);

        nodeAdded[seg[j]] = true;
        addedNode2Index[seg[j]] = vid;
      }

      newSeg[j] = addedNode2Index[seg[j]];
    }

    segments.emplace_back(newSeg);
  }

  Eigen::MatrixXd P(addedNodeIds.size(), 3);
  for (int i = 0; i < addedNodeIds.size(); i++) {
    auto p = addedNodes[addedNodeIds[i]];
    P.row(i) << p.x, p.y, p.z;
  }

  Eigen::MatrixXd B;
  Eigen::VectorXi I;

  {
    Eigen::VectorXd sqrD;
    Eigen::MatrixXd C;

    tree.squared_distance(V, F, P, sqrD, I, C);

    Eigen::MatrixXd _A(I.rows(), 3), _B(I.rows(), 3), _C(I.rows(), 3);

    for (int i = 0; i < I.rows(); i++) {
      _A.row(i) = V.row(F(I(i), 0));
      _B.row(i) = V.row(F(I(i), 1));
      _C.row(i) = V.row(F(I(i), 2));
    }

    igl::barycentric_coordinates(C, _A, _B, _C, B);
  }

  // std::cout << B << std::endl;
  
  for (int i = 0; i < I.rows(); i++) {
    auto fid = I(i);
    Eigen::Vector3d _bcoord = B.row(i).transpose();
    Vector3 bcoord {_bcoord(0), _bcoord(1), _bcoord(2)};

    assert(_bcoord.sum() >= 1 - EPSILON);
    assert(_bcoord.minCoeff() >= -EPSILON);

    int maxIdx, minIdx;
    _bcoord.maxCoeff(&maxIdx);
    _bcoord.minCoeff(&minIdx);

    if (_bcoord(maxIdx) >= 1 - EPSILON) {
      // on a vertex
      nodes.emplace_back(SurfacePoint(mesh.vertex(F(fid, maxIdx))));
    } else if (_bcoord(minIdx) <= EPSILON) {
      // on an edge
      int v0 = F(fid, (maxIdx + 1) % 3), v1 = F(fid, (maxIdx + 2) % 3);
      double t0 = _bcoord((maxIdx + 1) % 3), t1 = _bcoord((maxIdx + 2) % 3);

      for (auto e : mesh.face(fid).adjacentEdges()) {
        if (e.halfedge().vertex().getIndex() == v0) {
          nodes.emplace_back(SurfacePoint(e, t0));
          break;
        } else if (e.halfedge().vertex().getIndex() == v1) {
          nodes.emplace_back(SurfacePoint(e, t1));
          break;
        }
      }
    } else {
      // on a face
      nodes.emplace_back(SurfacePoint(mesh.face(fid), bcoord));
    }
  }

  return {nodes, segments, isFixedNode};
}
}
