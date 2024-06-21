#include "modules/cut_mesh_with_curve.h"

#include <igl/boundary_loop.h>
#include <igl/triangle/cdt.h>
#include <igl/is_edge_manifold.h>
#include <igl/edge_topology.h>
#include <igl/boundary_loop.h>

#include "modules/surface_point_to_cartesian.h"
#include "modules/gc_igl_convert.h"

#include <queue>

namespace modules {
  std::tuple<
    std::vector<SurfacePoint>,
    std::vector<std::vector<std::array<int, 2>>>
  > construct_segments_on_face(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry,
    const std::vector<SurfacePoint> &nodes,
    const std::vector<std::array<int, 2>> &segments,
    const std::vector<std::vector<SurfacePoint>> &segmentSurfacePoints
  ) {
    std::vector<SurfacePoint> surfacePoints = {};
    std::vector<std::vector<std::array<int, 2>>> segmentsOnFace(mesh.nFaces());

    for (int i = 0; i < nodes.size(); i++) {
      auto sp = nodes[i];
      surfacePoints.push_back(sp);
    }

    for (int i = 0; i < segments.size(); i++) {
      int node0 = segments[i][0];
      int node1 = segments[i][1];

      SurfacePoint surfacePoint0 = nodes[node0];
      SurfacePoint surfacePoint1 = nodes[node1];

      std::vector<SurfacePoint> sps = {surfacePoint0};
      std::vector<int> spsIndex = {node0};

      for (auto sp : segmentSurfacePoints[i]) {
        sps.push_back(sp);

        surfacePoints.push_back(sp);
        spsIndex.push_back(surfacePoints.size() - 1);
      }
      sps.push_back(surfacePoint1);
      spsIndex.push_back(node1);

      for (int j = 0; j < sps.size() - 1; j++) {
        auto prev = sps[j];
        auto next = sps[j + 1];

        int prevIndex = spsIndex[j];
        int nextIndex = spsIndex[j + 1];

        std::array<int, 2> segment = {prevIndex, nextIndex};

        if (
          prev.type == SurfacePointType::Vertex &&
          next.type == SurfacePointType::Vertex
        ) {
          continue;
        }

        if (
          (
            prev.type == SurfacePointType::Vertex ||
            next.type == SurfacePointType::Vertex
          ) && (
            prev.type == SurfacePointType::Edge ||
            next.type == SurfacePointType::Edge
          )
        ) {
          auto edge = prev.type == SurfacePointType::Edge ? prev.edge : next.edge;
          auto vertex = prev.type == SurfacePointType::Vertex ? prev.vertex : next.vertex;

          if (edge.firstVertex() == vertex || edge.secondVertex() == vertex) {
            continue;
          }

          Face face;
          for (auto f : edge.adjacentFaces()) {
            for (auto v : f.adjacentVertices()) {
              if (v == vertex) {
                face = f;
                break;
              }
            }
          }

          assert(face != Face());

          int faceId = face.getIndex();

          segmentsOnFace[faceId].push_back(segment);
        } else if (
          prev.type == SurfacePointType::Face &&
          next.type == SurfacePointType::Face
        ) {
          assert(prev.face == next.face);

          int faceId = prev.face.getIndex();

          segmentsOnFace[faceId].push_back(segment);
        } else if (
          prev.type == SurfacePointType::Edge &&
          next.type == SurfacePointType::Edge
        ) {
          auto e0 = prev.edge;
          auto e1 = next.edge;

          Face f;
          for (auto f0 : e0.adjacentFaces()) {
            for (auto f1 : e1.adjacentFaces()) {
              if (f0 == f1) {
                f = f0;
                break;
              }
            }
          }
          assert(f != Face());

          int faceId = f.getIndex();

          segmentsOnFace[faceId].push_back(segment);
        } else if (
          prev.type == SurfacePointType::Edge ||
          next.type == SurfacePointType::Edge
        ) {
          int faceId = prev.type == SurfacePointType::Face ? prev.face.getIndex() : next.face.getIndex();

          segmentsOnFace[faceId].push_back(segment);

          {
            auto edge = prev.type == SurfacePointType::Edge ? prev.edge : next.edge;
            auto face = prev.type == SurfacePointType::Face ? prev.face : next.face;
            bool isEdgeInFace = false;
            for (auto e : face.adjacentEdges()) {
              if (e == edge) {
                isEdgeInFace = true;
                break;
              }
            }
            assert(isEdgeInFace);
          }
        } else if (
          prev.type == SurfacePointType::Vertex ||
          next.type == SurfacePointType::Vertex
        ) {
          int faceId = prev.type == SurfacePointType::Face ? prev.face.getIndex() : next.face.getIndex();

          segmentsOnFace[faceId].push_back(segment);

          {
            auto vertex = prev.type == SurfacePointType::Vertex ? prev.vertex : next.vertex;
            auto face = prev.type == SurfacePointType::Face ? prev.face : next.face;
            bool isVertexInFace = false;
            for (auto v : face.adjacentVertices()) {
              if (v == vertex) {
                isVertexInFace = true;
                break;
              }
            }
            assert(isVertexInFace);
          }
        }
      }
    }

    return std::make_tuple(surfacePoints, segmentsOnFace);
  }

  void vector_decomposition_3d(Eigen::Vector3d &x, Eigen::Vector3d &v_a,
                             Eigen::Vector3d &v_b,
                             Eigen::Vector2d &result) {
    Eigen::Matrix2d matrix;
    matrix << v_a(0), v_b(0), v_a(1), v_b(1);
    Eigen::Vector2d _x;
    _x << x(0), x(1);
    result = matrix.inverse() * _x;
  };

  void subdivide_by_path(const Eigen::MatrixXd &V, 
                       const Eigen::MatrixXi &F, 
                       const std::vector<std::vector<std::vector<int>>> &pathOnFaces,
                       Eigen::MatrixXi &NF,
                       std::vector<std::vector<int>> &F2NF,
                       std::vector<int> &NF2F) {
    std::vector<std::array<int, 3>> newFacesList = {};

    assert(igl::is_edge_manifold(F));

    F2NF.resize(F.rows(), std::vector<int>());
    NF2F = {};

    for (int i = 0; i < pathOnFaces.size(); i++) {
      auto paths = pathOnFaces[i];
      auto pathsNum = paths.size();
      int fid = i;

      if (paths.size() == 0) {
        std::array<int, 3> f = {F(fid, 0), F(fid, 1), F(fid, 2)};
        newFacesList.emplace_back(f);
        NF2F.emplace_back(fid);

        int nfid = newFacesList.size() - 1;
        F2NF[fid].emplace_back(nfid);

        continue;
      }

      std::vector<int> vertexList = {};
      std::map<int, bool> isVertexInList;
      std::map<int, int> v2nv;
      std::vector<std::pair<int, int>> edgesList = {};

      for (int j = 0; j < 3; j++) {
        int vid = F(fid, j);
        vertexList.emplace_back(vid);
        isVertexInList[vid] = true;
        v2nv[vid] = vertexList.size() - 1;
      }

      // edge 0 is v0->v1, edge 1 is v1->v2, edge 2 is v2->v0
      std::array<std::vector<std::pair<int, double>>, 3> pointsOnEdges = {
        std::vector<std::pair<int, double>>(),
        std::vector<std::pair<int, double>>(),
        std::vector<std::pair<int, double>>()
      };

      for (auto path : paths) {
        int prev = -1;
        for (int j = 0; j < path.size(); j++) {
          auto p = path[j];

          if (!isVertexInList[p]) {
            vertexList.emplace_back(p);
            isVertexInList[p] = true;
            v2nv[p] = vertexList.size() - 1;
          }

          if (j > 0) {
            int curr = v2nv[p];
            edgesList.emplace_back(std::pair<int, int>(curr, prev));
          }

          prev = v2nv[p];
        }

        // for the first and last element, constrain the face to outside
        for (int j = 0; j < path.size(); j++) {
          if (j != 0 && j != path.size() - 1) {
            continue;
          }

          auto p = path[j];

          for (int k = 0; k < 3; k++) {
            int e0 = vertexList[k], e1 = vertexList[(k + 1) % 3];
            double d0 = (V.row(e0) - V.row(p)).norm();
            double d1 = (V.row(e1) - V.row(p)).norm();
            double d2 = (V.row(e1) - V.row(e0)).norm();

            if (std::abs(d2 - d1 - d0) < .00001) {
              double ratio = d0 / d2;
              pointsOnEdges[k].emplace_back(std::pair<int, double>(v2nv[p], ratio));
              break;
            }
          }
        }
      }

      // re-order the vertices on the edges
      for (int i = 0; i < 3; i++) {
        std::sort(pointsOnEdges[i].begin(), pointsOnEdges[i].end(), [&](std::pair<int, double> a, std::pair<int, double> b) {
          return a.second < b.second;
        });
      }

      // add all the edges
      for (int i = 0; i < 3; i++) {
        std::vector<int> pointsOnEdge = {i};

        for (int j = 0; j < pointsOnEdges[i].size(); j++) {
          pointsOnEdge.emplace_back(pointsOnEdges[i][j].first);
        }

        pointsOnEdge.emplace_back((i + 1) % 3);
        
        for (int j = 1; j < pointsOnEdge.size(); j++) {
          edgesList.emplace_back(pointsOnEdge[j - 1], pointsOnEdge[j]);
        }
      }

      Eigen::Vector3d v0 = V.row(F(fid, 0)).transpose();
      Eigen::Vector3d v_a = (V.row(F(fid, 1)).transpose() - v0).normalized();
      Eigen::Vector3d v_b = (V.row(F(fid, 2)).transpose() - v0).normalized();
      Eigen::Vector3d v_n = v_a.cross(v_b).normalized();
      Eigen::Vector3d v_n2 = v_n.cross(v_a).normalized();

      Eigen::MatrixXd vertices(vertexList.size(), 2);
      for (int i = 0; i < vertexList.size(); i++) {
        // project the vertex onto the face and get the xy coordinate
        auto vid = vertexList[i];

        Eigen::Vector3d x = V.row(vid).transpose() - v0;

        Eigen::Vector2d result;
        vector_decomposition_3d(x, v_a, v_n2, result);

        vertices.row(i) = result.transpose();
      }

      Eigen::MatrixXd edges(edgesList.size(), 2);
      for (int j = 0; j < edgesList.size(); j++) {
        auto e = edgesList[j];
        edges(j, 0) = std::min(e.first, e.second);
        edges(j, 1) = std::max(e.first, e.second);
      }

      // std::cout << "vertices: " << vertices << std::endl;
      // std::cout << std::endl;
      // std::cout << "edges: " << edges << std::endl;
      // std::cout << std::endl;

      Eigen::MatrixXd newV;
      Eigen::MatrixXi newF;
      igl::triangle::triangulate(vertices, edges, Eigen::MatrixXi(0, 2), "-YYQ", newV, newF);

      // std::cout << "newV: " << newV << std::endl;
      // std::cout << std::endl;
      // std::cout << "newF: " << newF << std::endl;
      // std::cout << std::endl;

      assert(newV.rows() == vertices.rows());

      for (int j = 0; j < newF.rows(); j++) {
        std::array<int, 3> f = {vertexList[newF(j, 0)], vertexList[newF(j, 1)], vertexList[newF(j, 2)]};
        newFacesList.emplace_back(f);
        NF2F.emplace_back(fid);

        int nfid = newFacesList.size() - 1;
        F2NF[fid].emplace_back(nfid);

        for (int k = 0; k < 3; k++) {
          auto v0 = newV.row(newF(j, k)), v1 = newV.row(newF(j, (k + 1) % 3)), v2 = newV.row(newF(j, (k + 2) % 3));

          if (std::abs((v1 - v0).norm() + (v2 - v0).norm() - (v1 - v2).norm()) < .00001) {
            std::cout << "warn! the system produced thin triangles: " << j << ", " << k << ": " << (v1 - v0).norm() << " " << (v2 - v0).norm() << " " << (v1 - v2).norm() << std::endl;
          }
        }
      }

      // for debug
      {
        // NF.resize(newFacesList.size(), 3);

        // for (int i = 0; i < newFacesList.size(); i++) {
        //   for (int j = 0; j < 3; j++) {
        //     NF(i, j) = newFacesList[i][j];
        //   }
        // }

        // assert(igl::is_edge_manifold(NF));
      }
    }

    NF.resize(newFacesList.size(), 3);

    for (int i = 0; i < newFacesList.size(); i++) {
      for (int j = 0; j < 3; j++) {
        NF(i, j) = newFacesList[i][j];
      }
    }

    assert(igl::is_edge_manifold(NF));
  };

  std::tuple<
    Eigen::MatrixXd,
    Eigen::MatrixXi
  > cut_faces(
    const std::vector<SurfacePoint> &surfacePoints,
    const std::vector<std::vector<std::array<int, 2>>> &segmentsOnFace,
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry
  ) {
    // 1. add surfacePoints to the mesh and geometry
    auto cartesianCoords = modules::surface_point_to_cartesian(mesh, geometry, surfacePoints);

    std::map<int, int> sp2vid;
    std::vector<Vector3> newVertices = {};

    int count = mesh.nVertices();
    for (int i = 0; i < surfacePoints.size(); i++) {
      if (
        surfacePoints[i].type == SurfacePointType::Vertex
      ) {
        sp2vid[i] = surfacePoints[i].vertex.getIndex();
      } else {
        sp2vid[i] = count;

        auto coord = cartesianCoords[i];
        newVertices.emplace_back(coord);

        count++;
      }
    }

    std::vector<std::vector<std::vector<int>>> pathOnFaces(mesh.nFaces());

    for (int i = 0; i < segmentsOnFace.size(); i++) {
      for (int j = 0; j < segmentsOnFace[i].size(); j++) {
        std::vector<int> path = {
          sp2vid[segmentsOnFace[i][j][0]],
          sp2vid[segmentsOnFace[i][j][1]] 
        };
        pathOnFaces[i].emplace_back(path);
      }
    }

    auto [V, F, EF, FE, EV] = gc_to_igl(mesh, geometry);
    V.conservativeResize(V.rows() + newVertices.size(), 3);
    for (int i = 0; i < newVertices.size(); i++) {
      for (int j = 0; j < 3; j++) {
        int id = V.rows() - newVertices.size() + i;
        assert(id >= mesh.nVertices());
        V(id, j) = newVertices[i][j];
      }
    }

    std::vector<std::vector<int>> _;
    std::vector<int> __;
    Eigen::MatrixXi NF;
    subdivide_by_path(V, F, pathOnFaces, NF, _, __);
    F = NF;

    std::vector<std::vector<int>> lps;
    igl::boundary_loop(F, lps);

    if (lps.size() != mesh.nBoundaryLoops() && false) {
      // fill up the triangle
      // this is heuristics, sorry
      int loopnum = 0;
      for (auto l : lps) {
        if (l.size() == 3) {
          loopnum++;
        }
      }

      F.conservativeResize(F.rows() + loopnum, 3);
      int cnt = 0;
      for (int i = 0; i < lps.size(); i++) {
        if (lps[i].size() != 3) {
          continue;
        }
        assert(lps[i].size() == 3);
        for (int j = 0; j < lps[i].size(); j++) {
          F(F.rows() - lps.size() + cnt, 2 - j) = lps[i][j];
        }
        cnt++;
      }

      assert(cnt == loopnum);

      std::cout << "warn! the topology of the mesh is not preserved!" << std::endl;
      std::cout << "current numBoundaryLoops: " << lps.size() << std::endl;
      for (auto l : lps) {
        for (auto p : l) {
          std::cout << p << " ";
        }
        std::cout << std::endl;
      }
      std::cout << "previous numBoundaryLoops: " << mesh.nBoundaryLoops() << std::endl << std::endl;
    }
    // assert(_.size() == mesh.nBoundaryLoops());

    assert(igl::is_edge_manifold(F));

    return {V, F};
  }

  std::tuple<
    std::vector<Eigen::MatrixXd>,
    std::vector<Eigen::MatrixXi>
  > segment_faces(
    Eigen::MatrixXd &V,
    Eigen::MatrixXi &F
  ) {
    std::vector<Eigen::MatrixXd> Vs = {};
    std::vector<Eigen::MatrixXi> Fs = {};

    Eigen::MatrixXi EF, FE, EV;
    igl::edge_topology(V, F, EV, FE, EF);

    std::queue<int> q;
    std::map<int, bool> isFaceVisited;

    int itr = 0;

    while (isFaceVisited.size() < F.rows()) {
      itr++;

      // std::cout << std::endl << "itr: " << itr << std::endl;

      // assert(itr < 3);

      int f0 = -1;
      for (int i = 0; i < F.rows(); i++) {
        if (isFaceVisited[i]) {
          continue;
        }

        f0 = i;
        isFaceVisited[i] = true;
        // std::cout << f0 << std::endl;
        break;
      }

      assert(f0 != -1);

      for (int j = 0; j < 3; j++) {
        int e = FE(f0, j);

        assert(e != -1);

        // int v0 = EV(e, 0), v1 = EV(e, 1);
        // if (v0 >= maxVid && v1 >= maxVid) {
        //   continue;
        // }

        int otherf = EF(e, 0) == f0 ? EF(e, 1) : EF(e, 0);

        q.push(otherf);
      }

      std::vector<int> faceList = {f0};

      while (!q.empty()) {
        int f = q.front();
        q.pop();

        // std::cout << f << ", " << isFaceVisited[f] << std::endl;

        if (f == -1) {
          continue;
        }

        if (isFaceVisited[f]) {
          continue;
        }

        isFaceVisited[f] = true;
        faceList.emplace_back(f);

        for (int j = 0; j < 3; j++) {
          int e = FE(f, j);

          assert(e != -1);

          // int v0 = EV(e, 0), v1 = EV(e, 1);s
          // if (v0 >= maxVid && v1 >= maxVid) {
          //   continue;
          // }

          int otherf = EF(e, 0) == f ? EF(e, 1) : EF(e, 0);

          q.push(otherf);
        }
      }

      std::map<int, bool> isVertexVisited;
      std::map<int, int> v2nv;
      for (auto f : faceList) {
        for (int j = 0; j < 3; j++) {
          int v = F(f, j);
          if (isVertexVisited[v]) {
            continue;
          }

          isVertexVisited[v] = true;
          int nvid = v2nv.size();
          v2nv[v] = nvid;
        }
      }

      Eigen::MatrixXi _F(faceList.size(), 3);

      for (int i = 0; i < faceList.size(); i++) {
        for (int j = 0; j < 3; j++) {
          _F(i, j) = F(faceList[i], j);
        }
      }

      Vs.emplace_back(V);
      Fs.emplace_back(_F);
    }

    return {Vs, Fs};
  }

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
  ) {
    auto [surfacePoints, segmentsOnFace] = construct_segments_on_face(_mesh, _geometry, _nodes, _segments, _segmentSurfacePoints);

    // auto [V, F] = cut_faces(surfacePoints, segmentsOnFace, _mesh, _geometry);

    // cut_faces doesnt work, so I'll go the other way around
    std::map<int, bool> isDelete = {};
    for (int i = 0; i < _mesh.nFaces(); i++) {
      if (segmentsOnFace[i].size() > 0) {
        isDelete[i] = true;
      }
    }
    auto [V, F, EF, FE, EV] = gc_to_igl(_mesh, _geometry);
    Eigen::MatrixXi NF(F.rows() - isDelete.size(), 3);
    int cnt = 0;
    for (int i = 0; i < F.rows(); i++) {
      if (!isDelete[i]) {
        NF.row(cnt) = F.row(i);
        cnt++;
      }
    }
    F = NF;

    auto [Vs, Fs] = segment_faces(V, F);

    return {Vs, Fs, V, F};
  }
}
