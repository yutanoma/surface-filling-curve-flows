#include "modules/medial_axis.h"

#include <geometrycentral/surface/trace_geodesic.h>

#include <modules/exact_geodesics.h>
#include <modules/surface_point_to_cartesian.h>
#include <modules/get_tangent_basis.h>

#include <knn-cpp/knncpp.h>

namespace modules {
  int closestPointIndex(
    const Vector3 &p,
    const knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> &kdtree
  ) {
    Eigen::MatrixXd P(3, 1);
    Eigen::MatrixXd distances;
    knncpp::Matrixi indices;
    P << p.x, p.y, p.z;
    kdtree.query(P, 1, indices, distances);

    assert(indices.rows() == 1);

    return indices(0, 0);
  }

  double maximumBallRadius(
    const Vector3 &x,
    const Vector3 &b,
    const int i,
    const double maxRadius,
    const std::vector<Vector3> &cartesianCoords,
    knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> &kdtree
  ) {
    auto r = maxRadius;
    auto c = x + r * b;

    int nn = closestPointIndex(c, kdtree);
    bool finished = nn == i;

    double bsMax = 1., bsMin = 0.;
    int itrc = 0;

    // continue until the 
    while (!finished) {
      // std::cout << "  " << bsMax << ", " << bsMin << ", " << (bsMax + bsMin) / 2 << std::endl;

      itrc++;

      r = maxRadius * (bsMax + bsMin) / 2;

      auto c = x + r * b;
      int nn = closestPointIndex(c, kdtree);

      if (nn == i) {
        bsMin = (bsMax + bsMin) / 2;
      } else {
        // check if there is a point inside the ball
        Vector3 xy = cartesianCoords[nn] - cartesianCoords[i];
        r = pow(norm(xy), 2) / (2 * dot(xy, b));

        auto c = x + r * b;
        int _nn = closestPointIndex(c, kdtree);

        if (_nn == nn || _nn == i) {
          // if no one inside, break
          finished = true;
        } else {
          // else, continue the binary search
          bsMax = (bsMax + bsMin) / 2;
          assert(bsMax > bsMin);
        }
      }

      if (itrc > 100) {
        break;
      }
    }

    // std::cout << itrc << std::endl;

    return r;
  }

  std::vector<SurfacePoint> tracePath(
    VertexPositionGeometry &geometry,
    const SurfacePoint &node,
    const Vector3 &d,
    const double length
  ) {
    auto [x, y] = modules::get_tangent_basis(geometry, node);

    Vector2 traceVec {dot(d, x), dot(d, y)};
    TraceOptions opts = {
      includePath: true,
      maxIters: 1000
    };
    auto res = traceGeodesic(geometry, node, traceVec, opts);

    return res.pathPoints;
  }

  bool isIdenticalSurfacePoint(
    const SurfacePoint &p1,
    const SurfacePoint &p2
  ) {
    if (p1.type != p2.type) {
      return false;
    }

    switch (p1.type) {
      case SurfacePointType::Vertex:
        return p1.vertex == p2.vertex;
      case SurfacePointType::Edge:
        return p1.edge == p2.edge && std::abs(p1.tEdge - p2.tEdge) < 1e-6;
      case SurfacePointType::Face:
        return p1.face == p2.face && std::abs(norm(p1.faceCoords - p2.faceCoords)) < 1e-6;
      default:
        throw std::runtime_error("Invalid SurfacePointType");
    }
  }

  std::vector<SurfacePoint> cutTracePath(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &meshGeom,
    const std::vector<SurfacePoint> &path,
    const std::vector<Vector3> &cartesianPath,
    const double _length
  ) {
    double length = .0;

    assert(path.size() == cartesianPath.size());

    std::vector<SurfacePoint> tracePath;

    for (int j = 1; j < cartesianPath.size(); j++) {
      double currentEdgeLen = (cartesianPath[j] - cartesianPath[j - 1]).norm();
      length += currentEdgeLen;

      tracePath.emplace_back(path[j-1]);

      if (length > _length) {
        double ratio = (_length - (length - currentEdgeLen)) / currentEdgeLen;

        Face face = sharedFace(path[j - 1], path[j]);

        Vector3 vec0 = path[j - 1].inFace(face).faceCoords;
        Vector3 vec1 = path[j].inFace(face).faceCoords;

        auto sp = SurfacePoint(face, (vec0 * (1 - ratio) + vec1 * ratio));
        tracePath.emplace_back(sp);

        break;
      }

      if (j == path.size() - 1) {
        tracePath.emplace_back(path.back());
      }
    }

    return tracePath;
  }

  double maximumGeodesicBallRadius(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry,
    const SurfacePoint &x,
    const Vector3 &b,
    const int i,
    const double maxRadius,
    GeodesicAlgorithmExact &mmp
  ) {
    auto r = maxRadius;
    auto wholePath = tracePath(geometry, x, b, r);
    std::vector<Vector3> cartesianCoords = modules::surface_point_to_cartesian(mesh, geometry, wholePath);

    auto p = mmp.traceBack(wholePath.back());
    bool finished = isIdenticalSurfacePoint(p.back(), x);

    double bsMax = 1., bsMin = 0.;
    int itrc = 0;

    while (!finished) {
      itrc++;

      r = maxRadius * (bsMax + bsMin) / 2;

      auto path = cutTracePath(mesh, geometry, wholePath, cartesianCoords, r);
      auto p = mmp.traceBack(path.back()).back();

      if (isIdenticalSurfacePoint(p, x)) {
        bsMin = (bsMax + bsMin) / 2;
      } else {
        bsMax = (bsMax + bsMin) / 2;
      }

      if (itrc > 10) {
        break;
      }
    }

    return r;
  }

  std::vector<std::vector<Vector3>> medial_axis_euclidian(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry,
    const std::vector<SurfacePoint> &nodes,
    const std::vector<std::array<int, 2>> &segments,
    const std::vector<Vector3> &cartesianCoords,
    const std::vector<Vector3> &nodeTangents,
    const std::vector<Vector3> &nodeNormals,
    const std::vector<Vector3> &nodeBitangents,
    const double maxRadius
  ) {
    // create a closest point query
    Eigen::MatrixXd _V(3, nodes.size());
    for (int i = 0; i < cartesianCoords.size(); i++) {
      _V.col(i) << cartesianCoords[i].x, cartesianCoords[i].y, cartesianCoords[i].z;
    }
    knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> kdtree(_V);
    kdtree.build();

    std::vector<std::vector<Vector3>> nodeMedialAxis(nodes.size());

    #pragma omp parallel for
    for (int i = 0; i < nodes.size(); i++) {
      auto t = nodeTangents[i];
      auto n = nodeNormals[i];
      auto b = nodeBitangents[i];
      auto x = cartesianCoords[i];

      double r_min_minus = std::numeric_limits<double>::infinity();
      double r_min_plus = std::numeric_limits<double>::infinity();

      r_min_plus = std::min(maximumBallRadius(x, b, i, maxRadius, cartesianCoords, kdtree), maxRadius);
      r_min_minus = std::min(maximumBallRadius(x, -b, i, maxRadius, cartesianCoords, kdtree), maxRadius);

      // std::cout << "result: " << r_min_minus << ", " << r_min_plus << std::endl;

      nodeMedialAxis[i].emplace_back(cartesianCoords[i] - r_min_minus * b);
      nodeMedialAxis[i].emplace_back(cartesianCoords[i] + r_min_plus * b);
    }

    return nodeMedialAxis;
  }

  std::vector<std::vector<Vector3>> medial_axis_geodesic(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry,
    const std::vector<SurfacePoint> &nodes,
    const std::vector<std::array<int, 2>> &segments,
    const std::vector<Vector3> &cartesianCoords,
    const std::vector<Vector3> &nodeTangents,
    const std::vector<Vector3> &nodeNormals,
    const std::vector<Vector3> &nodeBitangents,
    const double maxRadius
  ) {
    // create a geodesic closest point query
    GeodesicAlgorithmExact mmp(mesh, geometry);
    mmp.propagate(nodes);

    std::vector<std::vector<Vector3>> nodeMedialAxis(nodes.size());

    // for debug
    // Eigen::MatrixXd _V(3, nodes.size());
    // for (int i = 0; i < cartesianCoords.size(); i++) {
    //   _V.col(i) << cartesianCoords[i].x, cartesianCoords[i].y, cartesianCoords[i].z;
    // }
    // knncpp::KDTreeMinkowskiX<double, knncpp::EuclideanDistance<double>> kdtree(_V);
    // kdtree.build();

    for (int i = 0; i < nodes.size(); i++) {
      auto t = nodeTangents[i];
      auto n = nodeNormals[i];
      auto b = nodeBitangents[i];
      auto node = nodes[i];

      double r_min_plus = maximumGeodesicBallRadius(mesh, geometry, node, b, i, maxRadius, mmp);
      double r_min_minus = maximumGeodesicBallRadius(mesh, geometry, node, -b, i, maxRadius, mmp);

      {
        // for debug

        // double r_euclidian_plus = maximumBallRadius(cartesianCoords[i], b, i, maxRadius, cartesianCoords, kdtree);
        // double r_euclidian_minus = maximumBallRadius(cartesianCoords[i], -b, i, maxRadius, cartesianCoords, kdtree);

        // std::cout << "euclidian: " << r_euclidian_minus << ", " << r_euclidian_plus << std::endl;
        // std::cout << "geodesic: " << r_min_minus << ", " << r_min_plus << std::endl;
      }

      r_min_plus = std::min(r_min_plus, maxRadius);
      r_min_minus = std::min(r_min_minus, maxRadius);

      // std::cout << "result: " << r_min_minus << ", " << r_min_plus << std::endl;

      nodeMedialAxis[i].emplace_back(cartesianCoords[i] - r_min_minus * b);
      nodeMedialAxis[i].emplace_back(cartesianCoords[i] + r_min_plus * b);
    }

    return nodeMedialAxis;
  }
  
  std::vector<std::vector<Vector3>> medial_axis(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry,
    const std::vector<SurfacePoint> &nodes,
    const std::vector<std::array<int, 2>> &segments,
    const std::vector<Vector3> &cartesianCoords,
    const std::vector<Vector3> &nodeTangents,
    const std::vector<Vector3> &nodeNormals,
    const std::vector<Vector3> &nodeBitangents,
    const double maxRadius,
    const bool isGeodesic
  ) {
    if (isGeodesic) {
      return medial_axis_geodesic(mesh, geometry, nodes, segments, cartesianCoords, nodeTangents, nodeNormals, nodeBitangents, maxRadius);
    } else {
      return medial_axis_euclidian(mesh, geometry, nodes, segments, cartesianCoords, nodeTangents, nodeNormals, nodeBitangents, maxRadius);
    }
  }
}