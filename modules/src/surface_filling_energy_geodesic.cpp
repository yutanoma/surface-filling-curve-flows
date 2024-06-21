#include "modules/surface_filling_energy_geodesic.h"

#include "modules/tet_topology.h"
#include "modules/surface_point_tangent_basis.h"
#include "modules/surface_point_to_cartesian.h"
#include "modules/surface_point_vector_field.h"
#include "modules/surface_point_value.h"
#include "modules/medial_axis.h"

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>

#include <igl/octree.h>
#include <igl/knn.h>
#include <igl/PI.h>

namespace modules {
const double branchRatio = std::sqrt(std::sqrt(2));

std::tuple<
  std::vector<Vector3> /* descent direction */,
  std::vector<Vector3> /* gradient direction */,
  double /* energy */,
  std::vector<std::vector<Vector3>> /* medial axis */
> surface_filling_energy_geodesic(
  ManifoldSurfaceMesh &mesh,
  VertexPositionGeometry &geometry,
  const std::vector<SurfacePoint> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<std::vector<SurfacePoint>> &segmentSurfacePoints,
  const std::vector<double> &segmentLengths,
  const std::vector<bool> &isFixedNode,
  const std::vector<Vector3> &cartesianCoords,
  const double radius,
  const double maxRadius,
  const SurfaceFillingEnergy::Options &options
) {
  auto start = std::chrono::high_resolution_clock::now();

  double p = options.p;
  double q = options.q;
  double branchRadius = radius * branchRatio;
  double alpha = 4 / (pow(branchRadius, 2));

  std::cout << "alpha: " << alpha << ", p: " << p << ", q: " << q  << std::endl;

  std::vector<Vector3> descent(nodes.size(), Vector3 {.0, .0, .0}), gradient(nodes.size(), Vector3 {.0, .0, .0});

  assert(nodes.size() == isFixedNode.size());
  assert(nodes.size() == cartesianCoords.size());

  std::map<int, int> activeNode2NodeIdx, node2ActiveNodeIdx;
  for (int i = 0; i < nodes.size(); i++) {
    if (!isFixedNode[i]) {
      activeNode2NodeIdx[activeNode2NodeIdx.size()] = i;
      node2ActiveNodeIdx[i] = node2ActiveNodeIdx.size();
    }
  }
  int numActiveNodes = activeNode2NodeIdx.size();

  double totalCurveLength = .0;
  for (int i = 0; i < segments.size(); i++) {
    totalCurveLength += segmentLengths[i];
  }

  auto func = TinyAD::scalar_function<3>(TinyAD::range(numActiveNodes));

  // ===== add dirichlet energy =====
  std::vector<std::array<Eigen::Matrix3d, 2>> rotationMatrix(segments.size());
  std::vector<std::array<Vector3, 2>> segmentTangent(segments.size());
  std::vector<std::array<Vector3, 2>> segmentNormal(segments.size());
  std::vector<std::array<Vector3, 2>> segmentBitangent(segments.size());

  for (int i = 0; i < segments.size(); i++) {
    std::vector<SurfacePoint> pointsOnSegment = {nodes[segments[i][0]]};
    for (auto sp : segmentSurfacePoints[i]) {
      pointsOnSegment.push_back(sp);
    }
    pointsOnSegment.emplace_back(nodes[segments[i][1]]);

    auto _edgeCartesians = modules::surface_point_to_cartesian(mesh, geometry, pointsOnSegment);
    std::vector<Vector3> edgeCartesians = {_edgeCartesians[0]};
    for (int j = 1; j < _edgeCartesians.size(); j++) {
      auto prev = edgeCartesians[edgeCartesians.size() - 1];
      if ((_edgeCartesians[j] - prev).norm() > 1e-6) {
        edgeCartesians.emplace_back(_edgeCartesians[j]);
      }
    }
    if (edgeCartesians.size() == 1) {
      edgeCartesians.emplace_back(_edgeCartesians[_edgeCartesians.size() - 1]);
    }

    assert(edgeCartesians.size() > 1);

    std::array<Vector3, 2> tangents;
    tangents[0] = normalize(edgeCartesians[1] - edgeCartesians[0]);
    tangents[1] = normalize(edgeCartesians[edgeCartesians.size() - 1] - edgeCartesians[edgeCartesians.size() - 2]);

    for (int j = 0; j < 2; j++) {
      int v = segments[i][j];
      SurfacePoint sp = nodes[v];

      auto [x, y] = surface_point_tangent_basis(geometry, sp);

      Vector3 t = normalize(dot(x, tangents[j]) * x + dot(y, tangents[j]) * y);
      Vector3 n = cross(x, y);
      Vector3 b = cross(n, t);

      segmentTangent[i][j] = t;
      segmentNormal[i][j] = n;
      segmentBitangent[i][j] = b;

      Eigen::Matrix3d R;
      R << t.x, b.x, n.x,
            t.y, b.y, n.y,
            t.z, b.z, n.z;

      rotationMatrix[i][j] = R;
    }
  }

  std::vector<std::pair<int, int>> segmentsWith2ActiveNodes = {};
  std::vector<int> activeTwoSegment2SegmentIdx = {};
  // first node is active while the second node is fixed
  std::vector<std::pair<int, int>> segmentsWith1ActiveNode = {};
  std::vector<int> activeOneSegment2SegmentIdx = {};
  std::vector<double> activeOneSegmentSigns = {};

  for (int i = 0; i < segments.size(); i++) {
    int v0 = segments[i][0], v1 = segments[i][1];

    if (isFixedNode[v0] && isFixedNode[v1]) {
      continue;
    } else if (!isFixedNode[v0] && !isFixedNode[v1]) {
      int newId = segmentsWith2ActiveNodes.size();
      segmentsWith2ActiveNodes.emplace_back(std::pair<int, int> {v0, v1});
      activeTwoSegment2SegmentIdx.emplace_back(i);
    } else {
      int activeNode = isFixedNode[v0] ? v1 : v0;
      int fixedNode = isFixedNode[v0] ? v0 : v1;
      segmentsWith1ActiveNode.emplace_back(std::pair<int, int> {activeNode, fixedNode});
      activeOneSegment2SegmentIdx.emplace_back(i);

      if (isFixedNode[v0]) {
        activeOneSegmentSigns.emplace_back(-1);
      } else {
        activeOneSegmentSigns.emplace_back(1);
      }
    }
  }

  assert(segmentsWith2ActiveNodes.size() == activeTwoSegment2SegmentIdx.size());
  assert(segmentsWith1ActiveNode.size() == activeOneSegment2SegmentIdx.size());
  assert(segmentsWith1ActiveNode.size() == activeOneSegmentSigns.size());

  std::cout << "field alignedness: " << options.w_fieldAlignedness << ", curvature alignedness: " << options.w_curvatureAlignedness << ", bilaplacian weight: " << options.w_bilaplacian << std::endl;

  // ===== add field-aligned energy =====
  std::vector<Vector3> vectorFieldOnNode(nodes.size());
  if (options.w_fieldAlignedness > 0) {
    for (int i = 0; i < nodes.size(); i++) {
      auto [vf] = modules::surface_point_vector_field(geometry, options.vectorField, nodes[i]);
      vectorFieldOnNode[i] = vf;
    }
  }

  std::vector<Vector3> principalCurvatureOnNode(nodes.size());
  if (options.w_curvatureAlignedness > 0) {
    for (int i = 0; i < nodes.size(); i++) {
      auto [vf] = modules::surface_point_vector_field(geometry, geometry.vertexPrincipalCurvatureDirections, nodes[i]);
      principalCurvatureOnNode[i] = vf;
    }
  }

  func.add_elements<2>(TinyAD::range(segmentsWith2ActiveNodes.size()), [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int e_id = element.handle;

    int _v0 = segmentsWith2ActiveNodes[e_id].first;
    int _v1 = segmentsWith2ActiveNodes[e_id].second;

    int segmentId = activeTwoSegment2SegmentIdx[e_id];

    Eigen::Vector3<T> x0 = element.variables(node2ActiveNodeIdx[_v0]);
    Eigen::Vector3<T> x1 = element.variables(node2ActiveNodeIdx[_v1]);

    Eigen::Vector3d v0(
      cartesianCoords[_v0].x,
      cartesianCoords[_v0].y,
      cartesianCoords[_v0].z
    );
    Eigen::Vector3d v1(
      cartesianCoords[_v1].x,
      cartesianCoords[_v1].y,
      cartesianCoords[_v1].z
    );

    double edgeLen = segmentLengths[segmentId];

    Eigen::Matrix3d R0 = rotationMatrix[segmentId][0];
    Eigen::Matrix3d R1 = rotationMatrix[segmentId][1];

    Eigen::Vector3<T> p0 = R0.transpose() * (x0 - v0);
    Eigen::Vector3<T> p1 = R1.transpose() * (x1 - v1) + Eigen::Vector3d(edgeLen, 0, 0);

    // std::cout << "id: " << _v0 << ", " << _v1 << std::endl;
    // std::cout << "v0(" << segments[segmentId][0] << "): " << v0 << std::endl << "R0: " << R0.transpose() << std::endl;
    // std::cout << "v1(" << segments[segmentId][1] << "): " << v1 << std::endl << "R1: " << R1.transpose() << std::endl;

    T dx = abs(p0(0) - p1(0));
    T dy = abs(p0(1) - p1(1));
    T dz = abs(p0(2) - p1(2));

    Eigen::Vector3d vf0(
      vectorFieldOnNode[_v0].x,
      vectorFieldOnNode[_v0].y,
      vectorFieldOnNode[_v0].z
    );
    Eigen::Vector3d vf1(
      vectorFieldOnNode[_v1].x,
      vectorFieldOnNode[_v1].y,
      vectorFieldOnNode[_v1].z
    );

    Eigen::Vector3d vec0 = R0.transpose() * vf0;
    Eigen::Vector3d vec1 = R1.transpose() * vf1;

    T crs0 = options.w_fieldAlignedness * (p1 - p0).cross(vec0).squaredNorm() / 2;
    T crs1 = options.w_fieldAlignedness * (p0 - p1).cross(vec1).squaredNorm() / 2;

    Eigen::Vector3d pc0(
      principalCurvatureOnNode[_v0].x,
      principalCurvatureOnNode[_v0].y,
      principalCurvatureOnNode[_v0].z
    );
    Eigen::Vector3d pc1(
      principalCurvatureOnNode[_v1].x,
      principalCurvatureOnNode[_v1].y,
      principalCurvatureOnNode[_v1].z
    );

    Eigen::Vector3d pc0_ = R0.transpose() * pc0;
    Eigen::Vector3d pc1_ = R1.transpose() * pc1;

    T dot0 = options.w_curvatureAlignedness * pow((p1 - p0).dot(pc0_), 2) / 2;
    T dot1 = options.w_curvatureAlignedness * pow((p0 - p1).dot(pc1_), 2) / 2;

    // std::cout << "[" << segmentId << "]: dx: " << dx << ", dy: " << dy << ", dz: " << dz << std::endl;

    return (
      pow(pow(dx, 2) + pow(dy, 2) + pow(dz, 2), p / 2) +
      crs0 + crs1 + dot0 + dot1
    ) / (edgeLen * totalCurveLength);
  });

  func.add_elements<1>(TinyAD::range(segmentsWith1ActiveNode.size()), [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int e_id = element.handle;
    int segmentId = activeOneSegment2SegmentIdx[e_id];

    int _v0 = segmentsWith1ActiveNode[e_id].first;
    int _v1 = segmentsWith1ActiveNode[e_id].second;

    double sgn = activeOneSegmentSigns[e_id];

    Eigen::Vector3<T> x0 = element.variables(node2ActiveNodeIdx[_v0]);

    Eigen::Vector3d v0(
      cartesianCoords[_v0].x,
      cartesianCoords[_v0].y,
      cartesianCoords[_v0].z
    );

    double edgeLen = segmentLengths[activeOneSegment2SegmentIdx[e_id]];

    Eigen::Matrix3d R0 = sgn > 0 ? rotationMatrix[segmentId][0] : rotationMatrix[segmentId][1];
    Eigen::Matrix3d R1 = sgn > 0 ? rotationMatrix[segmentId][1] : rotationMatrix[segmentId][0];
    Eigen::Vector3<T> p0 = sgn * R0.transpose() * (x0 - v0);
    Eigen::Vector3d p1 = Eigen::Vector3d(edgeLen, 0, 0);

    T dx = abs(p0(0) - p1(0));
    T dy = abs(p0(1) - p1(1));
    T dz = abs(p0(2) - p1(2));

    Eigen::Vector3d vf0(
      vectorFieldOnNode[_v0].x,
      vectorFieldOnNode[_v0].y,
      vectorFieldOnNode[_v0].z
    );
    Eigen::Vector3d vf1(
      vectorFieldOnNode[_v1].x,
      vectorFieldOnNode[_v1].y,
      vectorFieldOnNode[_v1].z
    );

    Eigen::Vector3d vec0 = R0.transpose() * vf0;
    Eigen::Vector3d vec1 = R1.transpose() * vf1;

    T crs0 = options.w_fieldAlignedness * (p1 - p0).cross(vec0).squaredNorm() / 2;
    T crs1 = options.w_fieldAlignedness * (p0 - p1).cross(vec1).squaredNorm() / 2;

    Eigen::Vector3d pc0(
      principalCurvatureOnNode[_v0].x,
      principalCurvatureOnNode[_v0].y,
      principalCurvatureOnNode[_v0].z
    );
    Eigen::Vector3d pc1(
      principalCurvatureOnNode[_v1].x,
      principalCurvatureOnNode[_v1].y,
      principalCurvatureOnNode[_v1].z
    );
    
    Eigen::Vector3d pc0_ = R0.transpose() * pc0;
    Eigen::Vector3d pc1_ = R1.transpose() * pc1;

    T dot0 = options.w_curvatureAlignedness * pow((p1 - p0).dot(pc0_), 2) / 2;
    T dot1 = options.w_curvatureAlignedness * pow((p0 - p1).dot(pc1_), 2) / 2;

    return (
      pow(pow(dx, 2) + pow(dy, 2) + pow(dz, 2), p / 2) +
      crs0 + crs1 + dot0 + dot1
    ) / (edgeLen * totalCurveLength);
  });

  // ===== add bilaplacian energy =====
  std::vector<std::array<int, 3>> vertexTriplets = {};
  std::vector<std::array<Eigen::Matrix3d, 3>> rotationMatrices = {};
  std::vector<std::array<Eigen::Vector3d, 3>> mappedPoints = {};
  std::vector<std::array<double, 2>> mappedLengths = {};

  std::vector<std::array<int, 3>> vertexTripletsWith2ActiveNodes = {};
  std::vector<std::array<Eigen::Matrix3d, 3>> rotationMatricesWith2ActiveNodes = {};
  std::vector<std::array<Eigen::Vector3d, 3>> mappedPointsWith2ActiveNodes = {};
  std::vector<std::array<double, 2>> mappedLengthsWith2ActiveNodes = {};

  std::vector<std::array<int, 3>> vertexTripletsWith1ActiveNode = {};
  std::vector<std::array<Eigen::Matrix3d, 3>> rotationMatricesWith1ActiveNode = {};
  std::vector<std::array<Eigen::Vector3d, 3>> mappedPointsWith1ActiveNode = {};
  std::vector<std::array<double, 2>> mappedLengthsWith1ActiveNode = {};

  if (options.w_bilaplacian > .0) {
    std::vector<std::vector<int>> node2Segment(nodes.size());
    for (int i = 0; i < segments.size(); i++) {
      for (int j = 0; j < 2; j++) {
        node2Segment[segments[i][j]].emplace_back(i);
      }
    }

    for (int i = 0; i < nodes.size(); i++) {
      if (isFixedNode[i]) {
        continue;
      }

      auto incidentSegments = node2Segment[i];

      if (incidentSegments.size() != 2) {
        // do not manupulate non-manifold points or loose-end points
        continue;
      }

      int s0 = incidentSegments[0];
      int s1 = incidentSegments[1];

      int v1 = i;

      int v0 = segments[s0][0] == v1 ? segments[s0][1] : segments[s0][0];
      int v2 = segments[s1][0] == v1 ? segments[s1][1] : segments[s1][0];

      double sgn0 = segments[s0][0] == v1 ? 1 : -1;
      double sgn1 = segments[s1][0] == v1 ? 1 : -1;

      Vector3 t0 = sgn0 > 0 ? segmentTangent[s0][1] : -segmentTangent[s0][0];
      Vector3 b0 = sgn0 > 0 ? segmentBitangent[s0][1] : -segmentBitangent[s0][0];
      Vector3 n0 = cross(t0, b0);

      Vector3 t1 = sgn0 > 0 ? segmentTangent[s0][0] : -segmentTangent[s0][1];
      Vector3 b1 = sgn0 > 0 ? segmentBitangent[s0][0] : -segmentBitangent[s0][1];
      Vector3 n1 = cross(t1, b1);

      Vector3 t2 = sgn1 > 0 ? -segmentTangent[s1][1] : segmentTangent[s1][0];
      Vector3 b2 = sgn1 > 0 ? -segmentBitangent[s1][1] : segmentBitangent[s1][0];
      Vector3 n2 = cross(t2, b2);

      Vector3 t12 = sgn1 > 0 ? segmentTangent[s1][0] : -segmentTangent[s1][1];
      double x = dot(t1, t12);
      double y = dot(t12, b1);
      // angle from 1->0 to 1->2
      double angle = std::atan2(y, x);

      Eigen::Matrix3d R0;
      R0 << t0.x, t0.y, t0.z,
            b0.x, b0.y, b0.z,
            n0.x, n0.y, n0.z;

      Eigen::Matrix3d R1;
      R1 << t1.x, t1.y, t1.z,
            b1.x, b1.y, b1.z,
            n1.x, n1.y, n1.z;

      Eigen::Matrix3d R2_rot;
      double rot = angle - igl::PI;
      R2_rot << std::cos(rot), -std::sin(rot), 0,
                std::sin(rot), std::cos(rot), 0,
                0, 0, 1;

      Eigen::Matrix3d R2;
      R2 << t2.x, t2.y, t2.z,
            b2.x, b2.y, b2.z,
            n2.x, n2.y, n2.z;
      
      R2 = R2_rot * R2;

      double l0 = segmentLengths[s0];
      double l1 = segmentLengths[s1];

      Eigen::Vector3d p0 {l0, 0, 0};
      Eigen::Vector3d p1 {0, 0, 0};
      Eigen::Vector3d p2 {l1 * std::cos(angle), l1 * std::sin(angle), 0};

      // std::cout << "x: " << x << ", y: " << y << ", angle: "<< angle << std::endl;
      // std::cout << R2 * Eigen::Vector3d(t2.x, t2.y, t2.z) << std::endl;

      assert((R2 * Eigen::Vector3d(t2.x, t2.y, t2.z) - Eigen::Vector3d(-x, -y, 0)).norm() < 1e-6);

      bool is1Active = isFixedNode[v0] && isFixedNode[v2];
      bool is2Active = isFixedNode[v0] || isFixedNode[v2];

      // make sure to add it in the order of v0 -> v1 -> v2
      if (is1Active) {
        vertexTripletsWith1ActiveNode.emplace_back(std::array<int, 3> {v0, v1, v2});
        rotationMatricesWith1ActiveNode.emplace_back(std::array<Eigen::Matrix3d, 3> {
          R0, R1, R2
        });
        mappedPointsWith1ActiveNode.emplace_back(std::array<Eigen::Vector3d, 3> {
          p0, p1, p2
        });
        mappedLengthsWith1ActiveNode.emplace_back(std::array<double, 2> {l0, l1});
      } else if (is2Active) {
        vertexTripletsWith2ActiveNodes.emplace_back(std::array<int, 3> {v0, v1, v2});
        rotationMatricesWith2ActiveNodes.emplace_back(std::array<Eigen::Matrix3d, 3> {
          R0, R1, R2
        });
        mappedPointsWith2ActiveNodes.emplace_back(std::array<Eigen::Vector3d, 3> {
          p0, p1, p2
        });
        mappedLengthsWith2ActiveNodes.emplace_back(std::array<double, 2> {l0, l1});
      } else {
        vertexTriplets.emplace_back(std::array<int, 3> {v0, v1, v2});
        rotationMatrices.emplace_back(std::array<Eigen::Matrix3d, 3> {
          R0, R1, R2
        });
        mappedPoints.emplace_back(std::array<Eigen::Vector3d, 3> {
          p0, p1, p2
        });
        mappedLengths.emplace_back(std::array<double, 2> {l0, l1});
      }
    }
  }

  func.add_elements<3>(TinyAD::range(vertexTriplets.size()), [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int e_id = element.handle;

    int _v0 = vertexTriplets[e_id][0];
    int _v1 = vertexTriplets[e_id][1];
    int _v2 = vertexTriplets[e_id][2];

    Eigen::Vector3<T> x0 = element.variables(node2ActiveNodeIdx[_v0]);
    Eigen::Vector3<T> x1 = element.variables(node2ActiveNodeIdx[_v1]);
    Eigen::Vector3<T> x2 = element.variables(node2ActiveNodeIdx[_v2]);

    Eigen::Vector3d v0(
      cartesianCoords[_v0].x,
      cartesianCoords[_v0].y,
      cartesianCoords[_v0].z
    );
    Eigen::Vector3d v1(
      cartesianCoords[_v1].x,
      cartesianCoords[_v1].y,
      cartesianCoords[_v1].z
    );
    Eigen::Vector3d v2(
      cartesianCoords[_v2].x,
      cartesianCoords[_v2].y,
      cartesianCoords[_v2].z
    );

    // flatten three vertices on a single plane
    Eigen::Matrix3d R0 = rotationMatrices[e_id][0];
    Eigen::Matrix3d R1 = rotationMatrices[e_id][1];
    Eigen::Matrix3d R2 = rotationMatrices[e_id][2];

    Eigen::Vector3<T> p0 = R0 * (x0 - v0) + mappedPoints[e_id][0];
    Eigen::Vector3<T> p1 = R1 * (x1 - v1) + mappedPoints[e_id][1];
    Eigen::Vector3<T> p2 = R2 * (x2 - v2) + mappedPoints[e_id][2];

    // std::cout << "id: " << _v0 << ", " << _v1 << ", " << _v2 << std::endl;
    // std::cout << "bilaplacian v0: " << v0 << std::endl << "R0: " << R0 << std::endl;
    // std::cout << "bilaplacian v1: " << v1 << std::endl << "R1: " << R1 << std::endl;

    double l0 = mappedLengths[e_id][0];
    double l1 = mappedLengths[e_id][1];

    // std::cout << "l0: " << l0 << ", l1: " << l1 << std::endl;

    Eigen::Vector3<T> d0 = (p0 - p1) / l0;
    Eigen::Vector3<T> d1 = (p1 - p2) / l1;

    // for debug
    // T _d = (l0 * pow(d0.norm(), 2) + l1 * pow(d1.norm(), 2)) / 2;
    // // T _d = (pow(d0, 2) / l0 + pow(d1, 2) / l1) / 2;
    // return (options.w_bilaplacian * _d) / totalCurveLength;

    T d = pow((d0 - d1).norm(), 2);
    double l = (l0 + l1) / 2;

    return (options.w_bilaplacian * l * d) / totalCurveLength;
  });

  auto dirichletEnd = std::chrono::high_resolution_clock::now();

  // ===== add medial axis energy =====

  auto tetTopologyEnd = std::chrono::high_resolution_clock::now();

  std::vector<Vector3> nodeTangents(nodes.size());
  std::vector<Vector3> nodeNormals(nodes.size());
  std::vector<Vector3> nodeBitangents(nodes.size());
  std::vector<double> nodeWeight(nodes.size(), .0);
  for (int i = 0; i < segments.size(); i++) {
    for (int j = 0; j < 2; j++) {
      auto t = segmentTangent[i][j];
      auto n = segmentNormal[i][j];
      auto b = segmentBitangent[i][j];

      int v = segments[i][j];

      nodeTangents[v] += normalize(t);
      nodeNormals[v] += normalize(n);
      nodeBitangents[v] += normalize(b);

      // for extremely challenging cases, we need to add a small perturbation
      if (norm(nodeTangents[v]) < .00001) {
        nodeTangents[v] -= .01 * normalize(n);
      }
      if (norm(nodeNormals[v]) < .00001) {
        nodeNormals[v] -= .01 * normalize(t);
      }
      if (norm(nodeBitangents[v]) < .00001) {
        nodeBitangents[v] -= .01 * normalize(t);
      }

      nodeWeight[v] += segmentLengths[i] / 2;
    }
  }

  for (int i = 0; i < nodes.size(); i++) {
    nodeTangents[i] = normalize(nodeTangents[i]);
    nodeNormals[i] = normalize(nodeNormals[i]);
    nodeBitangents[i] = normalize(nodeBitangents[i]);
  }

  auto tangentEvalEnd = std::chrono::high_resolution_clock::now();

  auto nodeMedialAxis = modules::medial_axis(mesh, geometry, nodes, segments, cartesianCoords, nodeTangents, nodeNormals, nodeBitangents, maxRadius, options.useGeodesicMedialAxis);

  auto medialAxisEnd = std::chrono::high_resolution_clock::now();
  
  std::vector<double> alphas(nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    if (options.useAnisotropicAlphaOnNodes) {
      alphas[i] = options.alphaRatioOnNodes[i] * alpha;
    } else if (options.useAnisotropicAlphaOnMesh) {
      alphas[i] = modules::surface_point_value(geometry, options.alphaRatioOnMesh, nodes[i]) * alpha;
      // std::cout << alphas[i] << std::endl;
    } else {
      alphas[i] = alpha;
    }
  }

  func.add_elements<1>(TinyAD::range(numActiveNodes), [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int e_id = element.handle;
    int nodeId = activeNode2NodeIdx[e_id];

    auto _c_0 = nodeMedialAxis[nodeId][0];
    auto _c_1 = nodeMedialAxis[nodeId][1];

    Eigen::Vector3d c_0(_c_0.x, _c_0.y, _c_0.z);
    Eigen::Vector3d c_1(_c_1.x, _c_1.y, _c_1.z);

    T l_0 = (element.variables(e_id) - c_0).norm();
    T l_1 = (element.variables(e_id) - c_1).norm();

    // std::cout << "[" << nodeId << "]: l0: " << l_0 << ", l1: " << l_1 << ", c_0: " << c_0 << ", c_1: " << c_1 << std::endl;

    return alphas[nodeId] * nodeWeight[nodeId] * (
      pow(pow(l_0, 2) + pow(l_1, 2), q / 2)
    ) / totalCurveLength;
  });

  auto repulsiveEnd = std::chrono::high_resolution_clock::now();

  auto x = func.x_from_data([&](int v_idx) {
    auto v = cartesianCoords[activeNode2NodeIdx[v_idx]];

    return Eigen::Vector3d(v.x, v.y, v.z);
  });

  TinyAD::LinearSolver solver;
  auto [f, g, H_proj] = func.eval_with_hessian_proj(x);

  auto hessianEnd = std::chrono::high_resolution_clock::now();

  auto d = TinyAD::newton_direction(g, H_proj, solver);

  auto newtonEnd = std::chrono::high_resolution_clock::now();

  func.x_to_data(d, [&] (int v_idx, const Eigen::Vector3d& p) {
    auto idx = activeNode2NodeIdx[v_idx];

    descent[idx] = Vector3 {p(0), p(1), p(2)};
  });

  func.x_to_data(g, [&] (int v_idx, const Eigen::Vector3d& p) {
    auto idx = activeNode2NodeIdx[v_idx];

    gradient[idx] = Vector3 {p(0), p(1), p(2)};
  });

  auto end = std::chrono::high_resolution_clock::now();

  if (true) {
    // std::ofstream file("./hessian.txt");

    // for (int i = 0; i < H_proj.rows(); i++) {
    //   for (int j = 0; j < H_proj.cols(); j++) {
    //     file << H_proj.coeff(i, j);

    //     if (j < H_proj.cols() - 1) {
    //       file << " ";
    //     }
    //   }
    //   file << std::endl;
    // }

    // file.close();
  }

  std::cout << "optimization: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
  std::cout << "  rigid transformation comp: " << std::chrono::duration_cast<std::chrono::milliseconds>(dirichletEnd - start).count() << "ms" << std::endl;
  // std::cout << "  tet topology: " << std::chrono::duration_cast<std::chrono::milliseconds>(tetTopologyEnd - dirichletEnd).count() << "ms" << std::endl;
  std::cout << "  tangent space comp: " << std::chrono::duration_cast<std::chrono::milliseconds>(tangentEvalEnd - tetTopologyEnd).count() << "ms" << std::endl;
  std::cout << "  medial axis comp: " << std::chrono::duration_cast<std::chrono::milliseconds>(medialAxisEnd - tangentEvalEnd).count() << "ms" << std::endl;
  std::cout << "  energy comp: " << std::chrono::duration_cast<std::chrono::milliseconds>(repulsiveEnd - medialAxisEnd).count() << "ms" << std::endl;
  std::cout << "  hessian evaluation: " << std::chrono::duration_cast<std::chrono::milliseconds>(hessianEnd - repulsiveEnd).count() << "ms" << std::endl;
  std::cout << "  linear solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(newtonEnd - hessianEnd).count() << "ms" << std::endl;
  std::cout << "    d norm: " << d.norm() << std::endl;
  std::cout << "    g norm: " << g.norm() << std::endl;
  std::cout << "    energy: " << f << std::endl;

  return {
    descent,
    gradient,
    f,
    nodeMedialAxis
  };
}
}
