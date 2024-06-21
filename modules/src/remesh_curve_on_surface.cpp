#include <modules/remesh_curve_on_surface.h>

#include <modules/surface_point_to_cartesian.h>
#include <modules/point_point_geodesic.h>
#include <modules/exact_geodesics.h>

namespace modules {
std::tuple<
  std::vector<SurfacePoint>,
  std::vector<std::array<int, 2>>,
  std::vector<std::vector<SurfacePoint>>,
  std::vector<double>,
  std::vector<bool>
> removeShortEdges(
  ManifoldSurfaceMesh &mesh,
  VertexPositionGeometry &geometry,
  const std::vector<SurfacePoint> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<std::vector<SurfacePoint>> &segmentSurfacePoints,
  const std::vector<double> &segmentLengths,
  const std::vector<bool> &isFixedNode,
  const double &h
) {
  std::map<int, bool> deletingNodes;

  // 1. delete nodes that are too close to others
  std::vector<std::vector<int>> node2Segments(nodes.size(), std::vector<int>{});
  for (int i = 0; i < segments.size(); i++) {
    node2Segments[segments[i][0]].emplace_back(i);
    node2Segments[segments[i][1]].emplace_back(i);
  }

  for (int i = 0; i < segments.size(); i++) {
    double edgeLen = segmentLengths[i];

    if (isFixedNode[segments[i][0]] && isFixedNode[segments[i][1]]) {
      continue;
    }

    if (edgeLen < h) {
      double shortestLen = INFINITY;
      int shortestV = -1;

      for (int j = 0; j < 2; j++) {
        int v = segments[i][j];

        if (isFixedNode[v]) {
          continue;
        }

        if (node2Segments[v].size() != 2) {
          continue;
        }

        int otherS = node2Segments[v][0] == i ? node2Segments[v][1] : node2Segments[v][0];

        if (segmentLengths[otherS] < shortestLen && !deletingNodes[v]) {
          shortestLen = segmentLengths[otherS];
          shortestV = v;
        }
      }

      if (shortestV != -1) {
        deletingNodes[shortestV] = true;
      }
    }
  }

  std::map<int, int> node2NewNode;
  std::vector<SurfacePoint> newNodes = {};
  std::vector<bool> newNodeIsFixed = {};

  for (int i = 0; i < nodes.size(); i++) {
    if (deletingNodes[i]) {
      continue;
    }

    node2NewNode[i] = newNodes.size();
    newNodes.emplace_back(nodes[i]);
    newNodeIsFixed.emplace_back(isFixedNode[i]);
  }

  // 2. update segments
  std::vector<std::array<int, 2>> newSegments = {};
  std::vector<std::vector<SurfacePoint>> newSegmentSurfacePoints = {};
  std::vector<double> newSegmentLengths = {};

  modules::GeodesicAlgorithmExact mmp(mesh, geometry);

  for (int i = 0; i < segments.size(); i++) {
    int v0 = segments[i][0], v1 = segments[i][1];

    if (!deletingNodes[v0] && !deletingNodes[v1]) {
      newSegments.emplace_back(std::array<int, 2>{node2NewNode[v0], node2NewNode[v1]});
      newSegmentSurfacePoints.emplace_back(segmentSurfacePoints[i]);
      newSegmentLengths.emplace_back(segmentLengths[i]);
    } else if (!deletingNodes[v0] && deletingNodes[v1]) {
      // traverse the curve until we find a node that is not deleted
      int v = v1;
      int currentSegment = i;
      while (deletingNodes[v]) {
        assert(node2Segments[v].size() == 2);

        currentSegment = node2Segments[v][0] == currentSegment ? node2Segments[v][1] : node2Segments[v][0];

        v = segments[currentSegment] [0] == v ? segments[currentSegment][1] : segments[currentSegment][0];
      }

      newSegments.emplace_back(std::array<int, 2>{node2NewNode[v0], node2NewNode[v]});

      auto [path] = point_point_geodesic(mesh, geometry, mmp, nodes[v0], nodes[v]);

      double length = .0;

      auto cartesianCoord = surface_point_to_cartesian(mesh, geometry, path);
      std::vector<SurfacePoint> edgeSurfacePoints = {};

      for (int j = 0; j < path.size(); j++) {
        if (j != 0) {
          length += (cartesianCoord[j] - cartesianCoord[j - 1]).norm();
        }

        if (j != 0 && j != path.size() - 1) {
          edgeSurfacePoints.emplace_back(path[j]);
        }
      }

      newSegmentSurfacePoints.emplace_back(edgeSurfacePoints);
      newSegmentLengths.emplace_back(length);
    } else {
      // do nothing
    }
  }

  assert(newSegmentSurfacePoints.size() == newSegments.size());
  assert(newSegmentLengths.size() == newSegments.size());
  
  return {
    newNodes,
    newSegments,
    newSegmentSurfacePoints,
    newSegmentLengths,
    newNodeIsFixed
  };
}

std::tuple<
  std::vector<SurfacePoint>,
  std::vector<std::array<int, 2>>,
  std::vector<std::vector<SurfacePoint>>,
  std::vector<double>,
  std::vector<bool>
> subdivideSegments(
  ManifoldSurfaceMesh &mesh,
  VertexPositionGeometry &geometry,
  const std::vector<SurfacePoint> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<std::vector<SurfacePoint>> &segmentSurfacePoints,
  const std::vector<double> &segmentLengths,
  const std::vector<bool> &isFixedNode,
  const double &h
) {
  auto newNodes = nodes;
  auto newSegments = segments;
  auto newSegmentSurfacePoints = segmentSurfacePoints;
  auto newSegmentLengths = segmentLengths;
  auto newNodeIsFixed = isFixedNode;

  for (int i = 0; i < segments.size(); i++) {
    double edgeLen = segmentLengths[i];

    if (isFixedNode[segments[i][0]] && isFixedNode[segments[i][1]]) {
      continue;
    }

    if (edgeLen > 2 * h) {
      int divisionNum = std::ceil(edgeLen / (2 * h));
      double lenPerDivision = edgeLen / divisionNum;

      std::vector<SurfacePoint> segmentPoints = {};

      {
        segmentPoints.emplace_back(nodes[segments[i][0]]);
        for (auto sp : segmentSurfacePoints[i]) {
          segmentPoints.emplace_back(sp);
        }
        segmentPoints.emplace_back(nodes[segments[i][1]]);
      }
      
      auto cartesianCoords = surface_point_to_cartesian(mesh, geometry, segmentPoints);
      double len = .0;

      std::vector<SurfacePoint> newSurfacePoints = {};
      std::vector<std::vector<SurfacePoint>> newSegmentPoints = {{}};

      for (int j = 1; j < cartesianCoords.size(); j++) {
        double lenBefore = len;
        len += (cartesianCoords[j] - cartesianCoords[j - 1]).norm();

        int numBefore = std::floor(lenBefore / lenPerDivision);
        int num = std::floor(len / lenPerDivision);

        for (int k = 0; k < num - numBefore; k++) {
          int _n = numBefore + 1 + k;

          if (_n == divisionNum) {
            break;
          }

          double ratio = (_n * lenPerDivision - lenBefore) / (len - lenBefore);

          Face face = sharedFace(segmentPoints[j - 1], segmentPoints[j]);

          Vector3 vec0 = segmentPoints[j - 1].inFace(face).faceCoords;
          Vector3 vec1 = segmentPoints[j].inFace(face).faceCoords;

          // std::cout << "vec0 " << vec0 << ", vec1 " << vec1 << ", dvd: " << (vec0 * (1 - ratio) + vec1 * ratio) << ", t: " << ratio << std::endl;

          auto nsp = SurfacePoint(face, (vec0 * (1 - ratio) + vec1 * ratio));

          newSurfacePoints.emplace_back(nsp);
          newSegmentPoints.emplace_back(std::vector<SurfacePoint> {});
        }

        if (j != cartesianCoords.size() - 1) {
          newSegmentPoints[newSegmentPoints.size() - 1].emplace_back(segmentPoints[j]);
        }
      }

      assert(newSegmentPoints.size() == divisionNum);
      assert(newSurfacePoints.size() == divisionNum - 1);

      std::vector<int> newNodeIds = {};

      for (int j = 0; j < newSurfacePoints.size(); j++) {
        newNodes.emplace_back(newSurfacePoints[j]);
        newNodeIsFixed.emplace_back(false);
        newNodeIds.emplace_back(newNodes.size() - 1);
      }

      for (int j = 0; j <= newNodeIds.size(); j++) {
        if (j == 0) {
          newSegments[i][1] = newNodeIds[j];
          newSegmentLengths[i] = (edgeLen / divisionNum);
          newSegmentSurfacePoints[i] = newSegmentPoints[j];
        } else if (j == newNodeIds.size()) {
          newSegments.emplace_back(std::array<int, 2>{newNodeIds[j - 1], segments[i][1]});
          newSegmentLengths.emplace_back(edgeLen / divisionNum);
          newSegmentSurfacePoints.emplace_back(newSegmentPoints[j]);
        } else {
          newSegments.emplace_back(std::array<int, 2>{newNodeIds[j - 1], newNodeIds[j]});
          newSegmentLengths.emplace_back(edgeLen / divisionNum);
          newSegmentSurfacePoints.emplace_back(newSegmentPoints[j]);
        }
      }
    }
  }

  assert(newSegmentSurfacePoints.size() == newSegments.size());
  assert(newSegmentLengths.size() == newSegments.size());

  return {
    newNodes,
    newSegments,
    newSegmentSurfacePoints,
    newSegmentLengths,
    newNodeIsFixed
  };
}

std::tuple<
  std::vector<SurfacePoint>,
  std::vector<std::array<int, 2>>,
  std::vector<std::vector<SurfacePoint>>,
  std::vector<double>,
  std::vector<bool>
> remesh_curve_on_surface(
  ManifoldSurfaceMesh &mesh,
  VertexPositionGeometry &geometry,
  const std::vector<SurfacePoint> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<std::vector<SurfacePoint>> &segmentSurfacePoints,
  const std::vector<double> &segmentLengths,
  const std::vector<bool> &isFixedNode,
  const double &h
) {
  assert(segmentLengths.size() == segments.size());
  assert(segmentSurfacePoints.size() == segments.size());
  assert(isFixedNode.size() == nodes.size());

  // 1. remove nodes that are too close
  auto [newNodes, newSegments, newSegmentSurfacePoints, newSegmentLengths, newNodeIsFixed]
    = removeShortEdges(mesh, geometry, nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode, h);

  // 2. subdivide segments
  std::tie(newNodes, newSegments, newSegmentSurfacePoints, newSegmentLengths, newNodeIsFixed) = subdivideSegments(mesh, geometry, newNodes, newSegments, newSegmentSurfacePoints, newSegmentLengths, newNodeIsFixed, h);

  if (newNodes.size() < 3 || newSegments.size() < 3) {
    return {
      nodes,
      segments,
      segmentSurfacePoints,
      segmentLengths,
      isFixedNode
    };
  }

  return {
    newNodes,
    newSegments,
    newSegmentSurfacePoints,
    newSegmentLengths,
    newNodeIsFixed
  };
}
}
