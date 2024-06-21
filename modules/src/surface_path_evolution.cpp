#include "modules/surface_path_evolution.h"

#include "modules/get_tangent_basis.h"
#include "modules/check_intersection.h"
#include "modules/connect_surface_points.h"
#include "modules/remesh_curve_on_surface.h"

#include <chrono>

namespace modules {
const int max_iters = 100;
const double shrink = .8;

std::tuple<
  std::vector<SurfacePoint>,
  std::vector<std::vector<SurfacePoint>>
> get_surface_points(
  ManifoldSurfaceMesh& mesh,
  VertexPositionGeometry& meshGeom,
  const std::vector<SurfacePoint> &surfacePoints,
  const std::vector<TraceGeodesicResult> &tracePathResults,
  const std::vector<double> &tracePathLengths,
  const double _ratio
) {
  std::vector<SurfacePoint> newSurfacePoints = surfacePoints;
  std::vector<std::vector<SurfacePoint>> tracePaths(surfacePoints.size());

  for (int i = 0; i < tracePathResults.size(); i++) {
    auto res = tracePathResults[i];
    double _length = _ratio * tracePathLengths[i];

    double length = .0;

    std::vector<Vector3> pathPoints = surface_point_to_cartesian(mesh, meshGeom, res.pathPoints);

    for (int j = 1; j < pathPoints.size(); j++) {
      double currentEdgeLen = (pathPoints[j] - pathPoints[j - 1]).norm();
      length += currentEdgeLen;

      tracePaths[i].emplace_back(res.pathPoints[j-1]);

      if (length > _length) {
        double ratio = (_length - (length - currentEdgeLen)) / currentEdgeLen;

        Face face = sharedFace(res.pathPoints[j - 1], res.pathPoints[j]);

        Vector3 vec0 = res.pathPoints[j - 1].inFace(face).faceCoords;
        Vector3 vec1 = res.pathPoints[j].inFace(face).faceCoords;

        newSurfacePoints[i] = SurfacePoint(face, (vec0 * (1 - ratio) + vec1 * ratio));
        tracePaths[i].emplace_back(newSurfacePoints[i]);

        break;
      }

      if (j == pathPoints.size() - 1) {
        newSurfacePoints[i] = res.pathPoints[j];
        tracePaths[i].emplace_back(newSurfacePoints[i]);
      }
    }
  }

  return {newSurfacePoints, tracePaths};
}

std::tuple<
  std::vector<SurfacePoint> /* nodes */,
  std::vector<std::array<int, 2>> /* segments */,
  std::vector<std::vector<SurfacePoint>> /* segment paths */,
  std::vector<double> /* segment lengths */,
  std::vector<bool> /* isFixedNode */,
  std::vector<std::vector<SurfacePoint>> /* retraction paths */
> surface_path_evolution(
  ManifoldSurfaceMesh& mesh,
  VertexPositionGeometry& geometry,
  const std::vector<SurfacePoint> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<std::vector<SurfacePoint>> &segmentSurfacePoints,
  const std::vector<double> &segmentLengths,
  const std::vector<bool> &isFixedNode,
  const double h,
  const std::vector<Vector3> &direction
) {
  assert(direction.size() == nodes.size());
  assert(isFixedNode.size() == nodes.size());
  assert(segmentLengths.size() == segments.size());
  assert(segmentSurfacePoints.size() == segments.size());

  auto newNodes = nodes;
  auto newSegments = segments;
  auto newSegmentSurfacePoints = segmentSurfacePoints;
  auto newSegmentLengths = segmentLengths;
  auto newIsFixedNode = isFixedNode;

  std::vector<TraceGeodesicResult> traceResult(nodes.size());
  std::vector<double> traceLengths(nodes.size());
  std::vector<std::vector<SurfacePoint>> retractionPath(segments.size());
  
  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < nodes.size(); i++) {
    auto [x, y] = modules::get_tangent_basis(geometry, nodes[i]);
    auto d = direction[i];

    Vector2 traceVec {dot(d, x), dot(d, y)};
    TraceOptions opts = {
      includePath: true,
      maxIters: 1000
    };
    auto res = traceGeodesic(geometry, nodes[i], traceVec, opts);

    traceResult[i] = res;
    traceLengths[i] = norm(traceVec);
  }

  auto endTrace = std::chrono::high_resolution_clock::now();

  double step = 1.;

  for (int i = 0; i < max_iters; i++) {
    std::tie(newNodes, retractionPath) = get_surface_points(mesh, geometry, nodes, traceResult, traceLengths, step);

    assert(newNodes.size() == isFixedNode.size());

    auto endGetSurfacePoints = std::chrono::high_resolution_clock::now();

    std::tie(newSegmentSurfacePoints, newSegmentLengths) = modules::connect_surface_points(mesh, geometry, newNodes, segments);

    auto endConnectSurfacePoints = std::chrono::high_resolution_clock::now();

    std::tie(newNodes, newSegments, newSegmentSurfacePoints, newSegmentLengths, newIsFixedNode) = modules::remesh_curve_on_surface(mesh, geometry, newNodes, segments, newSegmentSurfacePoints, newSegmentLengths, isFixedNode, h);

    auto endRemesh = std::chrono::high_resolution_clock::now();

    // bool intersecting = false;
    bool intersecting = modules::check_intersection(mesh, geometry, newNodes, newSegments, newSegmentSurfacePoints, newSegmentLengths);

    auto endCheckIntersection = std::chrono::high_resolution_clock::now();

    std::cout << "path evolution: " << std::chrono::duration_cast<std::chrono::milliseconds>(endCheckIntersection - start).count() << "ms" << std::endl;
    std::cout << "  trace: " << std::chrono::duration_cast<std::chrono::milliseconds>(endTrace - start).count() << "ms" << std::endl;
    std::cout << "  get surface points: " << std::chrono::duration_cast<std::chrono::milliseconds>(endGetSurfacePoints - endTrace).count() << "ms" << std::endl;
    std::cout << "  connect surface points: " << std::chrono::duration_cast<std::chrono::milliseconds>(endConnectSurfacePoints - endGetSurfacePoints).count() << "ms" << std::endl;
    std::cout << "  remesh: " << std::chrono::duration_cast<std::chrono::milliseconds>(endRemesh - endConnectSurfacePoints).count() << "ms" << std::endl;
    std::cout << "  check intersection: " << std::chrono::duration_cast<std::chrono::milliseconds>(endCheckIntersection - endRemesh).count() << "ms" << std::endl;

    if (!intersecting) {
      return {
        newNodes,
        newSegments,
        newSegmentSurfacePoints,
        newSegmentLengths,
        newIsFixedNode,
        retractionPath
      };
    }

    step *= shrink;
  }
}
}
