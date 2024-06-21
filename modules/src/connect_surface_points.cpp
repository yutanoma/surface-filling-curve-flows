#include <modules/connect_surface_points.h>

#include <modules/surface_point_to_cartesian.h>
#include <modules/point_point_geodesic.h>

namespace modules {
  std::tuple<
    std::vector<std::vector<SurfacePoint>>,
    std::vector<double>
  > connect_surface_points(
    ManifoldSurfaceMesh& mesh,
    VertexPositionGeometry& geometry,
    const std::vector<SurfacePoint> &nodes,
    const std::vector<std::array<int, 2>> &segments
  ) {
    int numSegments = segments.size();

    std::vector<std::vector<SurfacePoint>> edgeSurfacePoints;
    edgeSurfacePoints.resize(numSegments, {});

    std::vector<double> edgeLengths(numSegments, 0.0);

    GeodesicAlgorithmExact mmp(mesh, geometry);

    for (int i = 0; i < numSegments; i++) {
      auto p0 = nodes[segments[i][0]], p1 = nodes[segments[i][1]];

      auto [path] = point_point_geodesic(mesh, geometry, mmp, p0, p1);

      // mmp.clear_data();
      // mmp.propagate({p0}, GEODESIC_INF, {p1});
      // std::vector<SurfacePoint> path = mmp.traceBack(p1);
      // std::reverse(path.begin(), path.end());

      double length = .0;

      auto cartesianCoord = surface_point_to_cartesian(mesh, geometry, path);

      for (int j = 0; j < path.size(); j++) {
        if (j != 0) {
          length += (cartesianCoord[j] - cartesianCoord[j - 1]).norm();
        }

        if (j != 0 && j != path.size() - 1) {
          edgeSurfacePoints[i].emplace_back(path[j]);
        }
      }

      edgeLengths[i] = length;
    }

    assert(edgeSurfacePoints.size() == numSegments);
    assert(edgeLengths.size() == numSegments);

    return {edgeSurfacePoints, edgeLengths};
  }

}