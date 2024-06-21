#include "modules/check_intersection.h"

namespace modules {
  const double EPSILON = 1e-6;

  bool equals(double a, double b) {
    return std::abs(a - b) < EPSILON;
  }

  bool check_intersection(
    const ManifoldSurfaceMesh& mesh,
    const VertexPositionGeometry& meshGeom,
    const std::vector<SurfacePoint> &nodes,
    const std::vector<std::array<int, 2>> &segments,
    const std::vector<std::vector<SurfacePoint>> &edgeSurfacePoints,
    const std::vector<double> &edgeLengths
  ) {
    double aveEdgelen = 0;
    for (auto l : edgeLengths) {
      aveEdgelen += l;
    }
    aveEdgelen /= edgeLengths.size();

    std::map<int, std::vector<std::vector<SurfacePoint>>> faceToSurfacePoints;
    std::map<int, bool> faceRegistered;

    for (int i = 0; i < segments.size(); i++) {
      auto v0 = nodes[segments[i][0]], v1 = nodes[segments[i][1]];

      auto path = edgeSurfacePoints[i];

      std::vector<int> faceIds = {};
      std::vector<std::vector<SurfacePoint>> pointsOnCurrentFace = {};

      for (int j = 0; j < path.size() + 1; j++) {
        auto v_prev = j == 0 ? v0 : path[j - 1];
        auto v_next = j == path.size() ? v1 : path[j];

        // assert(checkAdjacent(v_prev, v_next));

        if (!checkAdjacent(v_prev, v_next)) {
          continue;
        }

        Face sharingFace = sharedFace(v_prev, v_next);

        int faceId = sharingFace.getIndex();

        if (faceIds.size() == 0 || faceId != faceIds.back()) {
          faceIds.push_back(faceId);
          pointsOnCurrentFace.push_back({v_prev});
        }

        pointsOnCurrentFace.back().push_back(v_next);
      }

      assert(faceIds.size() == pointsOnCurrentFace.size());

      for (int j = 0; j < faceIds.size(); j++) {
        int faceId = faceIds[j];

        if (!faceRegistered[faceId]) {
          faceRegistered[faceId] = true;
          faceToSurfacePoints[faceId] = {};
        }

        faceToSurfacePoints[faceId].push_back(pointsOnCurrentFace[j]);
      }
    }

    for (auto segments : faceToSurfacePoints) {
      auto faceId = segments.first;
      auto paths = segments.second;

      if (paths.size() < 2) {
        continue;
      }

      for (int i = 0; i < paths.size(); i++) {
        for (int j = 0; j < paths.size(); j++) {
          if (i == j) {
            continue;
          }

          std::array<std::vector<SurfacePoint>, 2> _paths = {paths[i], paths[j]};

          std::array<std::vector<Vector3>, 2> _cartesianPaths = {
            modules::surface_point_to_cartesian(mesh, meshGeom, _paths[0]),
            modules::surface_point_to_cartesian(mesh, meshGeom, _paths[1])
          };

          // check if the two paths intersect
          for (int l = 0; l < _cartesianPaths[0].size() - 1; l++) {
            for (int m = 0; m < _cartesianPaths[1].size() - 1; m++) {
              auto p0 = _cartesianPaths[0][l], p1 = _cartesianPaths[0][l + 1];
              auto q0 = _cartesianPaths[1][m], q1 = _cartesianPaths[1][m + 1];

              Vector3 d1 = p1 - p0, d2 = q1 - q0;

              Vector3 crossProduct = cross(d1, d2);

              if (crossProduct.norm() < EPSILON * aveEdgelen * aveEdgelen) {
                continue;
              }

              Vector3 diff = q0 - p0;
              double denominator = dot(crossProduct, crossProduct);
              double t1 = dot(cross(diff, d2), crossProduct) / denominator;
              double t2 = dot(cross(diff, d1), crossProduct) / denominator;

              if ((equals(t1, 0) || equals(t1, 1)) && (equals(t2, 0) || equals(t2, 1))) {
                // connected segments with their endpoints are allowed
                continue;
              }

              if (t1 > -EPSILON && t1 < 1+EPSILON && t2 > -EPSILON && t2 < 1+EPSILON) {
                std::cout << "intersected! t1: " << t1 << ", t2: " << t2 << ", face: " << faceId << std::endl;
                return true;  // The line segments intersect
              }
            }
          }
        }
      }
    }

    return false;
  }
}
