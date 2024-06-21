#include <modules/point_point_geodesic.h>

#include <modules/surface_point_to_cartesian.h>

namespace modules {
bool _hasSharedFace(
  const SurfacePoint &p0,
  const SurfacePoint &p1
) {
  auto f = sharedFace(p0, p1);

  // std::cout << f << ", " << Face() << std::endl;

  return f != Face();
}

Edge _sharedEdge(
  const SurfacePoint &p0,
  const SurfacePoint &p1
) {
  // only consider points on faces
  if (p0.type != SurfacePointType::Face || p1.type != SurfacePointType::Face) {
    return Edge();
  }

  auto f0 = p0.face, f1 = p1.face;
  for (auto e0 : f0.adjacentEdges()) {
    for (auto e1 : f1.adjacentEdges()) {
      if (e0 == e1) {
        return e0;
      }
    }
  }

  return Edge();
}

std::tuple<
  std::vector<SurfacePoint>
> point_point_geodesic(
  ManifoldSurfaceMesh &mesh,
  VertexPositionGeometry &geometry,
  GeodesicAlgorithmExact &mmp,
  const SurfacePoint &p0,
  const SurfacePoint &p1
) {
  // std::cout << "[" << p0.face() << ", " << p0.edge() << ", " << p0.vertex() << std::endl;

  if (_hasSharedFace(p0, p1)) {
    return {{p0, p1}};
  }

  auto e = _sharedEdge(p0, p1);
  bool hasEdge = e != Edge();

  if (hasEdge) {
    assert(p0.type == SurfacePointType::Face && p1.type == SurfacePointType::Face);

    auto v0 = geometry.vertexPositions[e.firstVertex()];
    auto v1 = geometry.vertexPositions[e.secondVertex()];

    auto cc = surface_point_to_cartesian(mesh, geometry, {p0, p1});
    auto c0 = cc[0], c1 = cc[1];

    double l0_0 = norm(c0 - v0);
    double l0_1 = norm(c0 - v1);
    double l1_0 = norm(c1 - v0);
    double l1_1 = norm(c1 - v1);
    double l = norm(v0 - v1);

    double theta0_0 = acos(dot(c0 - v0, v1 - v0) / (l0_0 * l));
    double theta0_1 = acos(dot(c0 - v1, v0 - v1) / (l0_1 * l));
    double theta1_0 = acos(dot(c1 - v0, v1 - v0) / (l1_0 * l));
    double theta1_1 = acos(dot(c1 - v1, v0 - v1) / (l1_1 * l));

    double theta_0 = theta0_0 + theta1_0;
    double theta_1 = theta0_1 + theta1_1;

    double area_tri1 = l0_0 * l1_0 * sin(theta_0) / 2.0;
    double area_tri2 = l0_1 * l1_1 * sin(theta_1) / 2.0;

    double ratio = area_tri1 / (area_tri1 + area_tri2);

    if (0 < ratio && ratio < 1) {
      return {
        {p0, SurfacePoint(e, ratio), p1}
      };
    }
  }
  // else if (hasSharedEdge(mesh, geometry, p0, p1)) {

  // }
  // else if (hasSharedVertex()) {
    
  // }

  // default
  mmp.propagate({p0}, GEODESIC_INF, {p1});
  std::vector<SurfacePoint> path = mmp.traceBack(p1);
  std::reverse(path.begin(), path.end());

  mmp.clear_data();

  return {
    path
  };
}
}
