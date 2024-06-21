#include <modules/get_tangent_basis.h>

namespace modules {
  std::tuple<Vector3, Vector3> get_tangent_basis(const VertexPositionGeometry &meshGeom, const SurfacePoint &surfacePoint) {
    switch (surfacePoint.type) {
      case SurfacePointType::Vertex:
        {
          auto x = meshGeom.vertexTangentBasis[surfacePoint.vertex][0];
          auto y = meshGeom.vertexTangentBasis[surfacePoint.vertex][1];
          return {x.normalize(), y.normalize()};
        }
        break;
      case SurfacePointType::Edge:
        {
          auto faces = surfacePoint.edge.adjacentFaces();
          Vector3 normal{.0, .0, .0};
          for (auto f : faces) {
            normal += meshGeom.faceNormal(f);
          }
          normal = normal.normalize();

          Vector3 x = (meshGeom.vertexPositions[surfacePoint.edge.halfedge().tipVertex()] - meshGeom.vertexPositions[surfacePoint.edge.halfedge().tailVertex()]).normalize();
          Vector3 y = cross(normal, x).normalize();

          assert(std::abs(1 - dot(normal, cross(x, y))) < .0001 && "Normal is not perpendicular to the tangent basis");

          return {Vector3{x.x, x.y, x.z}, Vector3{y.x, y.y, y.z}};
        }
        break;
      case SurfacePointType::Face:
        {
          auto x = meshGeom.faceTangentBasis[surfacePoint.face][0];
          auto y = meshGeom.faceTangentBasis[surfacePoint.face][1];
          return {x.normalize(), y.normalize()};
        }
        break;
      default:
        throw std::runtime_error("Invalid SurfacePointType");
    }
  }
}
