#include <modules/surface_point_tangent_basis.h>

namespace modules {
  std::tuple<Vector3, Vector3> surface_point_tangent_basis(const VertexPositionGeometry &meshGeom, const SurfacePoint &surfacePoint) {
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
          auto n0 = meshGeom.vertexNormals[surfacePoint.edge.firstVertex()];
          auto n1 = meshGeom.vertexNormals[surfacePoint.edge.secondVertex()];

          auto ratio = surfacePoint.tEdge;
          auto n = normalize((1 - ratio) * n0 + ratio * n1);

          Vector3 _x = (meshGeom.vertexPositions[surfacePoint.edge.halfedge().tipVertex()] - meshGeom.vertexPositions[surfacePoint.edge.halfedge().tailVertex()]).normalize();

          auto x = cross(n, _x);
          auto y = cross(n, x);

          return {Vector3{x.x, x.y, x.z}, Vector3{y.x, y.y, y.z}};
        }
        break;
      case SurfacePointType::Face:
        {
          auto n0 = meshGeom.vertexNormals[surfacePoint.face.halfedge().vertex()];
          auto n1 = meshGeom.vertexNormals[surfacePoint.face.halfedge().next().vertex()];
          auto n2 = meshGeom.vertexNormals[surfacePoint.face.halfedge().next().next().vertex()];

          auto n = normalize(
              (surfacePoint.faceCoords.x * n0) + (surfacePoint.faceCoords.y * n1) + (surfacePoint.faceCoords.z * n2));

          auto _x = meshGeom.faceTangentBasis[surfacePoint.face][0];
          auto _y = meshGeom.faceTangentBasis[surfacePoint.face][1];

          auto x = cross(n, _x);
          auto y = cross(n, x);

          return {x.normalize(), y.normalize()};
        }
        break;
      default:
        throw std::runtime_error("Invalid SurfacePointType");
    }
  }
}
