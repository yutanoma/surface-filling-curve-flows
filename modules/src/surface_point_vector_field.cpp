#include <modules/surface_point_vector_field.h>

namespace modules {
std::tuple<Vector3> surface_point_vector_field(
  const VertexPositionGeometry &meshGeom,
  const VertexData<Vector2> &vectorField,
  const SurfacePoint &surfacePoint
) {
  assert(meshGeom.vertexNormals.size() == meshGeom.vertexPositions.size());
  assert(vectorField.size() == meshGeom.vertexPositions.size());

  switch (surfacePoint.type) {
    case SurfacePointType::Vertex:
      {
        auto v = surfacePoint.vertex;

        auto v_tan = vectorField[v];
        auto v_euc = v_tan.x * meshGeom.vertexTangentBasis[v][0] + v_tan.y * meshGeom.vertexTangentBasis[v][1];

        return {v_euc};
      }
      break;
    case SurfacePointType::Edge:
      {
        auto v0 = surfacePoint.edge.firstVertex();
        auto v1 = surfacePoint.edge.secondVertex();

        // todo: do more sophisticated interpolation?
        auto ratio = surfacePoint.tEdge;

        auto vec0 = vectorField[v0].x * meshGeom.vertexTangentBasis[v0][0] + vectorField[v0].y * meshGeom.vertexTangentBasis[v0][1];
        auto vec1 = vectorField[v1].x * meshGeom.vertexTangentBasis[v1][0] + vectorField[v1].y * meshGeom.vertexTangentBasis[v1][1];

        auto r0 = vec0.norm();
        auto r1 = vec1.norm();

        auto v = ((1 - ratio) * vec0 + ratio * vec1);
        if (norm(v) > 1e-6) {
          v = v.normalize();
        }

        auto r = (1 - ratio) * r0 + ratio * r1;

        return {r * v};
      }
      break;
    case SurfacePointType::Face:
      {
        auto v0 = surfacePoint.face.halfedge().vertex();
        auto v1 = surfacePoint.face.halfedge().next().vertex();
        auto v2 = surfacePoint.face.halfedge().next().next().vertex();
        
        auto vec0 = vectorField[v0].x * meshGeom.vertexTangentBasis[v0][0] + vectorField[v0].y * meshGeom.vertexTangentBasis[v0][1];
        auto vec1 = vectorField[v1].x * meshGeom.vertexTangentBasis[v1][0] + vectorField[v1].y * meshGeom.vertexTangentBasis[v1][1];
        auto vec2 = vectorField[v2].x * meshGeom.vertexTangentBasis[v2][0] + vectorField[v2].y * meshGeom.vertexTangentBasis[v2][1];

        auto r0 = vec0.norm();
        auto r1 = vec1.norm();
        auto r2 = vec2.norm();

        auto v = (surfacePoint.faceCoords.x * vec0) + (surfacePoint.faceCoords.y * vec1) + (surfacePoint.faceCoords.z * vec2);
        if (norm(v) > 1e-6) {
          v = v.normalize();
        }

        auto r = (surfacePoint.faceCoords.x * r0) + (surfacePoint.faceCoords.y * r1) + (surfacePoint.faceCoords.z * r2);

        return {r * v};
      }
      break;
    default:
      throw std::runtime_error("Invalid SurfacePointType");
  }
}
}
