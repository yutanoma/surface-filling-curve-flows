#include <modules/surface_point_value.h>

namespace modules {
double surface_point_value(
  const VertexPositionGeometry &meshGeom,
  const VertexData<double> &function,
  const SurfacePoint &surfacePoint
) {
  assert(meshGeom.vertexNormals.size() == meshGeom.vertexPositions.size());

  switch (surfacePoint.type) {
    case SurfacePointType::Vertex:
      {
        return function[surfacePoint.vertex];
      }
      break;
    case SurfacePointType::Edge:
      {
        auto v0 = surfacePoint.edge.firstVertex();
        auto v1 = surfacePoint.edge.secondVertex();

        auto ratio = surfacePoint.tEdge;

        auto r0 = function[v0];
        auto r1 = function[v1];

        return (1 - ratio) * r0 + ratio * r1;
      }
      break;
    case SurfacePointType::Face:
      {
        auto v0 = surfacePoint.face.halfedge().vertex();
        auto v1 = surfacePoint.face.halfedge().next().vertex();
        auto v2 = surfacePoint.face.halfedge().next().next().vertex();
        
        auto r0 = function[v0];
        auto r1 = function[v1];
        auto r2 = function[v2];

        return surfacePoint.faceCoords.x * r0 + surfacePoint.faceCoords.y * r1 + surfacePoint.faceCoords.z * r2;
      }
      break;
    default:
      throw std::runtime_error("Invalid SurfacePointType");
  }
}
}
