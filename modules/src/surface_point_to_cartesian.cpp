#include <modules/surface_point_to_cartesian.h>

namespace modules {
  std::vector<Vector3> surface_point_to_cartesian(
    const ManifoldSurfaceMesh& mesh,
    const VertexPositionGeometry& meshGeom,
    const std::vector<SurfacePoint> &surfacePoints
  ) {
    std::vector<Vector3> cartesianPoints(surfacePoints.size());

    for (unsigned long i = 0; i < surfacePoints.size(); i++) {
      auto surfacePoint = surfacePoints[i];

      switch (surfacePoint.type) {
        case SurfacePointType::Vertex:
          {
            cartesianPoints[i] = meshGeom.vertexPositions[surfacePoint.vertex];
          }
          break;
        case SurfacePointType::Edge:
          {
            auto v0 = meshGeom.vertexPositions[surfacePoint.edge.halfedge().vertex()];
            auto v1 = meshGeom.vertexPositions[surfacePoint.edge.halfedge().twin().vertex()];
            double ratio = surfacePoint.tEdge;
            cartesianPoints[i] = (1 - ratio) * v0 + ratio * v1;
          }
          break;
        case SurfacePointType::Face:
          {
            auto faceCoords = surfacePoint.faceCoords;
            auto v0 = meshGeom.vertexPositions[surfacePoint.face.halfedge().vertex()];
            auto v1 = meshGeom.vertexPositions[surfacePoint.face.halfedge().next().vertex()];
            auto v2 = meshGeom.vertexPositions[surfacePoint.face.halfedge().next().next().vertex()];

            cartesianPoints[i] = (faceCoords.x * v0) + (faceCoords.y * v1) + (faceCoords.z * v2);
          }
          break;
        default:
          throw std::runtime_error("Invalid SurfacePointType");
      }

      // std::cout << cartesianPoints[i] << std::endl;
    }

    return cartesianPoints;
  }
}
