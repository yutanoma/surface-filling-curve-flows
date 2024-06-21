#include "modules/gc_igl_convert.h"

#include "geometrycentral/surface/surface_mesh_factories.h"

namespace modules {
  std::tuple<
    Eigen::MatrixXd, // V
    Eigen::MatrixXi, // F
    Eigen::MatrixXi, // EF
    Eigen::MatrixXi, // FE
    Eigen::MatrixXi // EV
  > gc_to_igl(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry
  ) {
    Eigen::MatrixXd V(mesh.nVertices(), 3);
    Eigen::MatrixXi F(mesh.nFaces(), 3);

    for (Vertex v : mesh.vertices()) {
      for (int i = 0; i < 3; i++) {
        V(v.getIndex(), i) = geometry.vertexPositions[v][i];
      }
    }

    for (Face f : mesh.faces()) {
      int i = 0;
      for (Vertex v : f.adjacentVertices()) {
        F(f.getIndex(), i) = v.getIndex();
        i++;
      }
    }

    Eigen::MatrixXi EV(mesh.nEdges(), 2);
    for (Edge e : mesh.edges()) {
      EV(e.getIndex(), 0) = e.firstVertex().getIndex();
      EV(e.getIndex(), 1) = e.secondVertex().getIndex();
    }

    Eigen::MatrixXi EF(mesh.nEdges(), 2);
    for (Edge e : mesh.edges()) {
      int i = 0;
      for (Face f : e.adjacentFaces()) {
        EF(e.getIndex(), i) = f.getIndex();
        i++;
      }
    }

    Eigen::MatrixXi FE(mesh.nFaces(), 3);
    for (Face f : mesh.faces()) {
      int i = 0;
      for (Edge e : f.adjacentEdges()) {
        FE(f.getIndex(), i) = e.getIndex();
        i++;
      }
    }

    return std::make_tuple(V, F, EF, FE, EV);
  }

  std::tuple<
    Eigen::MatrixXd, // V
    Eigen::MatrixXi // F
  > gc_to_igl_vf(
    ManifoldSurfaceMesh &mesh,
    VertexPositionGeometry &geometry
  ) {
    auto [V, F, EF, FE, EV] = gc_to_igl(mesh, geometry);

    return std::make_tuple(V, F);
  }
}
