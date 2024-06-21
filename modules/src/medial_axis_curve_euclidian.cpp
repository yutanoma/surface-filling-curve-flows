#include <modules/medial_axis_curve_euclidian.h>

#include <modules/circumcenter.h>
#include <modules/evolution_mode.h>
#include <modules/surface_point_to_cartesian.h>
#include <modules/get_tangent_basis.h>

#include <igl/copyleft/tetgen/cdt.h>
#include <igl/tet_tet_adjacency.h>
#include <igl/barycenter.h>
#include <igl/winding_number.h>
#include <igl/is_edge_manifold.h>

namespace modules {
std::tuple<
  std::vector<geometrycentral::Vector3>, // node list on medial axis
  std::vector<std::array<int, 2>>, // segment list on medial axis
  std::vector<double>, // radius list on medial axis
  std::vector<std::vector<std::vector<int>>>, // node list to the nearest medial axis node
  Eigen::MatrixXd, // V of underlying mesh
  Eigen::MatrixXi, // F of underlying mesh
  Eigen::MatrixXi // T of underlying mesh
> get_medial_axis_curve_euclidian(
  geometrycentral::surface::ManifoldSurfaceMesh &mesh,
  geometrycentral::surface::VertexPositionGeometry &geometry,
  const Eigen::MatrixXd &meshV,
  const Eigen::MatrixXi &meshF,
  const std::vector<std::vector<geometrycentral::surface::SurfacePoint>> &curves
) {
  geometry.requireFaceTangentBasis();
  geometry.requireFaceNormals();

  int numVertices = 0;
  std::vector<std::vector<int>> curveNodeList(curves.size(), std::vector<int> {});
  for (int i = 0; i < curves.size(); i++) {
    auto curve = curves[i];
    numVertices += curve.size();
    curveNodeList[i].resize(curve.size(), -1);
  }

  // 1. do delaunay triangulation
  Eigen::MatrixXi V2curveNode(numVertices, 2), _E(numVertices, 2);
  Eigen::MatrixXd T_V(numVertices, 3);

  std::vector<std::vector<Vector3>> cartesianCoords(curves.size());

  int _numVertices = 0;
  for (int i = 0; i < curves.size(); i++) {
    int firstIdx = _numVertices;

    cartesianCoords[i] = modules::surface_point_to_cartesian(mesh, geometry, curves[i]);

    for (int j = 0; j < curves[i].size(); j++) {
      auto p = cartesianCoords[i][j];
      T_V.row(_numVertices) << p.x, p.y, p.z;

      V2curveNode(_numVertices, 0) = i;
      V2curveNode(_numVertices, 1) = j;
      curveNodeList[i][j] = _numVertices;

      _numVertices++;
    }
  }

  // std::cout << _E << std::endl;

  Eigen::MatrixXd _T_V;
  Eigen::MatrixXi _T, _, _TF;

  igl::copyleft::tetgen::CDTParam param;
  param.flags = "S0";

  igl::copyleft::tetgen::cdt(T_V, _, param, _T_V, _T, _TF);

  // std::cout << TF << std::endl;

  assert(T_V.rows() == _T_V.rows());

  Eigen::MatrixXi CT;
  Eigen::MatrixXi F;

  // 0. just return the tet mesh
  {
    // CT = _T;
    // F = _TF;
  }

  // 1. judge using winding numbers
  {
    Eigen::MatrixXd BC;
    Eigen::VectorXd W;
    std::cout << "l85" << std::endl;
    igl::barycenter(T_V, _T, BC);
    std::cout << "l87" << std::endl;
    igl::winding_number(meshV, meshF, BC, W);
    std::cout << "l89" << std::endl;

    // Extract interior tets
    CT = Eigen::MatrixXi((W.array()>0.5).count(),4);
    {
      size_t k = 0;
      for(size_t t = 0;t<_T.rows();t++)
      {
        if(W(t)>0.5)
        {
          CT.row(k) = _T.row(t);
          k++;
        }
      }
    }
    // find bounary facets of interior tets
    igl::boundary_facets(CT, F);
  }

  // 2. judge using barycenter and normal
  {
    // std::vector<int> insideTets = {};

    // for (int i = 0; i < _T.rows(); i++) {
    //   int tid = i;

    //   Vector3 barycenter {.0, .0, .0};
    //   for (int j = 0; j < 4; j++) {
    //     int v = _T(tid, j);
    //     barycenter += Vector3 {T_V(v, 0), T_V(v, 1), T_V(v, 2)};
    //   }
    //   barycenter /= 4.0;

    //   std::array<int, 4> isInside = {false, false, false, false};
    //   int insideNum = 0;

    //   for (int j = 0; j < 4; j++) {
    //     int v = _T(tid, j);
    //     int curveId = V2curveNode(v, 0), curveNodeId = V2curveNode(v, 1);
    //     auto sp = curves[curveId][curveNodeId];
    //     auto p = cartesianCoords[curveId][curveNodeId];
    //     auto [x, y] = modules::get_tangent_basis(geometry, sp);
    //     auto normal = cross(x, y);

    //     std::cout << normal << ", " << barycenter - p << ", " << dot(normal, barycenter - p) << std::endl;

    //     bool isOutside = dot(normal, barycenter - p) > -__DBL_EPSILON__;

    //     isInside[j] = !isOutside;
    //     if (isInside[j]) {
    //       insideNum++;
    //     }
    //   }

    //   // note: maybe this is unnecessary?
    //   // assert(isInside[0] == isInside[1] && isInside[1] == isInside[2] && isInside[2] == isInside[3]);

    //   if (insideNum > 3) {
    //     insideTets.emplace_back(tid);
    //   }
    // }

    // CT.resize(insideTets.size(), 4);
    // for (int i = 0; i < insideTets.size(); i++) {
    //   CT.row(i) = _T.row(insideTets[i]);
    // }

    // igl::boundary_facets(CT, F);
  }

  // assert(igl::is_edge_manifold(F));

  return {
    {},
    {},
    {},
    {},
    T_V,
    F,
    CT
  };
  }
}
