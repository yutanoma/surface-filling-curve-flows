#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/direction_fields.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/volume_mesh.h>

#include <igl/predicates/delaunay_triangulation.h>
#include <igl/triangle/cdt.h>
#include <igl/edge_topology.h>
#include <igl/PI.h>
#include <igl/copyleft/tetgen/cdt.h>
#include <igl/copyleft/tetgen/mesh_to_tetgenio.h>
#include <igl/writeOBJ.h>
#include <igl/circulation.h>

#include <tetgen-src/tetgen.h>
// #include <igl/copyleft/cgal/delaunay_triangulation.h>

#include <queue>
#include <chrono>

#include "args/args.hxx"
#include "imgui.h"

#include "modules/scene_file.h"
#include "modules/curve_evolution.h"
#include "modules/medial_axis_2d.h"
#include "modules/evolution_mode.h"
#include "modules/curve_remesh_2d.h"
#include "modules/circumcenter.h"
#include "modules/surface_point_to_cartesian.h"
#include "modules/read_curve_on_mesh.h"
#include "modules/connect_surface_points.h"
#include "modules/remesh_curve_on_surface.h"
#include "modules/surface_filling_energy_geodesic.h"
#include "modules/surface_path_evolution.h"
#include "modules/get_tangent_basis.h"
#include "modules/write_curve.h"
#include "modules/curve_to_arclength.h"
#include "modules/cut_mesh_with_curve.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

using namespace modules;

float radius = 30;
float timestep = 1;
float h = igl::PI * radius / 25;
float rmax = 0.5;
// double branchRatio = 1.18;

// Polyscope visualization handle, to quickly add data to the surface
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

modules::SceneObject scene;

std::vector<SurfacePoint> nodes = {}, restNodes = {}, initialNodes = {};
std::vector<bool> isFixedNode = {}, isInitialFixedNode = {};
std::vector<std::array<int, 2>> segments = {}, restSegments = {}, initialSegments = {};

std::vector<std::vector<SurfacePoint>> segmentSurfacePoints = {}, restSegmentSurfacePoints = {}, initialSegmentSurfacePoints = {};
std::vector<double> segmentLengths = {}, restSegmentLengths = {}, initialSegmentLengths = {};

VertexData<Vector2> vField;
VertexData<double> smoothedFunction;

// Example computation function -- this one computes and registers a scalar
// quantity
int iteration = 0;
bool writeCurve = false;
bool writeData = false;
bool runLoop = false;

std::ofstream data_out;

std::tuple<
  std::vector<Vector3> /* nodes */,
  std::vector<std::array<int, 2>> /* segments */,
  std::vector<Vector3> /* cartesianCoords of the nodes */
> formatGeodesicCurve(
  ManifoldSurfaceMesh &_mesh,
  VertexPositionGeometry &_geometry,
  const std::vector<SurfacePoint> &_nodes,
  const std::vector<std::array<int, 2>> &_segments,
  const std::vector<std::vector<SurfacePoint>> &_segmentSurfacePoints
) {
  std::vector<Vector3> cartesianCoords = modules::surface_point_to_cartesian(_mesh, _geometry, _nodes);

  std::vector<Vector3> allNodes = cartesianCoords;
  std::vector<std::array<int, 2>> allSegments = {};

  for (int i = 0; i < _segments.size(); i++) {
    int v0 = _segments[i][0], v1 = _segments[i][1];

    auto cc = surface_point_to_cartesian(_mesh, _geometry, _segmentSurfacePoints[i]);

    std::vector<int> ccIdx = {};
    ccIdx.emplace_back(v0);
    for (int j = 0; j < cc.size(); j++) {
      allNodes.emplace_back(cc[j]);
      ccIdx.emplace_back(allNodes.size() - 1);
    }
    ccIdx.emplace_back(v1);

    for (int j = 0; j < ccIdx.size() - 1; j++) {
      std::array<int, 2> line = {ccIdx[j], ccIdx[j + 1]};
      allSegments.push_back(line);
    }
  }

  return {
    allNodes,
    allSegments,
    cartesianCoords
  };
}

std::tuple<
  std::vector<Vector3> /* nodes */,
  std::vector<std::array<int, 2>> /* segments */,
  std::vector<double> /* norm of the vector */
> formatVector(
  std::vector<Vector3> &cartesianCoords,
  std::vector<Vector3> &vector
) {
  assert(cartesianCoords.size() == vector.size());

  std::vector<Vector3> dNodes = cartesianCoords;
  std::vector<std::array<int, 2>> dSegments = {};
  std::vector<double> dNorm = {};

  for (int i = 0; i < vector.size(); i++) {
    dNodes.emplace_back(cartesianCoords[i] + vector[i]);
    dNorm.emplace_back(norm(vector[i]));
  }

  for (int i = 0; i < vector.size(); i++) {
    int oi = i + vector.size();
    std::array<int, 2> line = {i, oi};
    dSegments.push_back(line);
  }

  return {
    dNodes,
    dSegments,
    dNorm
  };
}

std::tuple<
  std::vector<Vector3> /* nodes */,
  std::vector<std::array<int, 2>> /* segments */
> formatGeodesicPaths(
  ManifoldSurfaceMesh &_mesh,
  VertexPositionGeometry &_geometry,
  std::vector<std::vector<SurfacePoint>> &paths
) {
  std::vector<Vector3> nodes = {};
  std::vector<std::array<int, 2>> segments = {};

  for (int i = 0; i < paths.size(); i++) {
    auto path = paths[i];

    if (path.size() == 0) {
      continue;
    }

    auto cc = surface_point_to_cartesian(_mesh, _geometry, path);

    for (int j = 0; j < path.size() - 1; j++) {
      int v0 = nodes.size();
      nodes.emplace_back(cc[j]);
      int v1 = nodes.size();
      nodes.emplace_back(cc[j + 1]);

      std::array<int, 2> line = {v0, v1};
      segments.push_back(line);
    }
  }

  return {
    nodes,
    segments
  };
}

std::tuple<
  std::vector<Vector3> /* nodes */,
  std::vector<double> /* radius */
> formatMedialAxis(
  const std::vector<std::vector<Vector3>> &medialAxis,
  const std::vector<Vector3> &cartesianCoords
){
  std::vector<Vector3> maNodes = {};
  std::vector<double> maRadius = {};

  for (int i = 0; i < medialAxis.size(); i++) {
    assert(medialAxis[i].size() == 2);

    for (int j = 0; j < medialAxis[i].size(); j++) {
      maNodes.emplace_back(medialAxis[i][j]);

      double r = norm(medialAxis[i][j] - cartesianCoords[i]);
      maRadius.emplace_back(r);
    }
  }

  return {
    maNodes,
    maRadius
  };
}

std::tuple<
  Eigen::MatrixXd /* V */,
  Eigen::MatrixXi /* T */
> tetrahedralization(
  ManifoldSurfaceMesh &_mesh,
  VertexPositionGeometry &_geometry,
  std::vector<Vector3> &cartesianPaths
) {
  Eigen::MatrixXd V(cartesianPaths.size(), 3);

  for (int i = 0; i < cartesianPaths.size(); i++) {
    V(i, 0) = cartesianPaths[i].x;
    V(i, 1) = cartesianPaths[i].y;
    V(i, 2) = cartesianPaths[i].z;
  }

  auto start = std::chrono::high_resolution_clock::now();

  Eigen::MatrixXi T;

  {
    // Eigen::MatrixXd T_V;
    // Eigen::MatrixXi _, __;
    // std::string flags = "S0cQ";
    // igl::copyleft::tetgen::tetrahedralize(V, _, flags, T_V, T, __);
  }

  {
    tetgenio in, out;

    in.firstnumber = 0;

    in.numberofpoints = V.rows();
    in.pointlist = new REAL[in.numberofpoints * 3];
    // loop over points
    for(int i = 0; i < V.rows(); i++)
    {
      in.pointlist[i*3+0] = V(i, 0);
      in.pointlist[i*3+1] = V(i, 1);
      in.pointlist[i*3+2] = V(i, 2);
    }
    in.numberoffacets = 0;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    char* flags = "S0cQ";

    tetrahedralize(flags, &in, &out);

    assert(out.numberofpoints == V.rows());
    assert(out.numberofcorners == 4);

    T.resize(out.numberoftetrahedra, 4);
    int min_index = 1e7;
    int max_index = -1e7;
    // loop over tetrahedra
    for(int i = 0; i < out.numberoftetrahedra; i++)
    {
      for(int j = 0; j<out.numberofcorners; j++)
      {
        int index = out.tetrahedronlist[i * out.numberofcorners + j];
        T(i, j) = index;
        min_index = (min_index > index ? index : min_index);
        max_index = (max_index < index ? index : max_index);
      }
    }
  }

  auto end = std::chrono::high_resolution_clock::now();

  std::cout << "tetrahedralize: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;

  return {
    V,
    T
  };
}

void doWork() {
  iteration++;

  restNodes = nodes;
  restSegments = segments;
  restSegmentSurfacePoints = segmentSurfacePoints;
  restSegmentLengths = segmentLengths;

  std::cout << "===== iteration: " << iteration << "=====" << std::endl;

  std::cout << "numNodes: " << nodes.size() << std::endl;
  std::cout << "numSegments: " << segments.size() << std::endl;
  std::cout << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  auto cartesianCoords = modules::surface_point_to_cartesian(*mesh, *geometry, nodes);

  // 1. get the descent direction
  modules::SurfaceFillingEnergy::Options options = {
    p: scene.p,
    q: scene.q,
    w_fieldAlignedness: scene.w_fieldAlignedness,
    vectorField: vField,
    w_curvatureAlignedness: scene.w_curvatureAlignedness,
    w_bilaplacian: scene.w_bilaplacian,
    useAnisotropicAlphaOnMesh: scene.varyingAlpha,
    alphaRatioOnMesh: smoothedFunction,
    useGeodesicMedialAxis: scene.useGeodesicMedialAxis
  };
  auto [d, g, f, medialAxis] = modules::surface_filling_energy_geodesic(*mesh, *geometry, nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode, cartesianCoords, radius, rmax, options);

  // 2. evolve the curve but be careful not to cause self-intersection
  std::vector<std::vector<SurfacePoint>> retractionPaths = {};
  std::tie(nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode, retractionPaths) = modules::surface_path_evolution(*mesh, *geometry, nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode, h, d);

  auto end = std::chrono::high_resolution_clock::now();

  int time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  std::cout << "total time: " << time << "ms" << std::endl;
  std::cout << "===== iteration end =====" << std::endl << std::endl;

  auto [allNodes, allSegments, newCartesianCoords] = formatGeodesicCurve(*mesh, *geometry, nodes, segments, segmentSurfacePoints);

  auto crv = polyscope::registerCurveNetwork("curve", allNodes, allSegments);

  // visualize
  {
    auto [prevNodes, prevSegments, prevCartesianCoords] = formatGeodesicCurve(*mesh, *geometry, restNodes, restSegments, restSegmentSurfacePoints);
    polyscope::registerCurveNetwork("prev curve", prevNodes, prevSegments)
      ->setEnabled(false)
    ;

    auto [allNodes, allSegments, newCartesianCoords] = formatGeodesicCurve(*mesh, *geometry, nodes, segments, segmentSurfacePoints);

    auto arclength = modules::curve_to_arclength(allNodes, allSegments);

    crv->addNodeScalarQuantity("arclength", arclength);

    // polyscope::registerPointCloud("initial nodes", allNodes);

    std::vector<Vector3> fixedNodes = {}, freeNodes = {};
    for (int i = 0; i < nodes.size(); i++) {
      if (isFixedNode[i]) {
        fixedNodes.emplace_back(newCartesianCoords[i]);
      } else {
        freeNodes.emplace_back(newCartesianCoords[i]);
      }
    }

    polyscope::registerPointCloud("free nodes", freeNodes)
      // ->setPointColor({1.0, 1.0, 1.0})
    ;
    polyscope::registerPointCloud("fixed nodes", fixedNodes)
      // ->setPointColor({.5, .5, .5})
    ;

    auto [dNodes, dSegments, dNorm] = formatVector(cartesianCoords, d);
    auto [gNodes, gSegments, gNorm] = formatVector(cartesianCoords, g);

    auto dn = polyscope::registerCurveNetwork("descent direction", dNodes, dSegments);
    dn->setEnabled(false);
    dn->addEdgeScalarQuantity("norm", dNorm);

    // auto gn = polyscope::registerCurveNetwork("gradient direction", gNodes, gSegments);
    // gn->setEnabled(false);
    // gn->addEdgeScalarQuantity("norm", gNorm);

    auto [maNodes, maRadius] = formatMedialAxis(medialAxis, cartesianCoords);
    auto pma = polyscope::registerPointCloud("medial axis", maNodes);
    pma->setEnabled(false);
    pma->addScalarQuantity("radius", maRadius);

    auto [rNodes, rSegments] = formatGeodesicPaths(*mesh, *geometry, retractionPaths);
    polyscope::registerCurveNetwork("retraction paths", rNodes, rSegments)
      ->setEnabled(false)
    ;

    std::vector<Vector3> normalDirections, tangentX, tangentY;
    for (auto node : restNodes) {
      auto [x, y] = modules::get_tangent_basis(*geometry, node);
      tangentX.emplace_back(x / 30);
      tangentY.emplace_back(y / 30);
      normalDirections.emplace_back(cross(x, y) / 30);
    }

    auto [nn, ns, _] = formatVector(cartesianCoords, normalDirections);
    polyscope::registerCurveNetwork("normal directions", nn, ns)->setEnabled(false);

    auto [txn, txs, __] = formatVector(cartesianCoords, tangentX);
    polyscope::registerCurveNetwork("tangent x", txn, txs)->setEnabled(false);

    auto [tyn, tys, ___] = formatVector(cartesianCoords, tangentY);
    polyscope::registerCurveNetwork("tangent y", tyn, tys)->setEnabled(false);

    auto finishVisualize = std::chrono::high_resolution_clock::now();

    std::cout << "visualize time: " << std::chrono::duration_cast<std::chrono::milliseconds>(finishVisualize - end).count() << "ms" << std::endl;

    if (writeCurve) {
      std::string fname = "./objs/curve_" + std::to_string(iteration) + ".obj";
      modules::write_curve(fname, allNodes, allSegments);

      std::string radiusfilename = "./radii/radii_" + std::to_string(iteration - 1) + ".csv";
      std::ofstream radiusfile(radiusfilename);

      radiusfile << "radius, weight, alpha" << std::endl;

      std::vector<double> nodeWeight(cartesianCoords.size(), 0);
      for (int i = 0; i < restSegments.size(); i++) {
        auto [v0, v1] = restSegments[i];
        nodeWeight[v0] += restSegmentLengths[i] / 2;
        nodeWeight[v1] += restSegmentLengths[i] / 2;
      }

      double branchRatio = std::sqrt(std::sqrt(2));
      double branchRadius = radius * branchRatio;
      double alpha = 4 / (branchRadius * branchRadius);

      for (int i = 0; i < cartesianCoords.size(); i++) {
        auto v = cartesianCoords[i];
        for (int j = 0; j < 2; j++) {
          auto c = medialAxis[i][j];
          double radius = norm(v - c);
          radiusfile << radius << ", " << nodeWeight[i] << ", " << alpha << std::endl;
        }
      }

      std::string _fname = "./fixed_points/fixed_points_" + std::to_string(iteration) + ".obj";
      std::ofstream _myfile(_fname);

      for (auto v : fixedNodes) {
        _myfile << "v " << v.x << " " << v.y << " " << v.z << std::endl;
      }

      // std::string _fname2 = "./medial_axis/medial_axis_" + std::to_string(iteration) + ".obj";
      // std::ofstream _myfile2(_fname2);

      // for (auto v : maNodes) {
      //   _myfile2 << "v " << v.x << " " << v.y << " " << v.z << std::endl;
      // }

      // std::string _fname3 = "./descent/descent_" + std::to_string(iteration) + ".obj";
      // std::ofstream _myfile3(_fname3);

      // for (auto d : dNodes) {
      //   _myfile3 << "v " << d.x << " " << d.y << " " << d.z << std::endl;
      // }
      // for (auto d : dSegments) {
      //   _myfile3 << "l " << d[0]+1 << " " << d[1]+1 << std::endl;
      // }

      // std::string _fname4 = "./trace_path/trace_path_" + std::to_string(iteration) + ".obj";
      // std::ofstream _myfile4(_fname4);

      // for (auto d : rNodes) {
      //   _myfile4 << "v " << d.x << " " << d.y << " " << d.z << std::endl;
      // }
      // for (auto d : rSegments) {
      //   _myfile4 << "l " << d[0]+1 << " " << d[1]+1 << std::endl;
      // }

      // auto [cutVs, cutFs, _V, _F] = modules::cut_mesh_with_curve(*mesh, *geometry, nodes, segments, segmentSurfacePoints);
      // for (int i = 0; i < cutVs.size(); i++) {
      //   std::string _fname = "./cut_mesh/cut_mesh_" + std::to_string(iteration) + "_" + std::to_string(i) + ".obj";
      //   std::ofstream _myfile4(_fname);

      //   auto _V = cutVs[i];
      //   auto _F = cutFs[i];

      //   for (int j = 0; j < _V.rows(); j++) {
      //     _myfile4 << "v " << _V(j, 0) << " " << _V(j, 1) << " " << _V(j, 2) << std::endl;
      //   }

      //   for (int j = 0; j < _F.rows(); j++) {
      //     _myfile4 << "f " << _F(j, 0) + 1 << " " << _F(j, 1) + 1 << " " << _F(j, 2) + 1 << std::endl;
      //   }
      // }
    }

    if (writeData) {
      Eigen::VectorXd _d(d.size() * 3);
      _d.setZero();
      for (int i = 0; i < d.size(); i++) {
        _d(i * 3 + 0) = d[i].x;
        _d(i * 3 + 1) = d[i].y;
        _d(i * 3 + 2) = d[i].z;
      }

      double l1 = _d.lpNorm<1>();
      double l2 = _d.lpNorm<2>();
      double linf = _d.lpNorm<Eigen::Infinity>();

      // iteration, time, numNodes, f, descent norm (L2), descent norm (L1), descent norm (L∞)
      data_out << iteration << ", " << time << ", " << nodes.size() << ", " << f << ", " << l2 << ", " << l1 << ", " << linf << std::endl;
    }
  }
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
  ImGui::Checkbox("write obj", &writeCurve);
  ImGui::Checkbox("write data", &writeData);
  ImGui::Checkbox("run loop", &runLoop);

  ImGui::SliderFloat("radius", &radius, scene.radius / 100, scene.radius * 100, "%.3f", ImGuiSliderFlags_Logarithmic);
  ImGui::SliderFloat("rmax", &rmax, scene.rmax / 100, scene.rmax * 100, "%.3f", ImGuiSliderFlags_Logarithmic);

  if (ImGui::IsKeyDown(' ') || runLoop) {
    doWork();
  }

  if (ImGui::Button("cut mesh with curve (expermental)")) {
    auto [cutVs, cutFs, _V, _F] = modules::cut_mesh_with_curve(*mesh, *geometry, nodes, segments, segmentSurfacePoints);
    polyscope::registerSurfaceMesh("cut mesh", _V, _F);
    for (int i = 0; i < cutVs.size(); i++) {
      std::cout << "cut mesh " << i << ": " << cutVs[i].rows() << ", " << cutFs[i].rows() << std::endl;
      polyscope::registerSurfaceMesh("cut mesh " + std::to_string(i), cutVs[i], cutFs[i]);
    }
  }

  if (iteration % 100 == 0) {
    runLoop = false;
  }
}

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

  // Parse args
  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help &h) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  polyscope::init();

  // Set the callback function
  polyscope::state::userCallback = myCallback;

  // Load scene
  scene = modules::read_scene(args::get(inputFilename));

  radius = scene.radius;
  h = scene.h;
  timestep = scene.timestep;
  rmax = scene.rmax == 0 ? radius * 10 : scene.rmax;

  if (scene.meshFileName == "") {
    std::cerr << "Please specify a mesh file in the scene file" << std::endl;
    return EXIT_FAILURE;
  }

  std::tie(mesh, geometry) = readManifoldSurfaceMesh(scene.meshFileName);

  auto psMesh = polyscope::registerSurfaceMesh("mesh", geometry->inputVertexPositions, mesh->getFaceVertexList());

  geometry->requireVertexTangentBasis();
  geometry->requireFaceTangentBasis();
  geometry->requireVertexNormals();
  geometry->requireVertexPrincipalCurvatureDirections();

  smoothedFunction = VertexData<double>(*mesh);

  if (scene.scalarFileName != "") {
    std::ifstream inFile;
    inFile.open(scene.scalarFileName);

    std::cout << "reading scalar file: " << scene.scalarFileName << std::endl;

    if (!inFile) {
      std::cerr << "Could not open file " << scene.scalarFileName << std::endl;
      exit(1);
    }

    int i = 0;
    for (std::string line; std::getline(inFile, line ); ) {
      if (line == "") {
        continue;
      }

      if (i >= mesh->nVertices()) {
        std::cerr << "error!: scalar file has more entries than the mesh" << std::endl;
        exit(1);
      }

      double d = std::stod(line);

      auto v = mesh->vertex(i);
      smoothedFunction[v] = d;

      // std::cout << i << ", " << d << std::endl;
      i++;
    }

    inFile.close();
  }
  else
  {
    for (Vertex v : mesh->vertices()) {
      auto coord = geometry->vertexPositions[v];
      // smoothedFunction[v] = std::pow(std::abs(coord.x), 1.5) + std::pow(std::abs(coord.y), 1.5);
      smoothedFunction[v] = 3*std::abs(geometry->vertexGaussianCurvature(v));
    }
  }

  {
    std::string scalarFilename = "scalar.txt";
    std::ofstream scalarfile(scalarFilename);

    for (int i = 0; i < smoothedFunction.size(); i++) {
      auto v = smoothedFunction[i];
      scalarfile << v << std::endl;
    }
  }
  

  VertexData<Vector3> vBasisX(*mesh);
  for (Vertex v : mesh->vertices()) {
    vBasisX[v] = geometry->vertexTangentBasis[v][0];
  }
  psMesh->setVertexTangentBasisX(vBasisX);
  
  if (scene.w_fieldAlignedness > 0) {
    vField =
      geometrycentral::surface::computeSmoothestVertexDirectionField(*geometry);
    polyscope::getSurfaceMesh("mesh")->addVertexIntrinsicVectorQuantity("vector field", vField);
  }

  polyscope::getSurfaceMesh("mesh")->addVertexIntrinsicVectorQuantity("principal curvature", geometry->vertexPrincipalCurvatureDirections);
  polyscope::getSurfaceMesh("mesh")->addVertexScalarQuantity("smooth function", smoothedFunction);

  {
    std::string vfieldFilename = "vfield.obj";
    std::ofstream vfieldfile(vfieldFilename);
    
    for (int i = 0; i < vField.size(); i++) {
      auto v = vField[i];
      auto x = geometry->vertexTangentBasis[i][0];
      auto y = geometry->vertexTangentBasis[i][1];

      auto v_cartesian = v.x * x + v.y * y;

      vfieldfile << "v " << v_cartesian.x << " " << v_cartesian.y << " " << v_cartesian.z << std::endl;
    }
  }

  if (scene.curveFileName != "") {
    std::tie(nodes, segments, isFixedNode) = modules::readCurveOnMesh(scene.curveFileName, *mesh, *geometry);
  }

  std::string datafilename = "data.csv";
  data_out = std::ofstream(datafilename);
  std::cout << data_out.is_open() << ", " << datafilename << std::endl;
  data_out << "iteration, time, numNodes, f, descent norm (L2), descent norm (L1), descent norm (L∞)" << std::endl;

  if (nodes.size() == 0) {
    for (auto v : mesh->face(0).adjacentVertices()) {
      std::cout << "initially added nodes: " << v << std::endl;
      auto sp = SurfacePoint(v);
      nodes.emplace_back(sp);
      isFixedNode.emplace_back(false);
    }

    for (int i = 0; i < nodes.size(); i++) {
      segments.emplace_back(std::array<int, 2>{i, int((i + 1) % nodes.size())});
    }
  }

  // add boundary loops to fixed nodes
  {
    std::vector<std::array<int, 2>> addedSegments = {};
    for (auto l : mesh->boundaryLoops()) {
      std::vector<int> addedNodes = {};
      std::cout << "adding loop: ";
      for (auto v : l.adjacentVertices())
      {
        std::cout << v << " ";
        nodes.emplace_back(SurfacePoint(v));
        isFixedNode.emplace_back(true);
        addedNodes.emplace_back(nodes.size() - 1);
      }
      std::cout << std::endl;

      for (int i = 0; i < addedNodes.size(); i++) {
        std::array<int, 2> newsegment = {addedNodes[i], addedNodes[(i + 1) % addedNodes.size()]};
        segments.emplace_back(newsegment);
        addedSegments.emplace_back(newsegment);
      }
    }

    auto ssp = std::vector<std::vector<SurfacePoint>>(addedSegments.size());
    auto [allNodes, allSegments, _] = formatGeodesicCurve(*mesh, *geometry, nodes, addedSegments, ssp);

    modules::write_curve("boundary_loops.obj", allNodes, allSegments);
  }

  std::tie(segmentSurfacePoints, segmentLengths) = modules::connect_surface_points(*mesh, *geometry, nodes, segments);

  std::tie(nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode) = modules::remesh_curve_on_surface(*mesh, *geometry, nodes, segments, segmentSurfacePoints, segmentLengths, isFixedNode, h);

  initialNodes = nodes;
  initialSegments = segments;
  initialSegmentSurfacePoints = segmentSurfacePoints;
  initialSegmentLengths = segmentLengths;
  isInitialFixedNode = isFixedNode;

  {
    // for (auto segment : initialSegments) {
    //   std::cout << "segment: " << segment[0] << ", " << segment[1] << std::endl;
    // }

    auto [allNodes, allSegments, cartesianCoords] = formatGeodesicCurve(*mesh, *geometry, initialNodes, initialSegments, initialSegmentSurfacePoints);

    polyscope::registerCurveNetwork("initial curve", allNodes, allSegments)
      // ->setColor({1.0, 1.0, 1.0})
      // ->setRadius(0.003)
    ;

    polyscope::registerCurveNetwork("curve", allNodes, allSegments);

    // polyscope::registerPointCloud("initial nodes", allNodes);

    std::vector<Vector3> fixedNodes = {}, freeNodes = {};
    for (int i = 0; i < initialNodes.size(); i++) {
      if (isInitialFixedNode[i]) {
        fixedNodes.emplace_back(cartesianCoords[i]);
      } else {
        freeNodes.emplace_back(cartesianCoords[i]);
      }
    }

    std::vector<Vector3> normalDirections, tangentX, tangentY;
    for (auto node : nodes) {
      auto [x, y] = modules::get_tangent_basis(*geometry, node);
      tangentX.emplace_back(x / 100);
      tangentY.emplace_back(y / 100);
      normalDirections.emplace_back(cross(x, y) / 100);
    }

    auto [nn, ns, _] = formatVector(cartesianCoords, normalDirections);
    polyscope::registerCurveNetwork("normal directions", nn, ns)->setEnabled(false);

    auto [txn, txs, __] = formatVector(cartesianCoords, tangentX);
    polyscope::registerCurveNetwork("tangent x", txn, txs)->setEnabled(false);

    auto [tyn, tys, ___] = formatVector(cartesianCoords, tangentY);
    polyscope::registerCurveNetwork("tangent y", tyn, tys)->setEnabled(false);

    polyscope::registerPointCloud("free nodes", freeNodes)
      // ->setPointColor({1.0, 1.0, 1.0})
    ;
    polyscope::registerPointCloud("fixed nodes", fixedNodes)
      // ->setPointColor({.5, .5, .5})
    ;

    modules::write_curve("./objs/curve_0.obj", allNodes, allSegments);

    {
      std::string fname = "./object.obj";
      std::ofstream objfile(fname);

      for (auto v : mesh->vertices()) {
        auto p = geometry->vertexPositions[v];
        objfile << "v " << p.x << " " << p.y << " " << p.z << std::endl;
      }
      objfile << std::endl;

      for (auto f : mesh->faces()) {
        auto vs = f.adjacentVertices();
        objfile << "f";

        for (auto v : vs) {
          objfile << " " << v.getIndex() + 1;
        }
        objfile << std::endl;
      }
    }

    // polyscope::registerPointCloud("free nodes", freeNodes)
    //   ->setPointColor({1.0, 1.0, 1.0})
    // ;
    // polyscope::registerPointCloud("fixed nodes", fixedNodes)
    //   ->setPointColor({.5, .5, .5})
    // ;
  }

  // Give control to the polyscope gui
  if (scene.excecuteOnly) {
    // for debug
    writeData = true;
    writeCurve = true;

    for (int i = 0; i < 100; i++) {
      doWork();
    }
  } else {
    polyscope::show();
  }

  return EXIT_SUCCESS;
}
