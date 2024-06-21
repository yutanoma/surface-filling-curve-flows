#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/surface_point.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace modules {
namespace SurfaceFillingEnergy {
struct Options {
  double p = 2;
  double q = 2;

  double w_fieldAlignedness = 0;
  VertexData<Vector2> vectorField;

  double w_curvatureAlignedness = 0;

  double w_bilaplacian = 0;

  // double w_alignedness = 0;
  // VertexData<double> alignedness;

  bool useAnisotropicAlphaOnMesh = false;
  VertexData<double> alphaRatioOnMesh;

  bool useAnisotropicAlphaOnNodes = false;
  std::vector<double> alphaRatioOnNodes;

  bool useGeodesicMedialAxis = false;
};
};

std::tuple<
  std::vector<Vector3> /* descent direction */,
  std::vector<Vector3> /* gradient direction */,
  double /* energy */,
  std::vector<std::vector<Vector3>> /* medial axis */
> surface_filling_energy_geodesic(
  ManifoldSurfaceMesh &mesh,
  VertexPositionGeometry &geometry,
  const std::vector<SurfacePoint> &nodes,
  const std::vector<std::array<int, 2>> &segments,
  const std::vector<std::vector<SurfacePoint>> &segmentSurfacePoints,
  const std::vector<double> &segmentLengths,
  const std::vector<bool> &isFixedNode,
  const std::vector<Vector3> &cartesianCoords,
  const double radius,
  const double maxRadius,
  const SurfaceFillingEnergy::Options &options = {}
);
}
