#include <modules/curve_remesh_2d.h>

namespace modules {
std::vector<std::vector<Vector2>> curve_remesh_2d(
  const std::vector<std::vector<Vector2>> &curves,
  const double &h
) {
  std::vector<std::vector<Vector2>> newCurves = {};
  for (int i = 0; i < curves.size(); i++) {
    newCurves.emplace_back(std::vector<Vector2> {});

    int lastAdded = curves[i].size() - 1;
    for (int j = 0; j < curves[i].size(); j++) {
      // delete nodes that are too close
      if ((curves[i][j] - curves[i][lastAdded]).norm() < h) {
        // std::cout << "l16 deleted!" << std::endl;
        continue;
      }

      newCurves[i].emplace_back(curves[i][j]);
      lastAdded = j;
    }
  }

  auto refinedCurves = newCurves;
  newCurves = {};

  for (int i = 0; i < refinedCurves.size(); i++) {
    newCurves.emplace_back(std::vector<Vector2> {});

    int lastAdded = refinedCurves[i].size() - 1;
    for (int j = 0; j < refinedCurves[i].size(); j++) {
      // std::cout << "refinedCurves[" << i << "][" << j << "]: " << refinedCurves[i][j] << std::endl;

      // refinement
      Vector2 d = (refinedCurves[i][j] - refinedCurves[i][lastAdded]);
      if (d.norm() > 2 * h) {
        int divisionNum = std::ceil(d.norm() / (2 * h));

        for (int k = 1; k < divisionNum; k++) {
          // std::cout << "l41 added: [" << i << "][" << newCurves[i].size() << "]" << std::endl;
          newCurves[i].emplace_back(refinedCurves[i][lastAdded] + d * k / divisionNum);
        }
      }

      newCurves[newCurves.size() - 1].emplace_back(refinedCurves[i][j]);
      lastAdded = j;
    }
  }

  return newCurves;
}
}
