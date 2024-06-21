#include "modules/curve_evolution.h"

namespace modules {
void local_area_backward_solve(
  std::vector<Vector2> &nodes,
  const std::vector<int> &localNodes,
  const std::vector<std::array<int, 2>> &fullSegments,
  const std::vector<std::array<int, 2>> &halfSegments,
  std::map<int, int> &activeId2Id,
  std::map<int, int> &id2ActiveId,
  std::map<int, double> &massPerVertex,
  const double timestep
) {
  int numActiveNodes = localNodes.size(), numFullSegments = fullSegments.size(), numHalfSegments = halfSegments.size();

  std::cout << "numActiveNodes: " << numActiveNodes << ", numFullSegments: " << numFullSegments << ", numHalfSegments: " << numHalfSegments << std::endl;

  auto func = TinyAD::scalar_function<2>(TinyAD::range(numActiveNodes));

  func.add_elements<2>(TinyAD::range(numFullSegments), [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int e_id = element.handle;

    Eigen::Vector2<T> v0 = element.variables(id2ActiveId[fullSegments[e_id][0]]);
    Eigen::Vector2<T> v1 = element.variables(id2ActiveId[fullSegments[e_id][1]]);

    Eigen::Vector2d _v0(nodes[fullSegments[e_id][0]].x, nodes[fullSegments[e_id][0]].y);
    Eigen::Vector2d _v1(nodes[fullSegments[e_id][1]].x, nodes[fullSegments[e_id][1]].y);

    Eigen::Vector2<T> grad = v1 - v0;

    double area = (_v1 - _v0).norm();

    return timestep * area * grad.squaredNorm();
  });

  func.add_elements<1>(TinyAD::range(numHalfSegments), [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int e_id = element.handle;

    Eigen::Vector2<T> v0 = element.variables(id2ActiveId[halfSegments[e_id][0]]);

    Eigen::Vector2d _v0(nodes[halfSegments[e_id][0]].x, nodes[halfSegments[e_id][0]].y);
    Eigen::Vector2d v1(nodes[halfSegments[e_id][1]].x, nodes[halfSegments[e_id][1]].y);

    double area = (v1 - _v0).norm();

    Eigen::Vector2<T> grad = v1 - v0;

    return timestep * area * grad.squaredNorm();
  });

  func.add_elements<1>(TinyAD::range(numActiveNodes), [&](auto &element) -> TINYAD_SCALAR_TYPE(element) {
    using T = TINYAD_SCALAR_TYPE(element);
    int v_idx = element.handle;

    Eigen::Vector2<T> v = element.variables(v_idx);
    Eigen::Vector2d v_prev(nodes[activeId2Id[v_idx]].x, nodes[activeId2Id[v_idx]].y);

    double mass = massPerVertex[activeId2Id[v_idx]];

    return (v - v_prev).squaredNorm() * mass;
  });

  auto x = func.x_from_data([&](int v_idx) {
    auto p = nodes[activeId2Id[v_idx]];

    return Eigen::Vector2d(p.x, p.y);
  });

  TinyAD::LinearSolver solver;
  auto [f, g, H_proj] = func.eval_with_hessian_proj(x);
  auto d = TinyAD::newton_direction(g, H_proj, solver);

  x = x + d;

  func.x_to_data(x, [&] (int v_idx, const Eigen::Vector2d& p) {
    nodes[activeId2Id[v_idx]] = Vector2{p.x(), p.y()};
  });
}

void local_area_timestep(
  std::vector<Vector2> &nodes,
  const std::vector<int> &localNodes,
  const double timeStep
) {
  std::vector<std::array<int, 2>> fullSegments = {}, halfSegments = {};
  std::map<int, bool> isInArea;
  std::map<int, int> id2LocalId, localId2id;
  std::map<int, double> massPerVertex;

  for (int i = 0; i < localNodes.size(); i++) {
    int id = localNodes[i];
    isInArea[id] = true;

    id2LocalId[id] = i;
    localId2id[i] = id;
  }

  for (int i = 0; i < localNodes.size(); i++) {
    int id = localNodes[i];
    int prev = (id + nodes.size() - 1) % nodes.size();
    int next = (id + 1) % nodes.size();

    if (isInArea[prev]) {
      // add full segment
      fullSegments.emplace_back(std::array<int, 2> {prev, id});
    } else {
      // add half segment
      halfSegments.emplace_back(std::array<int, 2> {id, prev});
    }

    if (!isInArea[next]) {
      // add half segment
      halfSegments.emplace_back(std::array<int, 2> {id, next});
    }

    double prevLen = (nodes[prev] - nodes[id]).norm();
    double nextLen = (nodes[next] - nodes[id]).norm();

    massPerVertex[id] = (prevLen + nextLen) / 2;
  }

  local_area_backward_solve(nodes, localNodes, fullSegments, halfSegments, localId2id, id2LocalId, massPerVertex, timeStep);
}

void global_area_timestep(
  std::vector<Vector2> &nodes,
  const std::vector<int> &globalNodes,
  const double timeStep
) {
  // todo
}

void curve_evolution(
  std::vector<Vector2> &nodes,
  const std::vector<int> &globalNodes,
  const std::vector<int> &localNodes,
  const double timeStep
) {
  if (localNodes.size() > 0) {
    local_area_timestep(nodes, localNodes, timeStep);
  }
  
  if (globalNodes.size() > 0) {
    global_area_timestep(nodes, globalNodes, timeStep);
  }
}
}
