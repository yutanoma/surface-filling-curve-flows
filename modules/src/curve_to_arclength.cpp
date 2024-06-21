#include "modules/curve_to_arclength.h"

namespace modules {
  std::vector<int> find_starting_node(
    const std::vector<Vector3>& nodes,
    const std::vector<std::vector<int>> &node2segment
  ) {
    std::vector<int> startingNodes = {};
    for (int i = 0; i < nodes.size(); i++) {
      if (node2segment[i].size() != 2) {
        startingNodes.push_back(i);
      }
    }

    if (startingNodes.size() == 0) {
      startingNodes.push_back(0);
    }

    return startingNodes;
  }

  std::vector<double> curve_to_arclength(
    const std::vector<Vector3>& nodes,
    const std::vector<std::array<int, 2>> &segments
  ) {
    std::vector<std::vector<int>> node2segment(nodes.size());
    for (int i = 0; i < segments.size(); i++) {
      node2segment[segments[i][0]].push_back(i);
      node2segment[segments[i][1]].push_back(i);
    }

    std::vector<int> startingNodes = find_starting_node(nodes, node2segment);
    std::vector<bool> visited = std::vector<bool>(segments.size(), false);

    std::vector<double> arclengths = std::vector<double>(nodes.size(), 0.0);

    for (int i = 0; i < startingNodes.size(); i++) {
      int currentNode = startingNodes[i];
      int currentSegment = node2segment[currentNode][0];
      int nextNode = segments[currentSegment][0] == currentNode ? segments[currentSegment][1] : segments[currentSegment][0];
      visited[currentSegment] = true;

      arclengths[currentNode] = 0.0;
      arclengths[nextNode] = (nodes[nextNode] - nodes[currentNode]).norm();

      double l = (nodes[nextNode] - nodes[currentNode]).norm();

      currentNode = nextNode;

      while (true) {
        std::vector<int> nextSegments = node2segment[currentNode];
        int nextSegment = -1;
        for (int j = 0; j < nextSegments.size(); j++) {
          if (!visited[nextSegments[j]]) {
            nextSegment = nextSegments[j];
            break;
          }
        }

        if (nextSegment == -1) {
          break;
        }

        visited[nextSegment] = true;
        nextNode = segments[nextSegment][0] == currentNode ? segments[nextSegment][1] : segments[nextSegment][0];

        arclengths[nextNode] = l + (nodes[nextNode] - nodes[currentNode]).norm();
        l = arclengths[nextNode];

        currentNode = nextNode;
      }
    }

    return arclengths;
  }
}