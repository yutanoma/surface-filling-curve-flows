#include "modules/global_area.h"

#include <queue>
#include <map>

#include <iostream>

namespace modules {
std::vector<int> extract_global_area_3d(
  const Eigen::MatrixXi &T,
  const Eigen::MatrixXi &T_FE,
  const std::vector<std::vector<int>> &T_EF,
  const std::vector<std::vector<int>> &T_ET,
  const std::vector<int> &innerEdges,
  const std::vector<int> &innerFaces,
  const std::vector<double> &radiusPerTet,
  const double _radius
) {
  std::map<int, bool> isInsideEdge;
  for (auto eid : innerEdges) {
    isInsideEdge[eid] = true;
  }

  std::map<int, bool> hasSmallRadius;
  for (int i = 0; i < T.rows(); i++) {
    if (radiusPerTet[i] < _radius) {
      hasSmallRadius[i] = true;
    }
  }

  // edges correspond to cells
  // faces correspond to segments
  // tets correspond to nodes

  // mark cells to put inside the graph
  std::map<int, bool> isInsideGraph;
  for (int i = 0; i < innerEdges.size(); i++) {
    int eid = innerEdges[i];

    auto tets = T_ET[eid];

    // this cell is in the graph if at least one node on the cell has small radius
    bool insideGraph = false;
    for (auto tid : tets) {
      if (hasSmallRadius[tid]) {
        insideGraph = true;
        break;
      }
    }

    isInsideGraph[eid] = insideGraph;
  }

  // build a graph of segments and nodes
  // cell to innerEdge is a 1-to-1 correpondence
  std::vector<std::vector<int>> segments2Cell(innerFaces.size(), std::vector<int> {});
  std::vector<std::vector<int>> cell2Segment(innerEdges.size(), std::vector<int> {});
  std::queue<int> q;

  std::cout << "inneredgesnum: " << innerEdges.size() << std::endl;

  std::vector<int> cell2Edge = innerEdges;
  std::map<int, int> edge2Cell;
  for (int i = 0; i < innerEdges.size(); i++) {
    edge2Cell[innerEdges[i]] = i;
  }

  std::vector<int> face2Segment = innerFaces;
  std::map<int, int> segment2Face;
  for (int i = 0; i < innerFaces.size(); i++) {
    segment2Face[innerFaces[i]] = i;
  }

  for (int i = 0; i < innerFaces.size(); i++) {
    int fid = innerFaces[i];
    int segmentId = i;

    int incidentEdgesNum = 0;

    for (int j = 0; j < 3; j++) {
      int eid = T_FE(fid, j);

      assert(eid != -1);

      if (isInsideEdge[eid]) {
        incidentEdgesNum++;
      }

      if (!isInsideGraph[eid]) {
        continue;
      }

      if (isInsideEdge[eid]) {
        int cellId = edge2Cell[eid];

        segments2Cell[segmentId].emplace_back(cellId);
        cell2Segment[cellId].emplace_back(segmentId);
      }
    }

    // add cell on the boundary
    if (incidentEdgesNum == 1 && segments2Cell[segmentId].size() == 1) {
      q.push(segments2Cell[segmentId][0]);
    }
  }

  // if not all the nodes on this cell have small radius and they can divided into several chunks, then you will have to divide the cell
  for (int i = 0; i < innerEdges.size(); i++) {
    if (!isInsideGraph[innerEdges[i]]) {
      continue;
    }

    auto tets = T_ET[innerEdges[i]];
    auto faces = T_EF[innerEdges[i]];

    // if this cell has updowns, then you will have to divide the cell
    std::vector<std::vector<int>> segmentGroups = {};

    bool currentMode = false;

    for (int j = 0; j < tets.size(); j++) {
      int t = tets[j];
      double _r = radiusPerTet[t];

      if (!currentMode) {
        if (_r < _radius) {
          // current mode is false and _r < _radius
          // switch to true
          currentMode = true;

          // create a new group
          segmentGroups.emplace_back(std::vector<int> {j});
        } else {
          // current mode is false and _r >= _radius
          // do nothing
        }
      } else {
        // current mode is true
        if (_r < _radius) {
          // current mode is true and _r < _radius
          // add to the end of the current group
          segmentGroups.back().emplace_back(j);
        } else {
          // current mode is true and _r >= _radius
          // switch to false
          currentMode = false;
        }
      }
    }

    assert(segmentGroups.size() != 0);

    if (segmentGroups.size() > 1) {
      // duplicate cells
      std::vector<int> group2Cells(segmentGroups.size());
      for (int j = 1; j < segmentGroups.size(); j++) {
        cell2Segment.emplace_back(std::vector<int> {});
        cell2Edge.emplace_back(innerEdges[i]);
        group2Cells[j] = cell2Edge.size() - 1;
      }

      std::map<int, int> node2group;
      std::map<int, bool> node2groupSet;
      for (int j = 0; j < segmentGroups.size(); j++) {
        for (auto _j : segmentGroups[j]) {
          node2group[tets[_j]] = j;
          node2groupSet[tets[_j]] = true;
        }
      }

      for (auto f : faces) {
        int t0 = T_FE(f, 0), t1 = T_FE(f, 1);

        assert(t0 != -1);
        assert(t1 != -1);

        int group;
        if (node2groupSet[t0]) {
          group = node2group[t0];
        } else if (node2groupSet[t1]) {
          group = node2group[t1];
        } else {
          continue;
        }

        if (group == 0) {
          // if group is zero then do not change
          continue;
        }

        // swap the cell to the new group
        for (int j = 0; j < segments2Cell[f].size(); j++) {
          auto c = segments2Cell[f][j];
          if (c == innerEdges[i]) {
            segments2Cell[f][j] = group2Cells[group];
            cell2Segment[group2Cells[group]].emplace_back(f);
            break;
          }
        }
      }
    }
  }

  {
    // for debug
    // q = std::queue<int> {};
    std::cout << "queue size: " << q.size() << ", isEmpty?: " << q.empty() << std::endl;
  }

  while (!q.empty()) {
    int eid = q.front();
    q.pop();

    for (auto fid : cell2Segment[eid]) {
      std::vector<int> newList = {};
      for (auto _eid : segments2Cell[fid]) {
        if (_eid != eid) {
          newList.emplace_back(_eid);
        }
      }

      assert(newList.size() == segments2Cell[fid].size() - 1);

      if (newList.size() == 1) {
        q.push(newList[0]);
      }

      segments2Cell[fid] = newList;
    }

    cell2Segment[eid] = {};
  }

  // mark the remaining cells as global areas
  std::vector<int> globalArea = {};
  std::map<int, bool> isGlobalArea;
  for (int i = 0; i < cell2Edge.size(); i++) {
    int eid = cell2Edge[i];

    std::cout << "eid: " << eid << ", size: " << cell2Segment[i].size() << ", " << isGlobalArea[eid] << std::endl;

    if (cell2Segment[i].size() > 0 && !isGlobalArea[eid]) {
      globalArea.emplace_back(eid);
      isGlobalArea[eid] = true;
    }
  }

  std::cout << globalArea.size() << std::endl;

  return globalArea;
}
}
