#include "modules/write_curve.h"

#include <iostream>

namespace modules {
void write_curve(
  const std::string &fname,
  std::vector<Vector3> &nodes,
  std::vector<std::array<int, 2>> &segments
) {
  std::ofstream myfile(fname);

  // write vertices
  for(const Vector3& v : nodes) {
      myfile << "v " << v.x << " " << v.y << " " << v.z << std::endl;
  }

  // write polylines
  for(const auto& c : segments) {
      myfile << "l";
      for( size_t i : c ) {
          myfile << " " << (1+i);
      }
      myfile << std::endl;
  }
  myfile.close();
}
}
