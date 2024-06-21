#include <modules/circumcenter.h>

namespace modules {
Vector2 circumcenter(Vector2 A, Vector2 B, Vector2 C) {
  double D = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));

  if (std::abs(D) < 1e-9) {
    std::cout << "A: " << A << ", B: " << B << ", C: " << C << std::endl;
    return Vector2 {INFINITY, INFINITY};
      throw std::runtime_error("Points are collinear. Circumcenter doesn't exist.");
  }

  double Ux = ((A.x * A.x + A.y * A.y) * (B.y - C.y) + (B.x * B.x + B.y * B.y) * (C.y - A.y) + (C.x * C.x + C.y * C.y) * (A.y - B.y)) / D;
  double Uy = ((A.x * A.x + A.y * A.y) * (C.x - B.x) + (B.x * B.x + B.y * B.y) * (A.x - C.x) + (C.x * C.x + C.y * C.y) * (B.x - A.x)) / D;

  return Vector2 {Ux, Uy};
}
}
