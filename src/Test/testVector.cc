#include "Rivet/Tools/Math/Vector.h"
#include "Rivet/Tools/Math/Poincare.h"
#include <iostream>

using namespace std;
using namespace Rivet;

int main() {
  Vec4D myvec = Vec4D::ZVEC;
  cout << "Before: " << myvec << endl;

  Vec4D myboostvec = 4.0 * Vec4D(0.0, 0.0, 0.0, 1.0);
  cout << "Boost vector: " << myboostvec << endl;
  cout << "Boost vector mass: " << myboostvec.Mass() << endl;

  Poincare lorentzBoost(myboostvec);
  lorentzBoost.Boost(myvec);
  cout << "After: " << myvec << endl;

  return EXIT_SUCCESS;
}
