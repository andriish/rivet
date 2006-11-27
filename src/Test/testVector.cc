#include "Rivet/Tools/Math/LorentzVector.hh"
#include <iostream>

using namespace std;
using namespace Rivet;

int main() {
  LorentzVector myvec = LorentzVector::Z;
  cout << "Before: " << myvec << endl;

  LorentzVector myboostvec = 4.0 * LorentzVector(2.0, 0.0, 0.0, 1.0);
  cout << "Boost vector: " << myboostvec << endl;
  cout << "Boost vector mass: " << myboostvec.lorentzInvariant() << endl;

  //Poincare lorentzBoost(myboostvec);
  //lorentzBoost.Boost(myvec);
  //cout << "After: " << myvec << endl;

  return EXIT_SUCCESS;
}
