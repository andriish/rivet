#include <iostream>
#include <limits>

#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Math/Matrices.hh"

using namespace std;
using namespace Rivet;

int main() {

  FourVector a(1,0,0,0);
  cout << a << ": interval = " << a.invariant() << endl;
  a.setZ(1);
  cout << a << ": interval = " << a.invariant() << endl;
  a.setY(2).setZ(3);
  cout << a << ": interval = " << a.invariant() << endl;
  cout << a << ": vector = " << a.vector3() << endl << endl;

  FourMomentum b(1,0,0,1);
  cout << b << ": mass = " << b.mass() << endl;
  b.setPz(1);
  cout << b << ": mass = " << b.mass() << endl;
  b.setPy(2).setPz(3).setE(6);
  cout << b << ": mass = " << b.mass2() << endl;
  cout << b << ": vector = " << b.vector3() << endl << endl;

  Matrix3 m;
  m.set(0, 0, 7/4.0);
  m.set(0, 1, 3 * sqrt(3)/4.0);
  m.set(1, 0, 3 * sqrt(3)/4.0);
  m.set(1, 1, 13/4.0);
  m.set(2, 2, 9);
  cout << m << endl << endl;
  EigenSystem<3> es = diagonalize(m);
  /// @todo Fix the EigenSystem operator<< and toString() function
  //cout << "Eigensolns = " << endl << toString(es) << endl << endl;

  cout << "Matrices:" << endl;
  cout << Matrix3() << endl;
  cout << Matrix3::mkIdentity() << endl;
  const Matrix3 I3 = Matrix3::mkIdentity();
  cout << Matrix3::mkIdentity() * m * I3 << endl;
  cout << "tr(0) & det(0): " << Matrix3().trace() << ", " << Matrix3().det() << endl;
  cout << "tr(I3) & det(I3): " << I3.trace() << ", " << I3.det() << endl;
  Matrix3 m1 = Matrix3::mkIdentity();
  Matrix3 m2 = m1;
  m1.setRow(1, Vector3(1,2,3));
  m2.setColumn(1, Vector3(3,2,1));
  Matrix3 m3 = Matrix3::mkZero();
  cout << m1 << " + " << m2 << " = " << m1 + m2 << endl;
  m3.setRow(0, Vector3(2,3,0)).setRow(1, Vector3(1,4,3)).setRow(2, Vector3(0,1,2));
  cout << m1+m2 << " == " << m3 << ": " << (m1+m2 == m3 ? "true" : "false") << endl;
  cout << endl;
  
  Vector3 v3(1,2,3);
  cout << "Vector: " << v3 << endl;
  cout << "Invert: " << v3 << " --> " << -v3 << endl;
  const Matrix3 rot90(Vector3(0,0,1), PI/2.0);
  cout << "Rot 90: " << v3 << " --90deg--> " << rot90*v3 << endl;
  cout << "Rot 2 x 90: " << v3 << " --90deg--> " << rot90*rot90*v3 << endl;
  cout << "Rot 90*-90: "<< v3 << " --90deg--> " << rot90*rot90.inverse()*v3 << endl;
  const Matrix3 rot1(Vector3(0,1,0), PI/180.0);
  cout << "Rot 0 x 45 x 1: " << v3 << endl;
  for (size_t i = 0; i < 8; ++i) {
    for (size_t j = 0; j < 45; ++j) {
      v3 = rot1*v3;
    }
    cout << "Rot " << i+1 << " x 45 x 1: " << v3 << endl;
  }
  assert(v3 == Vector3(1,2,3));
  cout << endl;

  cout << "Boosts:" << endl;
  LorentzTransform ltX(0.5,0,0);
  cout << "LTx: " << ltX << endl;
  cout << "I on LTx: " << ltX.rotate(Matrix3::mkIdentity()) << endl;
  cout << "Rot90 on LTx: " << ltX.rotate(rot90) << endl;
  cout << endl;  

  cout << "X-boosts:" << endl;
  const FourMomentum p1 = FourMomentum(10,0,0,1);
  const FourMomentum p2 = ltX.transform(p1);
  cout << p1 << " -> " << p2 << endl;
  cout << p2 << " -> " << ltX.inverse().transform(p2) << endl;
  //cout << p1.boostVector() << endl;
  const FourMomentum p3 = LorentzTransform(-p1.boostVector()).transform(p1);
  cout << p1 << " -> " << p3 << endl;
  cout << endl;

  LorentzTransform ltY(0,0.4,0);
  cout << FourMomentum(1,0,0,1) << " -> " //<< "\n  " 
       << (ltX * ltY).transform(FourMomentum(1,0,0,1)) << endl;
  cout << FourMomentum(1,0,0,1) << " -> " //<< "\n  " 
       << (ltY * ltX).transform(FourMomentum(1,0,0,1)) << endl;
  cout << (ltX * ltY).boost() << endl;
  cout << (ltY * ltX).boost() << endl;
  cout << (ltX * ltX.inverse()).boost() << endl;

  return EXIT_SUCCESS;
}
