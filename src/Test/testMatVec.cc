#include <iostream>
#include "Vectors.hh"
#include "Matrices.hh"

using namespace std;

int main() {
  Vector4X a(1,0,0,0);
  cout << a << ": interval = " << a.interval() << endl;
  a.z(1);
  cout << a << ": interval = " << a.interval() << endl;
  a.y(2).z(3);
  cout << a << ": interval = " << a.interval() << endl;
  cout << a << ": vector = " << a.vector3() << endl;


  Matrix<3> m;
  m.set(0, 0,  7/4.0);
  m.set(0, 1, 3*sqrt(3)/4.0);
  m.set(1, 0, 3*sqrt(3)/4.0);
  m.set(1, 1, 13/4.0);
  m.set(2, 2,  9);

  Eigensystem<3> es;
  es.solve(m).show();

  return EXIT_SUCCESS;
}
