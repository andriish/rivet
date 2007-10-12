#include <iostream>
#include "Vectors.hh"

using namespace std;

int main() {
  Vector4X a(1,0,0,0);
  cout << a << ": interval = " << a.interval() << endl;
  a.z(1);
  cout << a << ": interval = " << a.interval() << endl;
  a.y(2).z(3);
  cout << a << ": interval = " << a.interval() << endl;
  cout << a << ": vector = " << a.vector3() << endl;
  return EXIT_SUCCESS;
}
