#include <iostream>
#include <limits>

#include "Rivet/Tools/Cmp.hh"

using namespace std;

int main() {
  using namespace Rivet;

  CmpState cs = UNDEFINED;

  cs = cmp(0.5, 0.6);
  cout << "cmp(0.5, 0.6) = " << cs << '\n';
  assert(cs == ORDERED);

  cs = cmp(0.5, 0.5);
  cout << "cmp(0.5, 0.5) = " << cs << '\n';
  assert(cs == EQUIVALENT);

  cs = cmp(0.6, 0.5);
  cout << "cmp(0.6, 0.5) = " << cs << '\n';
  assert(cs == UNORDERED);

  return EXIT_SUCCESS;
}
