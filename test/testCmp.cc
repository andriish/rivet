#include <iostream>
#include <limits>

#include "Rivet/Tools/Cmp.hh"

using namespace std;

ostream & operator<<(ostream & os, Rivet::CmpState c) {
  string s;
  switch (c) {
    case Rivet::CmpState::UNDEF : s = "UNDEF"; break;
    case Rivet::CmpState::LT : s = "LT"; break;
    case Rivet::CmpState::EQ : s = "EQ"; break;
    case Rivet::CmpState::GT : s = "GT"; break;
  }
  os << s;
  return os;
}

int main() {
  using namespace Rivet;

  CmpState cs = CmpState::UNDEF;

  cs = cmp(0.5, 0.6);
  cout << "cmp(0.5, 0.6) = " << cs << '\n';
  assert(cs == CmpState::LT);

  cs = cmp(0.5, 0.5);
  cout << "cmp(0.5, 0.5) = " << cs << '\n';
  assert(cs == CmpState::EQ);

  cs = cmp(0.6, 0.5);
  cout << "cmp(0.6, 0.5) = " << cs << '\n';
  assert(cs == CmpState::GT);

  cs = cmp(1.,1.) || cmp(0.6, 0.5);
  cout << "cmp(1.,1.) || cmp(0.6, 0.5) = " << cs << '\n';
  assert(cs == CmpState::GT);

  return EXIT_SUCCESS;
}
