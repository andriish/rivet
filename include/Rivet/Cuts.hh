#ifndef RIVET_Cuts_HH
#define RIVET_Cuts_HH

#include <Rivet/Cuts_implementation.hh>

using namespace std;

namespace Rivet {

/// These pointers are used to build cuts

Cut ptGtr(double n);
Cut ptLess(double n);
Cut ptIn(double n, double m);

Cut massGtr(double n);
Cut massLess(double n);
Cut massIn(double n, double m);

Cut rapGtr(double n);
Cut rapLess(double n);
Cut rapIn(double n, double m);

Cut etaGtr(double n);
Cut etaLess(double n);
Cut etaIn(double n, double m);

Cut phiGtr(double n);
Cut phiLess(double n);
Cut phiIn(double n, double m);

}

#endif

