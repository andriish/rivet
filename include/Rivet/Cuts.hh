#ifndef RIVET_Cuts_HH
#define RIVET_Cuts_HH
#include <boost/smart_ptr.hpp>

namespace Rivet {

  class Cuttable;

  class CutBase {
  public:
    template <typename T>
    bool accept(const T & t);
  protected:
    virtual bool accept_(const Cuttable & o) const = 0;
  public:
    virtual ~CutBase() {}
  };

typedef boost::shared_ptr<CutBase> Cut;

/// These functions are used to build cuts


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




/// operator &, operator |, operator ~, and operator ^ overloads

Cut operator & (const Cut aptr, const Cut bptr);
Cut operator | (const Cut aptr, const Cut bptr);
Cut operator ~ (const Cut cptr);
Cut operator ^ (const Cut aptr, const Cut bptr);

}

#endif

