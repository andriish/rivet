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

  namespace Cuts {
    enum Token {
      pt, mass, rap, eta, phi
    };
  }

  Cut operator < (Cuts::Token, double);
  inline Cut operator < (Cuts::Token tk, int i) { return tk < double(i); }
  Cut operator >= (Cuts::Token, double);
  inline Cut operator >= (Cuts::Token tk, int i) { return tk >= double(i); }


Cut ptIn(double n, double m);

Cut massIn(double n, double m);

Cut rapIn(double n, double m);

Cut etaIn(double n, double m);

Cut phiIn(double n, double m);




/// operator &, operator |, operator ~, and operator ^ overloads

Cut operator & (const Cut aptr, const Cut bptr);
Cut operator | (const Cut aptr, const Cut bptr);
Cut operator ~ (const Cut cptr);
Cut operator ^ (const Cut aptr, const Cut bptr);

}

#endif

