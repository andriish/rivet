#ifndef RIVET_Cuts_HH
#define RIVET_Cuts_HH
#include <boost/smart_ptr.hpp>

namespace Rivet {

  class CuttableBase;

  class CutBase {
  public:
    template <typename T>
    bool accept(const T & t);
  protected:
    virtual bool accept_(const CuttableBase & o) const = 0;
  public:
    virtual ~CutBase() {}
  };

typedef boost::shared_ptr<CutBase> Cut;

/// These functions are used to build cuts

  namespace Cuts {
    enum Quantity {
      pT, mass, rap, eta, phi
    };
    
    const Cut & open();
  }

  Cut operator < (Cuts::Quantity, double);
  Cut operator > (Cuts::Quantity, double);
  Cut operator <= (Cuts::Quantity, double);
  Cut operator >= (Cuts::Quantity, double);

  Cut Range(Cuts::Quantity, double m, double n);

  // overload helpers
  inline Cut operator < (Cuts::Quantity qty, int i) { 
    return qty < double(i); 
  }
  inline Cut operator > (Cuts::Quantity qty, int i) { 
    return qty > double(i); 
  }
  inline Cut operator <= (Cuts::Quantity qty, int i) { 
    return qty <= double(i); 
  }
  inline Cut operator >= (Cuts::Quantity qty, int i) { 
    return qty >= double(i); 
  }


/// operator &, operator |, operator ~, and operator ^ overloads

Cut operator & (const Cut aptr, const Cut bptr);
Cut operator | (const Cut aptr, const Cut bptr);
Cut operator ~ (const Cut cptr);
Cut operator ^ (const Cut aptr, const Cut bptr);

}

#endif

