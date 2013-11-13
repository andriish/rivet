#include <Rivet/Particle.hh>
#include <Rivet/Jet.hh>
#include <Rivet/Math/Vectors.hh>
#include <fastjet/PseudoJet.hh>
#include <Rivet/Cuts.hh>
#include <boost/make_shared.hpp>

namespace Rivet {

class CuttableBase {
public:
  virtual double getValue(Cuts::Quantity) const = 0;
  virtual ~CuttableBase() {}
};


template <>
bool CutBase::accept<CuttableBase>(const CuttableBase & t) {
  return accept_(t);
}


class Open_Cut : public CutBase {
protected:
  bool accept_(const CuttableBase &) const { return true; }
};

const Cut & Cuts::open() {
  static const Cut open = boost::make_shared<Open_Cut>();
  return open;
}

class Cut_Gtr : public CutBase {
public:
  Cut_Gtr(const Cuts::Quantity qty, const double low) : qty_(qty), low_(low) {}
protected:
  bool accept_(const CuttableBase & o) const { return o.getValue(qty_) >= low_; }
private:
  Cuts::Quantity qty_;
  double low_;
};

class Cut_Less : public CutBase {
public:
  Cut_Less(const Cuts::Quantity qty, const double high) : qty_(qty), high_(high) {}
protected:
  bool accept_(const CuttableBase & o) const { return o.getValue(qty_) < high_; }
private:
  Cuts::Quantity qty_;
  double high_;
};


  template <typename T>
  inline Cut make_cut(T t) {
    return boost::make_shared<T>(t);
  }

  Cut operator < (Cuts::Quantity qty, double n) {
    return make_cut(Cut_Less(qty, n));
  }

  Cut operator >= (Cuts::Quantity qty, double n) {
    return make_cut(Cut_Gtr(qty, n));
  }

  Cut Range(Cuts::Quantity qty, double m, double n) {
    if (m > n) swap(m,n);
    return (qty >= m) & (qty < n);
  }


    //////////////
    /// Combiners

/// AND, OR, NOT, and XOR objects for combining cuts

class CutsOr : public CutBase {
public:
    CutsOr(const Cut c1, const Cut c2) : cut1(c1), cut2(c2) {}
protected:
    bool accept_(const CuttableBase & o) const {
      return cut1->accept(o) || cut2->accept(o);
    }
private:
    const Cut cut1;
    const Cut cut2;
};

class CutsAnd : public CutBase {
public:
    CutsAnd(const Cut c1, const Cut c2) : cut1(c1), cut2(c2) {}
protected:
    bool accept_(const CuttableBase & o) const {
      return cut1->accept(o) && cut2->accept(o);
    }
private:
    const Cut cut1;
    const Cut cut2;
};

class CutInvert : public CutBase {
public:
    CutInvert(const Cut c1) : cut(c1) {}
protected:
    bool accept_(const CuttableBase & o) const {
      return !cut->accept(o);
    }
private:
    const Cut cut;
};

class CutsXor : public CutBase {
public:
    CutsXor(const Cut c1, const Cut c2) : cut1(c1), cut2(c2) {}
protected:
    bool accept_(const CuttableBase & o) const {
      bool A_and_B = cut1->accept(o) && cut2->accept(o);
      bool A_or_B  = cut1->accept(o) || cut2->accept(o);
      return A_or_B && (! A_and_B);
    }
private:
    const Cut cut1;
    const Cut cut2;
};

    ////////////
    ///Operators

  Cut operator & (const Cut aptr, const Cut bptr) {
    return make_cut(CutsAnd(aptr,bptr));
  }

  Cut operator | (const Cut aptr, const Cut bptr) {
    return make_cut(CutsOr(aptr,bptr));
  }

  Cut operator ~ (const Cut cptr) {
    return make_cut(CutInvert(cptr));
  }

  Cut operator ^ (const Cut aptr, const Cut bptr) {
    return make_cut(CutsXor(aptr,bptr));
  }

  ///////////////////////
  /// Cuts


template <typename T>
class Cuttable : public CuttableBase {};


#define SPECIALISE_ACCEPT(TYPENAME) \
template <> \
bool CutBase::accept<TYPENAME>(const TYPENAME & t) { \
  return accept_(Cuttable<TYPENAME>(t)); \
} \


void qty_not_found() {
  throw Exception("Missing implementation for a Quantity.");
}


template<>
class Cuttable <Particle> : public CuttableBase {
    public:
    Cuttable(const Particle& p) : p_(p) {}
    double getValue(Cuts::Quantity qty) const {
      switch(qty) {
      case Cuts::pT:   return p_.momentum().pT();
      case Cuts::mass: return p_.momentum().mass();
      case Cuts::rap:  return p_.momentum().rapidity();
      case Cuts::eta:  return p_.momentum().pseudorapidity();
      case Cuts::phi:  return p_.momentum().phi();
      default: 
	qty_not_found();
      }
      return -999.;
    }
    private:
    const Particle & p_;
};
SPECIALISE_ACCEPT(Particle)


template<>
class Cuttable <FourMomentum> : public CuttableBase {
    public:
    Cuttable(const FourMomentum& fm) : fm_(fm) {}
    double getValue(Cuts::Quantity qty) const {
      switch(qty) {
      case Cuts::pT:   return fm_.pT();
      case Cuts::mass: return fm_.mass();
      case Cuts::rap:  return fm_.rapidity();
      case Cuts::eta:  return fm_.pseudorapidity();
      case Cuts::phi:  return fm_.phi();
      default: 
	qty_not_found();
      }
      return -999.;
    }
    private:
    const FourMomentum & fm_;
};
SPECIALISE_ACCEPT(FourMomentum)

template<>
class Cuttable <Jet> : public CuttableBase {
    public:
    Cuttable(const Jet& jet) : jet_(jet) {}
    double getValue(Cuts::Quantity qty) const {
      switch(qty) {
      case Cuts::pT:   return jet_.momentum().pT();
      case Cuts::mass: return jet_.momentum().mass();
      case Cuts::rap:  return jet_.momentum().rapidity();
      case Cuts::eta:  return jet_.momentum().pseudorapidity();
      case Cuts::phi:  return jet_.momentum().phi();
      default: 
	qty_not_found();
      }
      return -999.;
    }
    private:
    const Jet & jet_;
};
SPECIALISE_ACCEPT(Jet)

template<>
class Cuttable <fastjet::PseudoJet> : public CuttableBase {
    public:
    Cuttable(const fastjet::PseudoJet& pjet) : pjet_(pjet) {}
    double getValue(Cuts::Quantity qty) const {
      switch(qty) {
      case Cuts::pT:   return pjet_.perp();
      case Cuts::mass: return pjet_.m();
      case Cuts::rap:  return pjet_.rap();
      case Cuts::eta:  return pjet_.eta();
      case Cuts::phi:  return pjet_.phi();
      default: 
	qty_not_found();
      }
      return -999.;
    }
    private:
    const fastjet::PseudoJet & pjet_;
};
SPECIALISE_ACCEPT(fastjet::PseudoJet)

}
