#include <Rivet/Particle.hh>
#include <Rivet/Jet.hh>
#include <Rivet/Math/Vectors.hh>
#include <fastjet/PseudoJet.hh>
#include <Rivet/Cuts.hh>

namespace Rivet {

class Cuttable {
public:
    virtual double pT() const;
    virtual double m() const;
    virtual double y() const;
    virtual double eta() const;
    virtual double phi() const;
    virtual ~Cuttable() {}
};


template <>
bool CutBase::accept<Cuttable>(const Cuttable & t) {
  return accept_(t);
}



#define CUTOPS(CUTNAME,CUTFN) \
class Cut_ ## CUTNAME ## Gtr : public CutBase { \
public: \
    Cut_ ## CUTNAME ## Gtr(const double low) : low_(low) {} \
protected: \
    bool accept_(const Cuttable & o) const { return (CUTFN) >= low_; } \
    private: \
    double low_; \
}; \
\
class Cut_ ## CUTNAME ## Less : public CutBase { \
public: \
    Cut_ ## CUTNAME ## Less(const double high) : high_(high) {} \
protected: \
    bool accept_(const Cuttable & o) const { return (CUTFN) < high_; } \
    private: \
    double high_; \
}; \
\
 Cut CUTNAME ## In(double m, double n) { \
   if (m > n) swap(m,n); \
   return (Cuts::CUTNAME >= m) & (Cuts::CUTNAME < n);	\
} \




  template <typename T>
  Cut make_cut(T t) {
    return boost::shared_ptr<T>(new T(t));
  }



  CUTOPS(pt, o.pT())

  CUTOPS(mass, o.m())

  CUTOPS(rap, o.y())

  CUTOPS(eta, o.eta())

  CUTOPS(phi, o.phi())


  Cut operator < (Cuts::Token tk, double n) {
    switch (tk) {
    case Cuts::pt:  return make_cut(Cut_ptLess(n));
    case Cuts::mass:return make_cut(Cut_massLess(n));
    case Cuts::rap: return make_cut(Cut_rapLess(n));
    case Cuts::eta: return make_cut(Cut_etaLess(n));
    case Cuts::phi: return make_cut(Cut_phiLess(n));
    }
  }

  Cut operator >= (Cuts::Token tk, double n) {
    switch (tk) {
    case Cuts::pt:  return make_cut(Cut_ptGtr(n));
    case Cuts::mass:return make_cut(Cut_massGtr(n));
    case Cuts::rap: return make_cut(Cut_rapGtr(n));
    case Cuts::eta: return make_cut(Cut_etaGtr(n));
    case Cuts::phi: return make_cut(Cut_phiGtr(n));
    }
  }


    //////////////
    /// Combiners

/// AND, OR, NOT, and XOR objects for combining cuts

class CutsOr : public CutBase {
public:
    CutsOr(const Cut c1, const Cut c2) : cut1(c1), cut2(c2) {}
protected:
    bool accept_(const Cuttable & o) const {
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
    bool accept_(const Cuttable & o) const {
      return cut1->accept(o) && cut2->accept(o);
    }
private:
    const Cut cut1;
    const Cut cut2;
};

class CutInvert : public CutBase {
public:
    CutInvert(const Cut c1) : poscut(c1) {}
protected:
    bool accept_(const Cuttable & o) const {
      return !poscut->accept(o);
    }
private:
    const Cut poscut;
};

class CutsXor : public CutBase {
public:
    CutsXor(const Cut c1, const Cut c2) : cut1(c1), cut2(c2) {}
protected:
    bool accept_(const Cuttable & o) const {
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




  double Cuttable::pT() const { assert(false); }
  double Cuttable::m() const { assert(false); }
  double Cuttable::y() const { assert(false); }
  double Cuttable::eta() const { assert(false); }
  double Cuttable::phi() const { assert(false); }



template <typename T>
class MakeCuttable : public Cuttable {};


#define SPECIALISE_ACCEPT(TYPENAME) \
template <> \
bool CutBase::accept<TYPENAME>(const TYPENAME & t) { \
  return accept_(MakeCuttable<TYPENAME>(t)); \
} \



template<>
class MakeCuttable <Particle> : public Cuttable {
    public:
    MakeCuttable(const Particle& p) : p_(p) {}
    double pT() const {return p_.momentum().pT();}
    double m() const {return p_.momentum().mass();}
    double y() const {return p_.momentum().rapidity();}
    double eta() const {return p_.momentum().pseudorapidity();}
    double phi() const {return p_.momentum().phi();}
    private:
    const Particle & p_;
};
SPECIALISE_ACCEPT(Particle)


template<>
class MakeCuttable <FourMomentum> : public Cuttable {
    public:
    MakeCuttable(const FourMomentum& fm) : fm_(fm) {}
    double pT() const {return fm_.pT();}
    double m() const {return fm_.mass();}
    double y() const {return fm_.rapidity();}
    double eta() const {return fm_.pseudorapidity();}
    double phi() const {return fm_.phi();}
    private:
    const FourMomentum & fm_;
};
SPECIALISE_ACCEPT(FourMomentum)

template<>
class MakeCuttable <Jet> : public Cuttable {
    public:
    MakeCuttable(const Jet& jet) : jet_(jet) {}
    double pT() const {return jet_.momentum().pT();}
    double m() const {return jet_.momentum().mass();}
    double y() const {return jet_.momentum().rapidity();}
    double eta() const {return jet_.momentum().pseudorapidity();}
    double phi() const {return jet_.momentum().phi();}
    private:
    const Jet & jet_;
};
SPECIALISE_ACCEPT(Jet)

template<>
class MakeCuttable <fastjet::PseudoJet> : public Cuttable {
    public:
    MakeCuttable(const fastjet::PseudoJet& pjet) : pjet_(pjet) {}
    double pT() const {return pjet_.perp();}
    double m() const {return pjet_.m();}
    double y() const {return pjet_.rap();}
    double eta() const {return pjet_.eta();}
    double phi() const {return pjet_.phi();}
    private:
    const fastjet::PseudoJet & pjet_;
};
SPECIALISE_ACCEPT(fastjet::PseudoJet)

}
