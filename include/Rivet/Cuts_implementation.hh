#ifndef RIVET_Cuts_implementation_HH
#define RIVET_Cuts_implementation_HH

#include <iostream>
#include <algorithm>
#include <boost/smart_ptr.hpp>
#include <Rivet/Particle.hh>
#include <Rivet/Jet.hh>
#include <Rivet/Math/Vectors.hh>
#include <fastjet/PseudoJet.hh>

using namespace std;

namespace Rivet {

    /// Base classes for cuts
class Cuttable {
    public:
    virtual double pT() const;
    virtual double m() const;
    virtual double y() const;
    virtual double eta() const;
    virtual double phi() const;
    virtual ~Cuttable() {}
};

struct CutBase {
    virtual bool cut(const Cuttable& o) const = 0;
    virtual ~CutBase() {}
};

typedef boost::shared_ptr<CutBase> Cut;

template <typename T>
Cut make_cut(T t) {
    return boost::shared_ptr<T>(new T(t));
}

/// AND, OR, NOT, and XOR objects for combining cuts

struct CutsOr : public CutBase {
    CutsOr(const Cut c1, const Cut c2);
    bool cut(const Cuttable& o) const;
    const Cut cut1;
    const Cut cut2;
};

struct CutsAnd : public CutBase {
    CutsAnd(const Cut c1, const Cut c2);
    bool cut(const Cuttable& o) const;
    const Cut cut1;
    const Cut cut2;
};

struct CutInvert : public CutBase {
    CutInvert(const Cut c1);
    bool cut(const Cuttable& o) const;
    const Cut poscut;
};

struct CutsXor : public CutBase {
    CutsXor(const Cut c1, const Cut c2);
    bool cut(const Cuttable& o) const;
    const Cut cut1;
    const Cut cut2;
};

/// operator &, operator |, operator ~, and operator ^ overloads

Cut operator & (const Cut aptr, const Cut bptr);

Cut operator | (const Cut aptr, const Cut bptr);

Cut operator ~ (const Cut cptr);

Cut operator ^ (const Cut aptr, const Cut bptr);


    ////////////////////////////////////////////
    /// MakeCuttable - uses specialisation
    /// to deal with different access conventions

template <typename T>
class MakeCuttable : public Cuttable {};

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

    /// pT cut structs
struct CutPtGtr : public CutBase {
    CutPtGtr(const double pt_lowerlim);
    bool cut(const Cuttable & o) const;
    private:
    double low_;
};

struct CutPtLess : public CutBase {
    CutPtLess(const double pt_upperlim);
    bool cut(const Cuttable & o) const;
    private:
    double high_;
};


    /// Mass cut structs
struct CutMassGtr : public CutBase {
    CutMassGtr(const double mass_lowerlim);
    bool cut (const Cuttable & o) const;
    private:
    double low_;
};

struct CutMassLess : public CutBase {
    CutMassLess(const double mass_upperlim);
    bool cut (const Cuttable& o) const;
    private:
    double high_;
};


    /// Rapidity cut structs
struct CutRapGtr : public CutBase {
    CutRapGtr(const double rap_lowerlim);
    bool cut (const Cuttable & o) const;
    private:
    double low_;
};

struct CutRapLess : public CutBase {
    CutRapLess(const double rap_upperlim);
    bool cut (const Cuttable& o) const;
    private:
    double high_;
};


    /// Pseudorapidity cut structs
struct CutEtaGtr : public CutBase {
    CutEtaGtr(const double eta_lowerlim);
    bool cut (const Cuttable & o) const;
    private:
    double low_;
};

struct CutEtaLess : public CutBase {
    CutEtaLess(const double eta_upperlim);
    bool cut (const Cuttable& o) const;
    private:
    double high_;
};


    /// Azimuthal angle cut structs
struct CutPhiGtr : public CutBase {
    CutPhiGtr(const double phi_lowerlim);
    bool cut (const Cuttable & o) const;
    private:
    double low_;
};

struct CutPhiLess : public CutBase {
    CutPhiLess(const double phi_upperlim);
    bool cut (const Cuttable& o) const;
    private:
    double high_;
};


}

#endif
