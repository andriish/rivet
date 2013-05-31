#ifndef RIVET_Cuts_HH
#define RIVET_Cuts_HH

#include <iostream>
#include <algorithm>
#include <boost/smart_ptr.hpp>
#include <Rivet/Particle.hh>
#include <Rivet/Jet.hh>
#include <Rivet/Math/Vectors.hh>
#include "fastjet/PseudoJet.hh"

using namespace std;

namespace Rivet {

/// Base class for all cut functors

class Cuttable {
    public:
    virtual double pT() const;
    virtual double m() const;
    virtual ~Cuttable() {}
};


struct CutBase {
    virtual bool cut(const Cuttable& o) const = 0;
    virtual ~CutBase() {}
};

typedef boost::shared_ptr<CutBase> CutPtr; // Cut;

template <typename T>
CutPtr make_cut(T t) {
    return boost::shared_ptr<T>(new T(t));
}

/// These pointers are used to access the cuts

CutPtr PtGtr(double n);

CutPtr PtLess(double n);

CutPtr PtIn(double n, double m);

CutPtr MassGtr(double n);

CutPtr MassLess(double n);

CutPtr MassIn(double n, double m);


/// AND, OR, NOT, and XOR objects for combining cuts

struct CutsOr : public CutBase {
    CutsOr(const CutPtr c1, const CutPtr c2);
    bool cut(const Cuttable& o) const;
    const CutPtr cut1;
    const CutPtr cut2;
};

struct CutsAnd : public CutBase {
    CutsAnd(const CutPtr c1, const CutPtr c2);
    bool cut(const Cuttable& o) const;
    const CutPtr cut1;
    const CutPtr cut2;
};

struct CutInvert : public CutBase {
    CutInvert(const CutPtr c1);
    bool cut(const Cuttable& o) const;
    const CutPtr poscut;
};

struct CutsXor : public CutBase {
    CutsXor(const CutPtr c1, const CutPtr c2);
    bool cut(const Cuttable& o) const;
    const CutPtr cut1;
    const CutPtr cut2;
};

/// operator &, operator |, operator ~, and operator ^ overloads

CutPtr operator & (const CutPtr aptr, const CutPtr bptr);

CutPtr operator | (const CutPtr aptr, const CutPtr bptr);

CutPtr operator ~ (const CutPtr cptr);

CutPtr operator ^ (const CutPtr aptr, const CutPtr bptr);

    /////////////////////////////////
    /// MakeCuttable - allows specialisation
    /// to deal with different access conventions

template <typename T>
class MakeCuttable : public Cuttable {};

template<>
class MakeCuttable <Particle> : public Cuttable {
    public:
    MakeCuttable(const Particle& p) : p_(p) {}
    double pT() const {return p_.momentum().pT();}
    double m() const {return p_.momentum().mass();}
    private:
    const Particle & p_;
};

template<>
class MakeCuttable <FourMomentum> : public Cuttable {
    public:
    MakeCuttable(const FourMomentum& fm) : fm_(fm) {}
    double pT() const {return fm_.pT();}
    double m() const {return fm_.mass();}
    private:
    const FourMomentum & fm_;
};

template<>
class MakeCuttable <Jet> : public Cuttable {
    public:
    MakeCuttable(const Jet& jet) : jet_(jet) {}
    double pT() const {return jet_.momentum().pT();}
    double m() const {return jet_.momentum().mass();}
    private:
    const Jet & jet_;
};

template<>
class MakeCuttable <fastjet::PseudoJet> : public Cuttable {
    public:
    MakeCuttable(const fastjet::PseudoJet& pjet) : pjet_(pjet) {}
    double pT() const {return pjet_.perp();}
    double m() const {return pjet_.m();}
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

}

#endif

