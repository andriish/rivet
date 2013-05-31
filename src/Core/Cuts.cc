#include "Rivet/Cuts.hh"

using namespace std;
namespace Rivet {

    ///////////////////
    /// Access pointers

  Cut ptGtr(double n) {
    CutPtGtr pg(n);
    return make_cut(pg);
  }

  Cut ptLess(double n) {
    CutPtLess pl(n);
    return make_cut(pl);
  }

  Cut ptIn(double n, double m) {
    if(n > m) swap(n , m);
    return ptGtr(n) & ptLess(m);
  }

  Cut massGtr(double n) {
    CutMassGtr mg(n);
    return make_cut(mg);
  }

  Cut massLess(double n) {
    CutMassLess ml(n);
    return make_cut(ml);
  }

  Cut massIn(double n, double m) {
    if(n > m) swap(n , m);
    return massGtr(n) & massLess(m);
  }

  Cut rapGtr(double n) {
    CutRapGtr rg(n);
    return make_cut(rg);
  }

  Cut rapLess(double n) {
    CutRapLess rl(n);
    return make_cut(rl);
  }

  Cut rapIn(double n, double m) {
    if(n > m) swap(n , m);
    return rapGtr(n) & rapLess(m);
  }

  Cut etaGtr(double n) {
    CutEtaGtr eg(n);
    return make_cut(eg);
  }

  Cut etaLess(double n) {
    CutEtaLess el(n);
    return make_cut(el);
  }

  Cut etaIn(double n, double m) {
    if(n > m) swap(n , m);
    return etaGtr(n) & etaLess(m);
  }



    //////////////
    /// Combiners

  CutsOr::CutsOr(const Cut c1, const Cut c2)
    : cut1(c1), cut2(c2) { }

  bool CutsOr::cut(const Cuttable& o) const
    {return cut1->cut(o) || cut2->cut(o);}

  CutsAnd::CutsAnd(const Cut c1, const Cut c2)
    : cut1(c1), cut2(c2) { }

  bool CutsAnd::cut(const Cuttable& o) const
    {return cut1->cut(o) && cut2->cut(o);}

  CutInvert::CutInvert(const Cut c1)
    : poscut(c1) { }

  bool CutInvert::cut(const Cuttable& o) const
    {return !poscut->cut(o);}

  CutsXor::CutsXor(const Cut c1, const Cut c2)
    : cut1(c1), cut2(c2) { }

  bool CutsXor::cut(const Cuttable& o) const
    {return !(cut1->cut(o) && cut2->cut(o)) && (cut1->cut(o) || cut2->cut(o));}

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


    /// pT
  CutPtGtr::CutPtGtr(const double pt_lowerlim)
      : low_(pt_lowerlim) { }

  bool CutPtGtr::cut(const Cuttable & o) const
    {return o.pT() >= low_;}

  CutPtLess::CutPtLess(const double pt_upperlim)
      : high_(pt_upperlim) { }

  bool CutPtLess::cut(const Cuttable & o) const
    {return o.pT() < high_;}


    /// Mass
  CutMassGtr::CutMassGtr(const double mass_lowerlim)
      : low_(mass_lowerlim) { }

  bool CutMassGtr::cut(const Cuttable & o) const
    {return o.m() >= low_;}

  CutMassLess::CutMassLess(const double mass_upperlim)
      : high_(mass_upperlim) { }

  bool CutMassLess::cut(const Cuttable & o) const
    {return o.m() < high_;}


    /// Rapidity
  CutRapGtr::CutRapGtr(const double rap_lowerlim)
      : low_(rap_lowerlim) { }

  bool CutRapGtr::cut(const Cuttable & o) const
    {return o.y() >= low_;}

  CutRapLess::CutRapLess(const double rap_upperlim)
      : high_(rap_upperlim) { }

  bool CutRapLess::cut(const Cuttable & o) const
    {return o.y() < high_;}


    /// Pseudorapidity
  CutEtaGtr::CutEtaGtr(const double eta_lowerlim)
      : low_(eta_lowerlim) { }

  bool CutEtaGtr::cut(const Cuttable & o) const
    {return o.eta() >= low_;}

  CutEtaLess::CutEtaLess(const double eta_upperlim)
      : high_(eta_upperlim) { }

  bool CutEtaLess::cut(const Cuttable & o) const
    {return o.eta() < high_;}




}
