#include "Rivet/Cuts.hh"

using namespace std;
namespace Rivet {

    ///////////////////
    /// Access pointers

  CutPtr PtGtr(double n) {
    CutPtGtr pg(n);
    return make_cut(pg);
  }

  CutPtr PtLess(double n) {
    CutPtLess pl(n);
    return make_cut(pl);
  }

    //////////////
    /// Combiners

  CutsOr::CutsOr(const CutPtr c1, const CutPtr c2)
    : cut1(c1), cut2(c2) { }

  bool CutsOr::cut(const Cuttable& o) const
    {return cut1->cut(o) || cut2->cut(o);}

  CutsAnd::CutsAnd(const CutPtr c1, const CutPtr c2)
    : cut1(c1), cut2(c2) { }

  bool CutsAnd::cut(const Cuttable& o) const
    {return cut1->cut(o) && cut2->cut(o);}

  CutInvert::CutInvert(const CutPtr c1)
    : poscut(c1) { }

  bool CutInvert::cut(const Cuttable& o) const
    {return !poscut->cut(o);}

  CutsXor::CutsXor(const CutPtr c1, const CutPtr c2)
    : cut1(c1), cut2(c2) { }

  bool CutsXor::cut(const Cuttable& o) const
    {return !(cut1->cut(o) && cut2->cut(o)) && (cut1->cut(o) || cut2->cut(o));}

    ////////////
    ///Operators

  CutPtr operator & (const CutPtr aptr, const CutPtr bptr) {
    return make_cut(CutsAnd(aptr,bptr));
  }

  CutPtr operator | (const CutPtr aptr, const CutPtr bptr) {
    return make_cut(CutsOr(aptr,bptr));
  }

  CutPtr operator ~ (const CutPtr cptr) {
    return make_cut(CutInvert(cptr));
  }

  CutPtr operator ^ (const CutPtr aptr, const CutPtr bptr) {
    return make_cut(CutsXor(aptr,bptr));
  }

  double Cuttable::pT() const { assert(false); }

  CutPtGtr::CutPtGtr(const double pt_lowerlim)
      : low_(pt_lowerlim) { }

  bool CutPtGtr::cut(const Cuttable & o) const
    {return o.pT() > low_;}

  CutPtLess::CutPtLess(const double pt_upperlim)
      : high_(pt_upperlim) { }

  bool CutPtLess::cut(const Cuttable & o) const
    {return o.pT() < high_;}

}
