#include "Rivet/Cuts.hh"

using namespace std;
namespace Rivet {


    /////////////////////
    /// Low level structs


/*
  CutLess::CutLess(double n)
    :  val(n) {}

  bool CutLess::cut(double n) const
    {return n < val;}

  bool CutLess::cut(const Cuttable& o) const
    {return o.pT() < val;}

  CutLessEq::CutLessEq(double n)
    :  val(n) {}

  bool CutLessEq::cut(double n) const
    {return n <= val;}

  bool CutLessEq::cut(const Cuttable& o) const
    {return o.pT() <= val;}

  CutMore::CutMore(double n)
    :  val(n) {}

  bool CutMore::cut(double n) const
    {return n > val;}

  bool CutMore::cut(const Cuttable& o) const
    {return o.pT() > val;}

  CutMoreEq::CutMoreEq(double n)
    :  val(n) {}

  bool CutMoreEq::cut(double n) const
    {return n >= val;}

  bool CutMoreEq::cut(const Cuttable& o) const
    {return o.pT() >= val;}
*/


    ///////////////////
    /// Access pointers


/*
  CutPtr Less(double n) {
    CutLess l(n);
    return make_cut(l);
  }

  CutPtr More(double n) {
    CutMore m(n);
    return make_cut(m);
  }

  CutPtr LessEq(double n) {
    CutLessEq le(n);
    return make_cut(le);
  }

  CutPtr MoreEq(double n) {
    CutMoreEq me(n);
    return make_cut(me);
  }
*/

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

  //bool CutsOr::cut(double n) const
  //  {return cut1->cut(n) || cut2->cut(n);}

  bool CutsOr::cut(const Cuttable& o) const
    {return cut1->cut(o) || cut2->cut(o);}

  CutsAnd::CutsAnd(const CutPtr c1, const CutPtr c2)
    : cut1(c1), cut2(c2) { }

  //bool CutsAnd::cut(double n) const
  //  {return cut1->cut(n) && cut2->cut(n);}

  bool CutsAnd::cut(const Cuttable& o) const
    {return cut1->cut(o) && cut2->cut(o);}

  CutInvert::CutInvert(const CutPtr c1)
    : poscut(c1) { }

  //bool CutInvert::cut(double n) const
  //  {return !poscut->cut(n);}

  bool CutInvert::cut(const Cuttable& o) const
    {return !poscut->cut(o);}

  CutsXor::CutsXor(const CutPtr c1, const CutPtr c2)
    : cut1(c1), cut2(c2) { }

  //bool CutsXor::cut(double n) const
  //  {return !(cut1->cut(n) && cut2->cut(n)) && (cut1->cut(n) || cut2->cut(n));}

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

  //bool CutPtGtr::cut(double n) const
  //  {return n > low_;}

  CutPtLess::CutPtLess(const double pt_upperlim)
      : high_(pt_upperlim) { }

  bool CutPtLess::cut(const Cuttable & o) const
    {return o.pT() < high_;}

  //bool CutPtLess::cut(double n) const
  //  {return n < high_;}

}
