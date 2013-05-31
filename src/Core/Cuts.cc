#include "Rivet/Cuts.hh"

using namespace std;
namespace Rivet {

  CutLess::CutLess(double n)
    :  val(n) {}

  bool CutLess::cut(double n) const
    {return n < val;}

  CutLessEq::CutLessEq(double n)
    :  val(n) {}

  bool CutLessEq::cut(double n) const
    {return n <= val;}

  CutMore::CutMore(double n)
    :  val(n) {}

  bool CutMore::cut(double n) const
    {return n > val;}

  CutMoreEq::CutMoreEq(double n)
    :  val(n) {}

  bool CutMoreEq::cut(double n) const
    {return n >= val;}

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

  CutsOr::CutsOr(const CutPtr c1, const CutPtr c2)
    : cut1(c1), cut2(c2) { }

  bool CutsOr::cut(double n) const
    {return cut1->cut(n) || cut2->cut(n);}

  CutsAnd::CutsAnd(const CutPtr c1, const CutPtr c2)
    : cut1(c1), cut2(c2) { }

  bool CutsAnd::cut(double n) const
    {return cut1->cut(n) && cut2->cut(n);}

  CutInvert::CutInvert(const CutPtr c1)
    : poscut(c1) { }

  bool CutInvert::cut(double n) const
    {return !poscut->cut(n);}

  CutsXor::CutsXor(const CutPtr c1, const CutPtr c2)
    : cut1(c1), cut2(c2) { }

  bool CutsXor::cut(double n) const
    {return !(cut1->cut(n) && cut2->cut(n)) && (cut1->cut(n) || cut2->cut(n));}

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



}
