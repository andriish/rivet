#ifndef RIVET_Cuts_HH
#define RIVET_Cuts_HH

#include <iostream>
#include "boost/smart_ptr.hpp"

using namespace std;

namespace Rivet {

/// Base class for all cut functors

struct Cut {
    virtual bool cut(double n) const = 0;
    virtual ~Cut() {}
};

typedef boost::shared_ptr<Cut> CutPtr; // Cut;

template <typename T>
CutPtr make_cut(T t) {
    return boost::shared_ptr<T>(new T(t));
}

/// Various cut functors operating on doubles, not actually
/// called by user.

struct CutLess : public Cut {
    CutLess(double n);
    bool cut(double n) const;
    double val;
};
struct CutLessEq : public Cut {
    CutLessEq(double n);
    bool cut(double n) const;
    double val;
};
struct CutMore : public Cut {
    CutMore(double n);
    bool cut(double n) const;
    double val;
};
struct CutMoreEq : public Cut {
    CutMoreEq(double n);
    bool cut(double n) const;
    double val;
};

/// These pointers are used to access the cuts

CutPtr Less(double n);

CutPtr More(double n);

CutPtr LessEq(double n);

CutPtr MoreEq(double n);

/// AND, OR, NOT, and XOR objects for combining cuts

struct CutsOr : public Cut {
    CutsOr(const CutPtr c1, const CutPtr c2);
    bool cut(double n) const;
    const CutPtr cut1;
    const CutPtr cut2;
};

struct CutsAnd : public Cut {
    CutsAnd(const CutPtr c1, const CutPtr c2);
    bool cut(double n) const;
    const CutPtr cut1;
    const CutPtr cut2;
};

struct CutInvert : public Cut {
    CutInvert(const CutPtr c1);
    bool cut(double n) const;
    const CutPtr poscut;
};

struct CutsXor : public Cut {
    CutsXor(const CutPtr c1, const CutPtr c2);
    bool cut(double n) const;
    const CutPtr cut1;
    const CutPtr cut2;
};

/// operator &, operator |, operator ~, and operator ^ overloads

CutPtr operator & (const CutPtr aptr, const CutPtr bptr);

CutPtr operator | (const CutPtr aptr, const CutPtr bptr);

CutPtr operator ~ (const CutPtr cptr);

CutPtr operator ^ (const CutPtr aptr, const CutPtr bptr);

/////

/*

class CutHolder {
public:
    CutHolder(const CutPtr cut) : cut_(cut) {}
    bool apply(int n) const {
        return cut_->cut(n);
    }
private:
    const CutPtr cut_;
};

*/
}

#endif

