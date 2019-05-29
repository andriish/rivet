#ifndef RIVET_RIVETYODA_HH
#define RIVET_RIVETYODA_HH

#include "Rivet/Config/RivetCommon.hh"
#include "YODA/AnalysisObject.h"
#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/Scatter1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Scatter3D.h"

#include <map>
#include <valarray>

namespace YODA {
  typedef std::shared_ptr<YODA::AnalysisObject> AnalysisObjectPtr;
}


namespace Rivet {


  class AnalysisObjectWrapper {
  public:
    virtual ~AnalysisObjectWrapper() {}

    virtual YODA::AnalysisObject* operator->() = 0;
    virtual YODA::AnalysisObject* operator->() const = 0;
    virtual const YODA::AnalysisObject & operator*() const = 0;

    /// @todo Rename to setActive(idx)
    virtual void setActiveWeightIdx(unsigned int iWeight) = 0;

    bool operator ==(const AnalysisObjectWrapper& p) { return (this == &p); }
    bool operator !=(const AnalysisObjectWrapper& p) { return (this != &p); }

  protected:
    /// @todo do we need this?
    // virtual void reset() = 0;
  };

  /// @todo
  /// implement scatter1dptr and scatter2dptr here
  /// these need to be multi-weighted eventually.
  /*
    class Scatter1DPtr : public AnalysisObjectPtr {
    public:
    Scatter1DPtr() : _persistent() { }

    Scatter1DPtr(size_t len_of_weightvec, const YODA::Scatter1D& p) {
    for (size_t m = 0; m < len_of_weightvec; ++m)
    _persistent.push_back(make_shared<YODA::Scatter1D>(p));
    }

    bool operator!() const { return !_persistent; }
    explicit operator bool() const { return bool(_persistent); }

    YODA::Scatter1D* operator->() { return _persistent.get(); }

    YODA::Scatter1D* operator->() const { return _persistent.get(); }

    YODA::Scatter1D & operator*() { return *_persistent; }

    const YODA::Scatter1D & operator*() const { return *_persistent; }

    protected:
    vector<YODA::Scatter1DPtr> _persistent;
    };

    class Scatter2DPtr : public AnalysisObjectPtr {
    public:
    Scatter2DPtr(size_t len_of_weightvec, const YODA::Scatter2D& p) {
    for (size_t m = 0; m < len_of_weightvec; ++m)
    _persistent.push_back(make_shared<YODA::Scatter2D>(p));
    }

    Scatter2DPtr() : _persistent() { }

    bool operator!() { return !_persistent; }
    explicit operator bool() { return bool(_persistent); }

    YODA::Scatter2D* operator->() { return _persistent.get(); }

    YODA::Scatter2D* operator->() const { return _persistent.get(); }

    YODA::Scatter2D & operator*() { return *_persistent; }

    const YODA::Scatter2D & operator*() const { return *_persistent; }

    protected:
    vector<YODA::Scatter2DPtr> _persistent;
    };

    class Scatter3DPtr : public AnalysisObjectPtr {
    public:
    Scatter3DPtr(size_t len_of_weightvec, const YODA::Scatter3D& p) {
    for (size_t m = 0; m < len_of_weightvec; ++m)
    _persistent.push_back(make_shared<YODA::Scatter3D>(p));
    }

    Scatter3DPtr() : _persistent() { }

    bool operator!() { return !_persistent; }
    explicit operator bool() { return bool(_persistent); }

    YODA::Scatter3D* operator->() { return _persistent.get(); }

    YODA::Scatter3D* operator->() const { return _persistent.get(); }

    YODA::Scatter3D & operator*() { return *_persistent; }

    const YODA::Scatter3D & operator*() const { return *_persistent; }

    protected:
    vector<YODA::Scatter3DPtr> _persistent;
    };
  */


  class MultiweightAOWrapper : public AnalysisObjectWrapper {

  public:
    using Inner = YODA::AnalysisObject;

    virtual void newSubEvent() = 0;

    virtual void pushToPersistent(const vector<std::valarray<double> >& weight) = 0;

    virtual YODA::AnalysisObjectPtr activeYODAPtr() const = 0;
  };


  using Weight = double;

  template <class T>
  using Fill = pair<typename T::FillType, Weight>;

  template <class T>
  using Fills = multiset<Fill<T>>;


  // TODO TODO TODO
  // need to override the old fill method too!
  // otherwise we can't intercept existing fill calls in analysis code
  // TODO TODO TODO


  /// Wrappers for analysis objects to store all fills unaggregated, until collapsed
  template <class T>
  class TupleWrapper;

  template<>
  class TupleWrapper<YODA::Counter> : public YODA::Counter {
  public:
    typedef shared_ptr<TupleWrapper<YODA::Counter>> Ptr;
    TupleWrapper(const YODA::Counter & h) : YODA::Counter(h) {}
    // todo: do we need to deal with users using fractions directly?
    void fill( double weight=1.0, double fraction=1.0 ) {
      fills_.insert( {YODA::Counter::FillType(),weight} );
    }
    void reset() { fills_.clear(); }
    const Fills<YODA::Counter> & fills() const { return fills_; }
  private:
    // x / weight pairs
    Fills<YODA::Counter> fills_;
  };

  template<>
  class TupleWrapper<YODA::Histo1D> : public YODA::Histo1D {
  public:
    typedef shared_ptr<TupleWrapper<YODA::Histo1D>> Ptr;
    TupleWrapper(const YODA::Histo1D & h) : YODA::Histo1D(h) {}
    // todo: do we need to deal with users using fractions directly?
    void fill( double x, double weight=1.0, double fraction=1.0 ) {
      if ( std::isnan(x) ) throw YODA::RangeError("X is NaN");
      fills_.insert( { x , weight } );
    }
    void reset() { fills_.clear(); }
    const Fills<YODA::Histo1D> & fills() const { return fills_; }
  private:
    // x / weight pairs
    Fills<YODA::Histo1D> fills_;
  };

  template<>
  class TupleWrapper<YODA::Profile1D> : public YODA::Profile1D {
  public:
    typedef shared_ptr<TupleWrapper<YODA::Profile1D>> Ptr;
    TupleWrapper(const YODA::Profile1D & h) : YODA::Profile1D(h) {}
    // todo: do we need to deal with users using fractions directly?
    void fill( double x, double y, double weight=1.0, double fraction=1.0 ) {
      if ( std::isnan(x) ) throw YODA::RangeError("X is NaN");
      if ( std::isnan(y) ) throw YODA::RangeError("Y is NaN");
      fills_.insert( { YODA::Profile1D::FillType{x,y}, weight } );
    }
    void reset() { fills_.clear(); }
    const Fills<YODA::Profile1D> & fills() const { return fills_; }
  private:
    // x / weight pairs
    Fills<YODA::Profile1D> fills_;
  };


  template<>
  class TupleWrapper<YODA::Histo2D> : public YODA::Histo2D {
  public:
    typedef shared_ptr<TupleWrapper<YODA::Histo2D>> Ptr;
    TupleWrapper(const YODA::Histo2D & h) : YODA::Histo2D(h) {}
    // todo: do we need to deal with users using fractions directly?
    void fill( double x, double y, double weight=1.0, double fraction=1.0 ) {
      if ( std::isnan(x) ) throw YODA::RangeError("X is NaN");
      if ( std::isnan(y) ) throw YODA::RangeError("Y is NaN");
      fills_.insert( { YODA::Histo2D::FillType{x,y}, weight } );
    }
    void reset() { fills_.clear(); }
    const Fills<YODA::Histo2D> & fills() const { return fills_; }
  private:
    // x / weight pairs
    Fills<YODA::Histo2D> fills_;
  };

  template<>
  class TupleWrapper<YODA::Profile2D> : public YODA::Profile2D {
  public:
    typedef shared_ptr<TupleWrapper<YODA::Profile2D>> Ptr;
    TupleWrapper(const YODA::Profile2D & h) : YODA::Profile2D(h) {}
    // todo: do we need to deal with users using fractions directly?
    void fill( double x, double y, double z, double weight=1.0, double fraction=1.0 ) {
      if ( std::isnan(x) ) throw YODA::RangeError("X is NaN");
      if ( std::isnan(y) ) throw YODA::RangeError("Y is NaN");
      if ( std::isnan(z) ) throw YODA::RangeError("Z is NaN");
      fills_.insert( { YODA::Profile2D::FillType{x,y,z}, weight } );
    }
    void reset() { fills_.clear(); }
    const Fills<YODA::Profile2D> & fills() const { return fills_; }
  private:
    // x / weight pairs
    Fills<YODA::Profile2D> fills_;
  };

  template<>
  class TupleWrapper<YODA::Scatter1D> : public YODA::Scatter1D {
  public:
    typedef shared_ptr<TupleWrapper<YODA::Scatter1D>> Ptr;
    TupleWrapper(const YODA::Scatter1D & h) : YODA::Scatter1D(h) {}
  };

  template<>
  class TupleWrapper<YODA::Scatter2D> : public YODA::Scatter2D {
  public:
    typedef shared_ptr<TupleWrapper<YODA::Scatter2D>> Ptr;
    TupleWrapper(const YODA::Scatter2D & h) : YODA::Scatter2D(h) {}
  };

  template<>
  class TupleWrapper<YODA::Scatter3D> : public YODA::Scatter3D {
  public:
    typedef shared_ptr<TupleWrapper<YODA::Scatter3D>> Ptr;
    TupleWrapper(const YODA::Scatter3D & h) : YODA::Scatter3D(h) {}
  };



  template <class T>
  class Wrapper : public MultiweightAOWrapper {
    friend class Analysis;

  public:

    using Inner = T;
    /* @todo
     * some things are not really well-defined here
     * for instance: fill() in the finalize() method and integral() in
     * the analyze() method.
     */

    Wrapper() = default;

    Wrapper(const vector<string>& weightnames, const T & p);

    ~Wrapper();

    typename T::Ptr active() const;

    /* @todo this probably need to loop over all? */
    bool operator!() const { return !_active; } // Don't use active() here, assert will catch

    explicit operator bool() const { return static_cast<bool>(_active); } // Don't use active() here, assert will catch

    T * operator->() { return active().get(); }

    T * operator->() const { return active().get(); }

    T & operator*() { return *active(); }

    const T & operator*() const { return *active(); }
    /* @todo
     * these need to be re-thought out.

     void reset() { active()->reset(); }
    */

    /* @todo
     * these probably need to loop over all?
     * do we even want to provide equality?
     */
    /* @todo
     * how about no.
     friend bool operator==(Wrapper a, Wrapper b){
     if (a._persistent.size() != b._persistent.size())
     return false;

     for (size_t i = 0; i < a._persistent.size(); i++) {
     if (a._persistent.at(i) != b._persistent.at(i)) {
     return false;
     }
     }

     return true;
     }

     friend bool operator!=(Wrapper a, Wrapper b){
     return !(a == b);
     }


     friend bool operator<(Wrapper a, Wrapper b){
     if (a._persistent.size() >= b._persistent.size())
     return false;

     for (size_t i = 0; i < a._persistent.size(); i++) {
     if (*(a._persistent.at(i)) >= *(b._persistent.at(i))) {
     return false;
     }
     }

     return true;
     }
    */


  private:
    void setActiveWeightIdx(unsigned int iWeight) {
      _active = _persistent.at(iWeight);
    }

    /* this is for dev only---we shouldn't need this in real runs. */
    void unsetActiveWeight() { _active.reset(); }

    void newSubEvent();

    virtual YODA::AnalysisObjectPtr activeYODAPtr() const { return _active; }

    const vector<typename T::Ptr> & persistent() const { return _persistent; }

    /* to be implemented for each type */
    void pushToPersistent(const vector<std::valarray<double> >& weight);

    /* M of these, one for each weight */
    vector<typename T::Ptr> _persistent;

    /* N of these, one for each event in evgroup */
    vector<typename TupleWrapper<T>::Ptr> _evgroup;

    typename T::Ptr _active;

    // do we need implicit cast?
    // operator typename T::Ptr () {
    //     return _active;
    // }

    friend class AnalysisHandler;
  };


  /// We need our own shared_ptr class, so we can dispatch -> and *
  /// all the way down to the inner YODA analysis objects
  ///
  /// TODO: provide remaining functionality that shared_ptr has (not needed right now)
  ///
  template <typename T>
  class rivet_shared_ptr {
  public:
    typedef T value_type;

    rivet_shared_ptr() = default;

    rivet_shared_ptr(decltype(nullptr)) : _p(nullptr) {}

    /// Convenience constructor, pass through to the Wrapper constructor
    rivet_shared_ptr(const vector<string>& weightNames, const typename T::Inner & p)
      : _p( make_shared<T>(weightNames, p) )
    {}

    template <typename U>
    rivet_shared_ptr(const shared_ptr<U> & p)
      : _p(p)
    {}

    template <typename U>
    rivet_shared_ptr(const rivet_shared_ptr<U> & p)
      : _p(p.get())
    {}

    // Goes right through to the active YODA object's members
    T & operator->()                            { return  *_p; }
    const T & operator->() const                { return  *_p; }

    // The active YODA object
    typename T::Inner & operator*()             { return **_p; }
    const typename T::Inner & operator*() const { return **_p; }

    bool operator!() const { return !_p || !(*_p);   }
    explicit operator bool()  const { return _p && bool(*_p); }

    template <typename U>
    bool operator==(const rivet_shared_ptr<U> & other) const {
      return _p == other._p;
    }

    template <typename U>
    bool operator!=(const rivet_shared_ptr<U> & other) const {
      return _p != other._p;
    }

    template <typename U>
    bool operator<(const rivet_shared_ptr<U> & other) const {
      return _p < other._p;
    }

    template <typename U>
    bool operator>(const rivet_shared_ptr<U> & other) const {
      return _p > other._p;
    }

    template <typename U>
    bool operator<=(const rivet_shared_ptr<U> & other) const {
      return _p <= other._p;
    }

    template <typename U>
    bool operator>=(const rivet_shared_ptr<U> & other) const {
      return _p >= other._p;
    }

    shared_ptr<T> get() const { return _p; }
  private:
    shared_ptr<T> _p;
  };



  // every object listed here needs a virtual fill method in YODA,
  // otherwise the Tuple fakery won't work.

  using MultiweightAOPtr = rivet_shared_ptr<MultiweightAOWrapper>;

  using Histo1DPtr   = rivet_shared_ptr<Wrapper<YODA::Histo1D>>;
  using Histo2DPtr   = rivet_shared_ptr<Wrapper<YODA::Histo2D>>;
  using Profile1DPtr = rivet_shared_ptr<Wrapper<YODA::Profile1D>>;
  using Profile2DPtr = rivet_shared_ptr<Wrapper<YODA::Profile2D>>;
  using CounterPtr   = rivet_shared_ptr<Wrapper<YODA::Counter>>;
  using Scatter1DPtr = rivet_shared_ptr<Wrapper<YODA::Scatter1D>>;
  using Scatter2DPtr = rivet_shared_ptr<Wrapper<YODA::Scatter2D>>;
  using Scatter3DPtr = rivet_shared_ptr<Wrapper<YODA::Scatter3D>>;

  using YODA::Counter;
  using YODA::Histo1D;
  using YODA::HistoBin1D;
  using YODA::Histo2D;
  using YODA::HistoBin2D;
  using YODA::Profile1D;
  using YODA::ProfileBin1D;
  using YODA::Profile2D;
  using YODA::ProfileBin2D;
  using YODA::Scatter1D;
  using YODA::Point1D;
  using YODA::Scatter2D;
  using YODA::Point2D;
  using YODA::Scatter3D;
  using YODA::Point3D;

  /// Function to get a map of all the refdata in a paper with the
  /// given @a papername.
  map<string, YODA::AnalysisObjectPtr> getRefData(const string& papername);

  /// @todo Also provide a Scatter3D getRefData() version?

  /// Get the file system path to the reference file for this paper.
  string getDatafilePath(const string& papername);


  /// Traits class to access the type of the AnalysisObject in the
  /// reference files.
  template<typename T> struct ReferenceTraits {};
  template<> struct ReferenceTraits<Counter> { typedef Counter RefT; };
  template<> struct ReferenceTraits<Scatter1D> { typedef Scatter1D RefT; };
  template<> struct ReferenceTraits<Histo1D> { typedef Scatter2D RefT; };
  template<> struct ReferenceTraits<Profile1D> { typedef Scatter2D RefT; };
  template<> struct ReferenceTraits<Scatter2D> { typedef Scatter2D RefT; };
  template<> struct ReferenceTraits<Histo2D> { typedef Scatter3D RefT; };
  template<> struct ReferenceTraits<Profile2D> { typedef Scatter3D RefT; };
  template<> struct ReferenceTraits<Scatter3D> { typedef Scatter3D RefT; };


  /// If @a dst and @a src both are of same subclass T, copy the
  /// contents of @a src into @a dst and return true. Otherwise return
  /// false.
  template <typename T>
  inline bool aocopy(YODA::AnalysisObjectPtr src, YODA::AnalysisObjectPtr dst) {
    shared_ptr<T> tsrc = dynamic_pointer_cast<T>(src);
    if ( !tsrc ) return false;
    shared_ptr<T> tdst = dynamic_pointer_cast<T>(dst);
    if ( !tdst ) return false;
    *tdst = *tsrc;
    return true;
  }

  /// If @a dst and @a src both are of same subclass T, add the
  /// contents of @a src into @a dst and return true. Otherwise return
  /// false.
  template <typename T>
  inline bool aoadd(YODA::AnalysisObjectPtr dst, YODA::AnalysisObjectPtr src, double scale) {
    shared_ptr<T> tsrc = dynamic_pointer_cast<T>(src);
    if ( !tsrc ) return false;
    shared_ptr<T> tdst = dynamic_pointer_cast<T>(dst);
    if ( !tdst ) return false;
    tsrc->scaleW(scale);
    *tdst += *tsrc;
    return true;
  }

  /// If @a dst is the same subclass as @a src, copy the contents of @a
  /// src into @a dst and return true. Otherwise return false.
  bool copyao(YODA::AnalysisObjectPtr src, YODA::AnalysisObjectPtr dst);

  /// If @a dst is the same subclass as @a src, scale the contents of
  /// @a src with @a scale and add it to @a dst and return true. Otherwise
  /// return false.
  bool addaos(YODA::AnalysisObjectPtr dst, YODA::AnalysisObjectPtr src, double scale);

  /// Check if two analysis objects have the same binning or, if not
  /// binned, are in other ways compatible.
  template <typename TPtr>
  inline bool bookingCompatible(TPtr a, TPtr b) {
    return a->sameBinning(*b);
  }
  inline bool bookingCompatible(CounterPtr, CounterPtr) {
    return true;
  }
  inline bool bookingCompatible(Scatter1DPtr a, Scatter1DPtr b) {
    return a->numPoints() == b->numPoints();
  }
  inline bool bookingCompatible(Scatter2DPtr a, Scatter2DPtr b) {
    return a->numPoints() == b->numPoints();
  }
  inline bool bookingCompatible(Scatter3DPtr a, Scatter3DPtr b) {
    return a->numPoints() == b->numPoints();
  }

}

#endif
