// -*- C++ -*-
#ifndef RIVET_JetAlg_HH
#define RIVET_JetAlg_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

namespace Rivet {


  inline bool cmpMomByPt(const FourMomentum& a, const FourMomentum& b) {
    return a.pT() > b.pT();
  }
  inline bool cmpMomByAscPt(const FourMomentum& a, const FourMomentum& b) {
    return a.pT() < b.pT();
  }
  inline bool cmpMomByEt(const FourMomentum& a, const FourMomentum& b) {
    return a.Et() > b.Et();
  }
  inline bool cmpMomByAscEt(const FourMomentum& a, const FourMomentum& b) {
    return a.Et() < b.Et();
  }
  inline bool cmpMomByE(const FourMomentum& a, const FourMomentum& b) {
    return a.E() > b.E();
  }
  inline bool cmpMomByAscE(const FourMomentum& a, const FourMomentum& b) {
    return a.E() < b.E();
  }
  inline bool cmpMomByDescPseudorapidity(const FourMomentum& a, const FourMomentum& b) {
    return a.pseudorapidity() > b.pseudorapidity();
  }
  inline bool cmpMomByAscPseudorapidity(const FourMomentum& a, const FourMomentum& b) {
    return a.pseudorapidity() < b.pseudorapidity();
  }
  inline bool cmpMomByDescAbsPseudorapidity(const FourMomentum& a, const FourMomentum& b) {
    return fabs(a.pseudorapidity()) > fabs(b.pseudorapidity());
  }
  inline bool cmpMomByAscAbsPseudorapidity(const FourMomentum& a, const FourMomentum& b) {
    return fabs(a.pseudorapidity()) < fabs(b.pseudorapidity());
  }
  inline bool cmpMomByDescRapidity(const FourMomentum& a, const FourMomentum& b) {
    return a.rapidity() > b.rapidity();
  }
  inline bool cmpMomByAscRapidity(const FourMomentum& a, const FourMomentum& b) {
    return a.rapidity() < b.rapidity();
  }
  inline bool cmpMomByDescAbsRapidity(const FourMomentum& a, const FourMomentum& b) {
    return fabs(a.rapidity()) > fabs(b.rapidity());
  }
  inline bool cmpMomByAscAbsRapidity(const FourMomentum& a, const FourMomentum& b) {
    return fabs(a.rapidity()) < fabs(b.rapidity());
  }

  
  /// Abstract base class for projections which can return a set of {@link Jet}s.
  class JetAlg : public Projection {
    
  public:

    /// Clone on the heap.
    virtual const Projection* clone() const = 0;

    /// Destructor
    virtual ~JetAlg() { }

    virtual Jets jets(double ptmin=0.0) const = 0;

    /// Get the jets, ordered by supplied sorting function object.
    template <typename F>
      Jets jets(F sorter, double ptmin=0.0) const {
      Jets js = jets(ptmin);
      if (sorter != 0) {
        std::sort(js.begin(), js.end(), sorter);
      }
      return js;
    }

    /// Get the jets, ordered by \f$ p_T \f$.
    Jets jetsByPt(double ptmin=0.0) const {
      return jets(cmpJetsByPt, ptmin);
    }

    /// Get the jets, ordered by \f$ E \f$.
    Jets jetsByE(double ptmin=0.0) const {
      return jets(cmpJetsByE, ptmin);
    }

    /// Get the jets, ordered by \f$ E_T \f$.
    Jets jetsByEt(double ptmin=0.0) const {
      return jets(cmpJetsByEt, ptmin);
    }


  public:

    /// Number of jets.
    virtual size_t size() const = 0;

    /// Clear the projection
    virtual void reset() = 0;

    typedef Jet entity_type;
    typedef Jets collection_type; 

    /// Template-usable interface common to FinalState.
    collection_type entities() const { return jets(); }

    /// Do the calculation locally (no caching).
    void calc(const ParticleVector& ps);


  protected:

    /// Perform the projection on the Event.
    virtual void project(const Event& e) = 0;

    /// Compare projections.
    virtual int compare(const Projection& p) const = 0;

  };

}

#endif
