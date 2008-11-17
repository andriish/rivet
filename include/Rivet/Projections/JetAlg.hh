// -*- C++ -*-
#ifndef RIVET_JetAlg_HH
#define RIVET_JetAlg_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

namespace Rivet {

  // Jets list sorting functions
  // Note that sorting is inverted, so that highest E is at the front of the list
  inline bool cmpJetsByPt(const Jet& a, const Jet& b) {
    return a.ptSum() > b.ptSum();
  }
  inline bool cmpJetsByEt(const Jet& a, const Jet& b) {
    return a.EtSum() > b.EtSum();
  }
  inline bool cmpJetsByE(const Jet& a, const Jet& b) {
    return a.momentum().E() > b.momentum().E();
  }


  
  /// Abstract base class for projections which can return a set of {@link Jet}s.
  class JetAlg : public Projection {
    
  public:

    /// Clone on the heap.
    virtual const Projection* clone() const = 0;

    /// Destructor
    virtual ~JetAlg() { }

    /// Get the jets (unordered).
    virtual Jets getJets(double ptmin=0.0) const = 0;

    /// Get the jets, ordered by supplied sorting function object.
    template <typename F>
    Jets jets(double ptmin=0.0, F sorter=0) const {
      Jets js = getJets(ptmin);
      if (sorter != 0) {
        std::sort(js.begin(), js.end(), sorter);
      }
      return js;
    }

    /// Get the jets, ordered by \f$ p_T \f$.
    Jets jetsByPt(double ptmin=0.0) const {
      return jets(ptmin, cmpJetsByPt);
    }

    /// Get the jets, ordered by \f$ E \f$.
    Jets jetsByE(double ptmin=0.0) const {
      return jets(ptmin, cmpJetsByE);
    }

    /// Get the jets, ordered by \f$ E_T \f$.
    Jets jetsByEt(double ptmin=0.0) const {
      return jets(ptmin, cmpJetsByEt);
    }


  public:

    typedef Jet entity_type;
    typedef Jets collection_type; 

    /// Template-usable interface common to FinalState.
    collection_type entities() const { return getJets(); }

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
