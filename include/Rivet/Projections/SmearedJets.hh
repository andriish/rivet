// -*- C++ -*-
#ifndef RIVET_SmearedJets_HH
#define RIVET_SmearedJets_HH

#include "Rivet/Jet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "boost/function.hpp"

namespace Rivet {


  void JET_SMEAR_IDENTITY(Jet&) { cout << "JSMEAR!" << endl; }
  void PARTICLE_SMEAR_IDENTITY(Particle&) { }
  double JET_FN0(const Jet& p) { return 0; }
  double JET_FN1(const Jet& p) { cout << "JEFF1" << endl; return 1; }
  double PARTICLE_FN0(const Particle& p) { return 0; }
  double PARTICLE_FN1(const Particle& p) { return 1; }
  double P4_FN0(const FourMomentum& p) { return 0; }
  double P4_FN1(const FourMomentum& p) { return 1; }


  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  class SmearedJets : public JetAlg {
  public:

    /// @name Constructors etc.
    //@{

    // SmearedJets(const JetAlg& ja,
    // boost::function<double(const Jet&)>& jetEffFn, boost::function<void(Jet&)>& jetSmearFn=jet_smear_identity)
    // boost::function<double(const Jet&)>& bTagEffFn=fn1, boost::function<double(const Jet&)>& bTagMistagFn=fn0,
    // boost::function<double(const Jet&)>& cTagEffFn=fn1, boost::function<double(const Jet&)>& cTagMistagFn=fn0)
    //: _jetEffFn(jetEffFn), _jetSmearFn(jetSmearFn)
    // _bTagEffFn(bTagEffFn), _bTagMistagFn(bTagMistagFn),
    // _cTagEffFn(cTagEffFn), _cTagMistagFn(cTagMistagFn)


    /// Constructor with jet efficiency function arg
    template <typename J2DFN>
    SmearedJets(const JetAlg& ja, J2DFN jetEffFn)
      : _jetEffFn(jetEffFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }


    /// Constructor with jet efficiency and smearing function args
    template <typename J2DFN, typename J2VFN>
    SmearedJets(const JetAlg& ja, J2DFN jetEffFn, J2VFN jetSmearFn)
      : _jetEffFn(jetEffFn), _jetSmearFn(jetSmearFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new SmearedJets(*this);
    }

    //@}


    int compare(const Projection& p) const {
      /// @todo Maybe do the comparison via void* pointer comparisons, extracted before the function<> wrapping?
      //const SmearedJets& other = dynamic_cast<const SmearedJets&>(p);
      //return cmp(_jetEffFn, other._jetEffFn) || cmp(_jetSmearFn, other._jetSmearFn);
      return UNDEFINED;
    }


    /// @todo Move!!! Improve random numbers...
    double rand01() const {
      return rand() / (double)RAND_MAX;
    }


    /// @todo Remove from JetAlg API? I *think* calc() doesn't work well on projections which chain others
    void calc(const Particles& constituents, const Particles& tagparticles=Particles()) {
      ///
    }


    /// Perform the jet finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Jets& truthjets = applyProjection<JetAlg>(e, "TruthJets").jetsByPt();
      _recojets.clear(); _recojets.reserve(truthjets.size());
      foreach (const Jet& j, truthjets) {
        const double jeff = (_jetEffFn) ? _jetEffFn(j) : 1;
        MSG_DEBUG("Efficiency of jet " << j.mom() << " = " << 100*jeff << "%");
        if (jeff == 0) continue; //< remove; no need to roll expensive dice
        if (jeff == 1 || jeff < rand01()) _recojets.push_back(j);
      }
      // Smearing
      if (_jetSmearFn) {
        foreach (Jet& j, _recojets) {
          _jetSmearFn(j); //< smear the jet momentum
        }
      }
    }


    /// Return the full jet list for the JetAlg methods to use
    Jets _jets() const { return _recojets; }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _recojets.clear(); }


    /// @todo Modified tagger?


  private:

    Jets _recojets;
    std::function<double(const Jet&)> _jetEffFn;
    //, _bTagEffFn, _bTagMistagFn, _cTagEffFn, _cTagMistagFn;
    std::function<void(Jet&)> _jetSmearFn;

  };


}


#endif
