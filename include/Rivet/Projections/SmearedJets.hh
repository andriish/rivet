// -*- C++ -*-
#ifndef RIVET_SmearedJets_HH
#define RIVET_SmearedJets_HH

#include "Rivet/Jet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "boost/function.hpp"

namespace Rivet {


  // double PARTICLE_FN0(const Particle& p) { return 0; }
  // double PARTICLE_FN1(const Particle& p) { return 1; }
  // double P4_FN0(const FourMomentum& p) { return 0; }
  // double P4_FN1(const FourMomentum& p) { return 1; }
  // Particle PARTICLE_SMEAR_IDENTITY(const Particle& p) { return p; }

  double JET_EFF_ZERO(const Jet& p) { return 0; }
  double JET_EFF_ONE(const Jet& p) { return 1; }
  Jet JET_SMEAR_IDENTITY(const Jet& j) { return j; }
  double rand01() { return rand() / (double)RAND_MAX; }

  double JET_BTAG_PERFECT(const Jet& j) { return j.bTagged() ? 1 : 0; }
  double JET_BTAG_ATLAS_RUN1C(const Jet& j) {
    if (j.bTagged()) return 0.80*tanh(0.003*j.pT()/GeV)*(30/(1+0.086*j.pT()/GeV));
    if (j.cTagged()) return 0.20*tanh(0.02*j.pT()/GeV)*(1/(1+0.0034*j.pT()/GeV));
    return 0.002 + 7.3e-6*j.pT()/GeV;
  }

  double JET_CTAG_PERFECT(const Jet& j) { return j.cTagged() ? 1 : 0; }


  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  class SmearedJets : public JetAlg {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with efficiency and smearing function args
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2DFN, typename J2JFN>
    SmearedJets(const JetAlg& ja, J2DFN jetEffFn,
                J2JFN jetSmearFn=JET_SMEAR_IDENTITY,
                J2DFN bTagEffFn=JET_BTAG_PERFECT, J2DFN cTagEffFn=JET_CTAG_PERFECT)
      : _jetEffFn(jetEffFn), _bTagEffFn(cTagEffFn), _cTagEffFn(cTagEffFn), _jetSmearFn(jetSmearFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }


    /// @brief Constructor with all-mandatory efficiency function args and no smearing
    /// The jet reconstruction efficiency is mandatory; the smearing and tagging functions are optional
    template <typename J2DFN>
    SmearedJets(const JetAlg& ja, J2DFN jetEffFn, J2DFN bTagEffFn, J2DFN cTagEffFn=JET_CTAG_PERFECT)
      : _jetEffFn(jetEffFn), _bTagEffFn(cTagEffFn), _cTagEffFn(cTagEffFn)
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


    // /// @todo Remove from JetAlg API? I *think* calc() doesn't work well on projections which chain others
    // void calc(const Particles& constituents, const Particles& tagparticles=Particles()) {
    //   ///
    // }


    /// Perform the jet finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Jets& truthjets = applyProjection<JetAlg>(e, "TruthJets").jetsByPt();
      _recojets.clear(); _recojets.reserve(truthjets.size());
      foreach (const Jet& j, truthjets) {
        const double jeff = (_jetEffFn) ? _jetEffFn(j) : 1;
        MSG_DEBUG("Efficiency of jet " << j.mom() << " = " << 100*jeff << "%");
        if (jeff == 0) continue; //< no need to roll expensive dice
        if (jeff == 1 || jeff < rand01()) {
          /// @todo Modify constituent particle vectors for consistency
          /// @todo Set a null PseudoJet if the Jet is smeared?
          _recojets.push_back(_jetSmearFn ? _jetSmearFn(j) : j); //< smearing
        }
      }
      // Tagging efficiencies
      foreach (Jet& j, _recojets) {
        const double beff = (_bTagEffFn) ? _bTagEffFn(j) : 1;
        const bool btag = beff == 1 || (beff != 0 && beff < rand01());
        /// @todo Remove b-tags if needed, and add a dummy one if needed
        //if (!btag && j.bTagged()) erase(remove_if(...)...);
        //if (btag && !j.bTagged()) j.tags().push_back(Particle(...));
        const double ceff = (_cTagEffFn) ? _cTagEffFn(j) : 1;
        const bool ctag = ceff == 1 || (ceff != 0 && ceff < rand01());
        /// @todo Remove c-tags if needed, and add a dummy one if needed
        //if (!ctag && j.cTagged()) erase(remove_if(...)...);
        //if (ctag && !j.cTagged()) j.tags().push_back(Particle(...));
      }
    }


    /// Return the full jet list for the JetAlg methods to use
    Jets _jets() const { return _recojets; }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _recojets.clear(); }


  private:

    Jets _recojets;
    boost::function<double(const Jet&)> _jetEffFn, _bTagEffFn, _cTagEffFn;
    boost::function<Jet(const Jet&)> _jetSmearFn;

  };


}


#endif
