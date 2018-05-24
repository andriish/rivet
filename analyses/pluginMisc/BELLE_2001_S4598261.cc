// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief BELLE pi0 spectrum at Upsilon(4S)
  /// @author Peter Richardson
  class BELLE_2001_S4598261 : public Analysis {
  public:

    BELLE_2001_S4598261()
      : Analysis("BELLE_2001_S4598261")
    { }


    void init() {
      declare(UnstableFinalState(), "UFS");
      book(_histdSigDp ,1, 1, 1); // spectrum
      book(_histMult   ,2, 1, 1); // multiplicity
      book(_weightSum, "TMP/weightSum");
    }


    void analyze(const Event& e) {
      // Find the upsilons
      Particles upsilons;
      // First in unstable final state
      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");
      for (const Particle& p : ufs.particles())
        if (p.pid()==300553) upsilons.push_back(p);
      // Then in whole event if fails
      if (upsilons.empty()) {
        for (const GenParticle* p : Rivet::particles(e.genEvent())) {
          if (p->pdg_id() != 300553) continue;
          const GenVertex* pv = p->production_vertex();
          bool passed = true;
          if (pv) {
            /// @todo Use better looping
            for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ; pp != pv->particles_in_const_end() ; ++pp) {
              if ( p->pdg_id() == (*pp)->pdg_id() ) {
                passed = false;
                break;
              }
            }
          }
          if (passed) upsilons.push_back(Particle(p));
        }
      }

      // Find upsilons
      for (const Particle& p : upsilons) {
        _weightSum->fill();
        // Find the neutral pions from the decay
        vector<GenParticle *> pions;
        findDecayProducts(p.genParticle(), pions);
        const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        for (size_t ix=0; ix<pions.size(); ++ix) {
          const double pcm = cms_boost.transform(FourMomentum(pions[ix]->momentum())).p();
          _histdSigDp->fill(pcm);
        }
        _histMult->fill(0., pions.size());
      }
    }


    void finalize() {
      scale(_histdSigDp, 1./ *_weightSum);
      scale(_histMult  , 1./ *_weightSum);
    }


  private:

    //@{
    // count of weights
    CounterPtr _weightSum;
    /// Histograms
    Histo1DPtr _histdSigDp;
    Histo1DPtr _histMult;
    //@}


    void findDecayProducts(const GenParticle* p, vector<GenParticle*>& pions) {
      const GenVertex* dv = p->end_vertex();
      /// @todo Use better looping
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        const int id = (*pp)->pdg_id();
        if (id == 111) {
          pions.push_back(*pp);
        } else if ((*pp)->end_vertex())
          findDecayProducts(*pp, pions);
      }
    }


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2001_S4598261);

}
