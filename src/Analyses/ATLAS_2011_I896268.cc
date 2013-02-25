// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief  J/Psi production at ATLAS
  class ATLAS_2011_I896268: public Analysis {
  public:

    /// Constructor
    ATLAS_2011_I896268()
      : Analysis("ATLAS_2011_I896268")
    {}


    /// @name Analysis methods
    //@{

    void init() {
      addProjection(UnstableFinalState(), "UFS");
      _nonPrRapHigh    = bookHistogram1D( 14, 1, 1);
      _nonPrRapMedHigh = bookHistogram1D( 13, 1, 1);
      _nonPrRapMedLow  = bookHistogram1D( 12, 1, 1);
      _nonPrRapLow     = bookHistogram1D( 11, 1, 1);
      _PrRapHigh    = bookHistogram1D( 18, 1, 1);
      _PrRapMedHigh = bookHistogram1D( 17, 1, 1);
      _PrRapMedLow  = bookHistogram1D( 16, 1, 1);
      _PrRapLow     = bookHistogram1D( 15, 1, 1);
      _IncRapHigh    = bookHistogram1D( 20, 1, 1);
      _IncRapMedHigh = bookHistogram1D( 21, 1, 1);
      _IncRapMedLow  = bookHistogram1D( 22, 1, 1);
      _IncRapLow     = bookHistogram1D( 23, 1, 1);
    }


    void analyze(const Event& e) {

      // Get event weight for histo filling
      const double weight = e.weight();


      // Final state of unstable particles to get particle spectra
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(e, "UFS");

      foreach (const Particle& p, ufs.particles()) {
        if(abs(p.pdgId())!=443) continue;
        HepMC::GenVertex* gv = p.genParticle().production_vertex();
        bool nonPrompt = false;
        if (gv) {
          foreach (const GenParticle* pi, Rivet::particles(gv, HepMC::ancestors)) {
            const PdgId pid2 = pi->pdg_id();
            if (PID::isHadron(pid2) && PID::hasBottom(pid2)) {
              nonPrompt = true;
              break;
            }
          }
        }
        double rapidity = p.momentum().rapidity();
        double xp = p.momentum().perp();

        if(rapidity<=2.4 and rapidity>2.){
          if(nonPrompt) _nonPrRapHigh->fill(xp, weight);
          else if(!nonPrompt) _PrRapHigh->fill(xp,weight);
          _IncRapHigh->fill(xp,weight);
        }
        else if(rapidity<=2. and rapidity>1.5){
          if(nonPrompt) _nonPrRapMedHigh->fill(xp,weight);
          else if(!nonPrompt) _PrRapMedHigh->fill(xp,weight);
          _IncRapMedHigh->fill(xp,weight);
        }
        else if(rapidity<=1.5 and rapidity>0.75){
          if(nonPrompt) _nonPrRapMedLow->fill(xp,weight);
          else if(!nonPrompt) _PrRapMedLow->fill(xp,weight);
          _IncRapMedLow->fill(xp,weight);
        }
        
        else if(rapidity<=0.75){
          if(nonPrompt) _nonPrRapLow->fill(xp,weight);
          else if(!nonPrompt) _PrRapLow->fill(xp,weight);
          _IncRapLow->fill(xp,weight);
        }
      }
    }


    /// Finalize
    void finalize() {
      double factor = crossSection()/nanobarn*0.0593;

      scale(_PrRapHigh     , factor/sumOfWeights());
      scale(_PrRapMedHigh  , factor/sumOfWeights());
      scale(_PrRapMedLow   , factor/sumOfWeights());
      scale(_PrRapLow      , factor/sumOfWeights());

      scale(_nonPrRapHigh     , factor/sumOfWeights());
      scale(_nonPrRapMedHigh  , factor/sumOfWeights());
      scale(_nonPrRapMedLow   , factor/sumOfWeights());
      scale(_nonPrRapLow      , factor/sumOfWeights());

      scale(_IncRapHigh     , 1000.*factor/sumOfWeights());
      scale(_IncRapMedHigh  , 1000.*factor/sumOfWeights());
      scale(_IncRapMedLow   , 1000.*factor/sumOfWeights());
      scale(_IncRapLow      , 1000.*factor/sumOfWeights());

    }

    //@}


  private:

    AIDA::IHistogram1D *_nonPrRapHigh;
    AIDA::IHistogram1D *_nonPrRapMedHigh;
    AIDA::IHistogram1D *_nonPrRapMedLow;
    AIDA::IHistogram1D *_nonPrRapLow;

    AIDA::IHistogram1D *_PrRapHigh;
    AIDA::IHistogram1D *_PrRapMedHigh;
    AIDA::IHistogram1D *_PrRapMedLow;
    AIDA::IHistogram1D *_PrRapLow;

    AIDA::IHistogram1D *_IncRapHigh;
    AIDA::IHistogram1D *_IncRapMedHigh;
    AIDA::IHistogram1D *_IncRapMedLow;
    AIDA::IHistogram1D *_IncRapLow;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I896268);

}
