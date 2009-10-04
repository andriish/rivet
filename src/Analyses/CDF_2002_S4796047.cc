// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerCDFRun0Run1.hh"

namespace Rivet {

  /*
   * @brief CDF Run I charged multiplicity measurement
   * @author Hendrik Hoeth
   * 
   * This analysis measures the charged multiplicity distribution
   * in minimum bias events at two different center-of-mass energies:
   * \f$ \sqrt{s} = \f$ 630 and 1800 GeV.
   * 
   * Particles with c*tau > 10 mm are considered stable, i.e. they
   * are reconstructed and their decay products removed. Selection
   * cuts are |eta|<1 and pT>0.4 GeV.
   * 
   * 
   * @par Run conditions
   * 
   * @arg Two different beam energies: \f$ \sqrt{s} = \$f 630 & 1800 GeV
   * @arg Run with generic QCD events.
   * @arg Set particles with c*tau > 10 mm stable
   * 
   */
  class CDF_2002_S4796047 : public Analysis {
  public:

    /// Constructor
    CDF_2002_S4796047()
      : Analysis("CDF_2002_S4796047")
    { 
      setBeams(PROTON, ANTIPROTON);
    }


    /// @name Analysis methods
    //@{
    
    /// Book projections and histograms
    void init() {
      addProjection(TriggerCDFRun0Run1(), "Trigger");
      addProjection(Beam(), "Beam");
      const ChargedFinalState cfs(-1.0, 1.0, 0.4*GeV);
      addProjection(cfs, "FS");

      _hist_multiplicity_630  = bookHistogram1D(1, 1, 1);
      _hist_multiplicity_1800 = bookHistogram1D(2, 1, 1);
      _hist_pt_vs_multiplicity_630  = bookProfile1D(3, 1, 1);
      _hist_pt_vs_multiplicity_1800 = bookProfile1D(4, 1, 1);
    }
    
    
    /// Do the analysis
    void analyze(const Event& evt) {
      // Trigger
      const bool trigger = applyProjection<TriggerCDFRun0Run1>(evt, "Trigger").minBiasDecision();
      if (!trigger) vetoEvent;
      const double weight = evt.weight();

      // Get beam energy and tracks
      const double sqrtS = applyProjection<Beam>(evt, "Beam").sqrtS();
      const ChargedFinalState& fs = applyProjection<ChargedFinalState>(evt, "FS");
      const size_t numParticles = fs.particles().size();

      // Fill histos of charged multiplicity distributions
      if (fuzzyEquals(sqrtS, 630/GeV)) {
        _hist_multiplicity_630->fill(numParticles, weight);
      } 
      else if (fuzzyEquals(sqrtS, 1800/GeV)) {
        _hist_multiplicity_1800->fill(numParticles, weight);
      }

      // Fill histos for <pT> vs. charged multiplicity
      foreach (const Particle& p, fs.particles()) {
        const double pT = p.momentum().pT();
        if (fuzzyEquals(sqrtS, 630/GeV)) {
          _hist_pt_vs_multiplicity_630->fill(numParticles, pT/GeV, weight);
        }
        else if (fuzzyEquals(sqrtS, 1800/GeV)) {
          _hist_pt_vs_multiplicity_1800->fill(numParticles, pT/GeV, weight);
        }
      }
      
    }
    

    void finalize() {
      /// @todo Get cross-section from the generator
      normalize(_hist_multiplicity_630, 3.21167);
      normalize(_hist_multiplicity_1800, 4.19121);
    }

    //@}


  private:

    AIDA::IHistogram1D *_hist_multiplicity_630;
    AIDA::IHistogram1D *_hist_multiplicity_1800;
    AIDA::IProfile1D   *_hist_pt_vs_multiplicity_630 ;
    AIDA::IProfile1D   *_hist_pt_vs_multiplicity_1800;

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_2002_S4796047> plugin_CDF_2002_S4796047;

}
