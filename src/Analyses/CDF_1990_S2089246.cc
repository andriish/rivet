// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/TriggerCDFRun0Run1.hh"

namespace Rivet {


  /* @brief CDF pseudorapidity analysis
   * @author Andy Buckley
   */
  class CDF_1990_S2089246 : public Analysis {
  public:

    /// Constructor
    CDF_1990_S2089246()
      : Analysis("CDF_1990_S2089246")
    {
      setBeams(PROTON, ANTIPROTON);
    }


    /// @name Analysis methods
    //@{

    void init() {
      addProjection(TriggerCDFRun0Run1(), "Trigger");
      addProjection(ChargedFinalState(-3.5, 3.5), "CFS");
      addProjection(Beam(), "Beam");

      _hist_eta1800 = bookHistogram1D(3, 1, 1);
      _hist_eta630 = bookHistogram1D(4, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Trigger
      const bool trigger = applyProjection<TriggerCDFRun0Run1>(event, "Trigger").minBiasDecision();
      if (!trigger) vetoEvent;

      // Get final state and energy
      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS(); 
      const FinalState& fs = applyProjection<FinalState>(event, "CFS");

      // Loop over final state charged particles to fill eta histos
      const double weight = event.weight();
      foreach (const Particle& p, fs.particles()) {
        const double eta = p.momentum().pseudorapidity();
        if (fuzzyEquals(sqrtS/GeV, 630)) {
          _hist_eta630->fill(fabs(eta), weight);
        } else if (fuzzyEquals(sqrtS/GeV, 1800)) {
          _hist_eta1800->fill(fabs(eta), weight);
        }
      }
    }
    
    
    
    /// Finalize
    void finalize() {
      // Divide through by num events to get d<N>/d(eta) in bins
      scale(_hist_eta630, 1/sumOfWeights());
      scale(_hist_eta1800, 1/sumOfWeights());
    }
   
    //@}


  private:

    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _hist_eta630;
    AIDA::IHistogram1D* _hist_eta1800;
    //@}

  };
 
    

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_1990_S2089246> plugin_CDF_1990_S2089246;

}
