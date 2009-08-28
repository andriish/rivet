// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"

namespace Rivet {


  /* @brief CDF pseudorapidity analysis
   * @author Andy Buckley
   */
  class CDF_1990_S2089246 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_1990_S2089246()
      : Analysis("CDF_1990_S2089246")
    {
      setBeams(PROTON, ANTIPROTON);
      addProjection(ChargedFinalState(-3.5, 3.5), "FS");
      addProjection(ChargedFinalState(-5.9, 5.9), "CFSAll");
      addProjection(Beam(), "Beam");
    }
    
    //@}


    /// @name Analysis methods
    //@{

    void init() {
      _hist_eta1800 = bookHistogram1D(3, 1, 1);
      _hist_eta630 = bookHistogram1D(4, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& event) {
      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
      const double weight = event.weight();
      
      // Minimum Bias trigger requirements from the BBC counters
      int n_trig_1 = 0;
      int n_trig_2 = 0;
      
      // Event selection based on tracks in VTPC (time projection chambers)
      // Require at least 4 tracks with at least one in each of the forward
      // and backward hemispheres
      int n_backward = 0;
      int n_forward = 0;
      
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFSAll");
      foreach (const Particle& p, cfs.particles()) {
        double eta = p.momentum().pseudorapidity();
        if (inRange(eta, -5.9, -3.2)) n_trig_1++;
        else if (inRange(eta, 3.2, 5.9)) n_trig_2++;
        
        if (inRange(eta, -3.0, 0.0)) n_backward++;
        else if (inRange(eta, 0.0, 3.0)) n_forward++;
      }
      
      // Require at least one coincidence hit in both BBC counters
      if (!n_trig_1 || !n_trig_2) vetoEvent; 
      getLog() << Log::DEBUG << "Trigger 1: " << n_trig_1 << " Trigger 2: " << n_trig_2 << endl;
      
      // Further event selection cut
      if ( (n_backward+n_forward < 4) || (n_backward*n_forward < 1) ) vetoEvent;
      getLog() << Log::DEBUG << " Num. forward: " << n_forward  << " Num. backward: " << n_backward << endl;
      
      // Loop over final state charged particles 
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
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
