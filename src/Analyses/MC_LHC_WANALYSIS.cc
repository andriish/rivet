// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class MC_LHC_WANALYSIS : public Analysis {
  public:
  
    /// Default constructor
    MC_LHC_WANALYSIS() : Analysis("MC_LHC_WANALYSIS") 
    {
      //
    }
    

    /// @name Analysis methods
    /// @todo change "Weights" to differential cross sections once histos normalised to cross-section.
    //@{

    void init() { 
      const ChargedFinalState cfs;
      addProjection(cfs, "CFS");
      /// @todo Handle muon-decay Ws as well
      const WFinder wf(-MAXRAPIDITY, MAXRAPIDITY, 0.0*GeV, ELECTRON, 30.0*GeV, 110.0*GeV, 0.2);
      addProjection(wf, "WF");
      FastJets fastjets(wf.remainingFinalState(), FastJets::KT, 0.7);
      addProjection(fastjets, "Jets");

      _hist_chargemulti = bookHistogram1D("n-ch", 50, -0.5, 199.5);
      _hist_chargept = bookHistogram1D("pt-ch", 25, 0, 25);
      /// @todo Use profile plots instead:
      _hist_chargemeanpt = bookHistogram1D("ptavg-ch", 25, 0, 10);
      /// @todo Use profile plots instead (and this isn't really an "RMS")
      _hist_chargermspt = bookHistogram1D("ptrms-ch", 25, 0, 10);
      _hist_wcount = bookHistogram1D("n-w", 16, -0.5, 15.5);
      _hist_wpt = bookHistogram1D("pt-w", 25, 0, 25);
      _hist_wlogpt = bookHistogram1D("logpt-w", 30, 0, 6);
      _hist_weta = bookHistogram1D("eta-w", 36, -6, 6);
      _hist_wphi = bookHistogram1D("phi-w", 25, 0, TWOPI);
      _hist_wmass = bookHistogram1D("m-w", 40, 60, 100);
      _hist_logmass = bookHistogram1D("logm-w", 20, 0, 10);
      _hist_jetcount = bookHistogram1D("n-jet", 16, -0.5, 15.5);
      _hist_jetpt = bookHistogram1D("pt-jet", 50, 20, 100);
      _hist_jetlogpt = bookHistogram1D("logpt-jet", 20, 0, 20);
    }
    
    
    void analyze(const Event& event) {
      const double weight = event.weight();
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      const WFinder& wf = applyProjection<WFinder>(event, "WF");
      const FastJets& fastjets = applyProjection<FastJets>(event, "Jets");
      const Jets jets = fastjets.jetsByPt();
    
      // Charged particles part
      _hist_chargemulti->fill(cfs.particles().size(), weight);
      double meanpt(0), rmspt(0);
      foreach (const Particle& p, cfs.particles()) {
        const double pT = p.momentum().pT();
        _hist_chargept->fill(pT/GeV, weight);
        meanpt += pT;
        rmspt += pT*pT;
      }
      meanpt = meanpt / cfs.particles().size();
      _hist_chargemeanpt->fill(meanpt/GeV, weight);
      rmspt = sqrt(rmspt / cfs.particles().size());
      _hist_chargermspt->fill(rmspt/GeV, weight);
      
      // W part
      _hist_wcount->fill(wf.particles().size(), weight);
      foreach (const Particle& wp, wf.particles()) {
        const double pT = wp.momentum().pT();
        _hist_wpt->fill(pT/GeV, weight);
        _hist_wlogpt->fill(log(pT/GeV), weight);
        _hist_weta->fill(wp.momentum().pseudorapidity(), weight);
        _hist_wphi->fill(wp.momentum().azimuthalAngle(), weight);
        const double m = wp.momentum().mass();
        _hist_wmass->fill(m/GeV, weight);
        _hist_wlogmass->fill(log(m/GeV), weight);	
      }
      
      // Jet part
      _hist_jetcount->fill(fastjets.size(), weight);
      foreach(const Jet& j, fastjets.jetsByPt()) {
        const double pT = j.momentum().pT();
        _hist_jetpt->fill(pT/GeV, weight);
        _hist_jetlogpt->fill(log(pT/GeV), weight);
      }
    }
    
    
    void finalize() {
      ///@todo Obtain cross-sections from generator and normalise histos to them.
    }
    
    //@}
    
  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _hist_chargemulti;
    AIDA::IHistogram1D* _hist_chargept;
    AIDA::IHistogram1D* _hist_chargemeanpt;
    AIDA::IHistogram1D* _hist_chargermspt;
    AIDA::IHistogram1D* _hist_wcount;
    AIDA::IHistogram1D* _hist_wpt;
    AIDA::IHistogram1D* _hist_wlogpt;
    //AIDA::IHistogram1D* _hist_wpthigh;
    //AIDA::IHistogram1D* _hist_wlogpthigh;
    AIDA::IHistogram1D* _hist_weta;
    AIDA::IHistogram1D* _hist_wphi;
    AIDA::IHistogram1D* _hist_wmass;
    AIDA::IHistogram1D* _hist_wlogmass;
    AIDA::IHistogram1D* _hist_jetcount;
    AIDA::IHistogram1D* _hist_jetpt;
    AIDA::IHistogram1D* _hist_jetlogpt;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_LHC_WANALYSIS> plugin_MC_LHC_WANALYSIS;

}
