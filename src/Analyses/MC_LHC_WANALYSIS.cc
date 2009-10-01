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
    MC_LHC_WANALYSIS()
      : Analysis("MC_LHC_WANALYSIS") 
    {
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

      _hist_chargemulti = bookHistogram1D("d01-x01-y01", 30, 0.5, 250.5);
      _hist_chargept = bookHistogram1D("d02-x01-y01", 32, 0., 25.);
      _hist_chargemeanpt = bookHistogram1D("d03-x01-y01", 25, 0., 10.);
      _hist_chargermspt = bookHistogram1D("d04-x01-y01", 32, 0., 10.);
      _hist_wcount = bookHistogram1D("d05-x01-y01", 30, 0., 15.);
      _hist_wpt = bookHistogram1D("d06-x01-y01", 32, 0., 25.);
      _hist_wlogpt = bookHistogram1D("d07-x01-y01", 32, 0., 6.);
      _hist_weta = bookHistogram1D("d08-x01-y01", 32, -6., 6.);
      _hist_wphi = bookHistogram1D("d09-x01-y01", 32, 0., 6.4);
      _hist_wmass = bookHistogram1D("d10-x01-y01", 40, 60., 100.);
      _hist_wlogmass = bookHistogram1D("d11-x01-y01", 32, 0., 10.);
      _hist_jetcount = bookHistogram1D("d12-x01-y01", 32, 0, 100);
      _hist_jetpt = bookHistogram1D("d13-x01-y01", 32, 20., 100.);
      _hist_jetlogpt = bookHistogram1D("d14-x01-y01", 32, 0., 20.);
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
    //AIDA::IHistogram1D* _hist_zpthigh;
    //AIDA::IHistogram1D* _hist_zlogpthigh;
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
