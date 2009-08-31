// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class MC_LHC_ZANALYSIS : public Analysis {
  public:
  
  /// Default constructor
    MC_LHC_ZANALYSIS()
      : Analysis("MC_LHC_ZANALYSIS") {
      const ChargedFinalState cfs;
      addProjection(cfs, "CFS");
      const ZFinder zf(-MaxRapidity, MaxRapidity, 0.0*GeV, ELECTRON, 30.0*GeV, 115.0*GeV, 0.2);
      addProjection(zf, "ZF");
      FastJets fastjets(zf.remainingFinalState(), FastJets::KT, 0.7);
      addProjection(fastjets, "Jets");
    }
    

    /// @name Analysis methods
    /// @todo change "Weights" to differential cross sections once histos normalised to cross-section.
    //@{

    void init() { 
      _hist_chargemulti = bookHistogram1D("d01-x01-y01", 30, 0.5, 250.5);
      _hist_chargept = bookHistogram1D("d02-x01-y01",  32, 0., 25.);
      _hist_chargemeanpt = bookHistogram1D("d03-x01-y01", 25, 0., 10.);
      _hist_chargermspt = bookHistogram1D("d04-x01-y01", 32, 0., 10.);
      _hist_zcount = bookHistogram1D("d05-x01-y01", 30, 0., 15.);
      _hist_zpt = bookHistogram1D("d06-x01-y01", 32, 0., 25.);
      _hist_zlogpt = bookHistogram1D("d07-x01-y01", 32, 0., 6.);
      _hist_zeta = bookHistogram1D("d08-x01-y01", 32, -6., 6.);
      _hist_zphi = bookHistogram1D("d09-x01-y01", 32, 0., 6.4);
      _hist_zmass = bookHistogram1D("d10-x01-y01", 40, 60., 100.);
      _hist_zlogmass = bookHistogram1D("d11-x01-y01",32, 0., 10.);
      _hist_jetcount = bookHistogram1D("d12-x01-y01", 32, 0, 100);
      _hist_jetpt = bookHistogram1D("d13-x01-y01",32, 20., 100.);
      _hist_jetlogpt = bookHistogram1D("d14-x01-y01", 32, 0., 20.);
    }
    
    
    void analyze(const Event& event) {
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      const ZFinder zf = applyProjection<ZFinder>(event, "ZF");
      const FastJets fastjets = applyProjection<FastJets>(event, "Jets");
      const Jets jets = fastjets.jetsByPt();
      const double weight = event.weight();
      double meanpt(0), rmspt(0);
    
      // Charged particles part    
      _hist_chargemulti->fill(cfs.particles().size(), weight);
      foreach(Particle p, cfs.particles()) {
        _hist_chargept->fill(p.momentum().pT(), weight);
        meanpt = meanpt + p.momentum().pT();
        rmspt = rmspt + p.momentum().pT()*p.momentum().pT();
      }
      
      meanpt = meanpt/ cfs.particles().size();
      _hist_chargemeanpt->fill(meanpt, weight);
      rmspt = sqrt(rmspt / cfs.particles().size());
      _hist_chargermspt->fill(rmspt, weight);
      
      // Z part
      _hist_zcount->fill(zf.particles().size(), weight);
      foreach (Particle zp, zf.particles()) {
        _hist_zpt->fill(zp.momentum().pT(), weight);
        _hist_zlogpt->fill(log(zp.momentum().pT()), weight);
        _hist_zeta->fill(zp.momentum().pseudorapidity(), weight);
        _hist_zphi->fill(zp.momentum().azimuthalAngle(), weight);
        _hist_zmass->fill(zp.momentum().mass(), weight);
        _hist_zlogmass->fill(log(zp.momentum().mass()), weight);	
      }
      
      // Jet part
      _hist_jetcount->fill(fastjets.size(), weight);
      foreach(Jet j, fastjets.jetsByPt()) {
        _hist_jetpt->fill(j.momentum().pT(), weight);
        _hist_jetlogpt->fill(log(j.momentum().pT()), weight);
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
    AIDA::IHistogram1D* _hist_zcount;
    AIDA::IHistogram1D* _hist_zpt;
    AIDA::IHistogram1D* _hist_zlogpt;
    //AIDA::IHistogram1D* _hist_zpthigh;
    //AIDA::IHistogram1D* _hist_zlogpthigh;
    AIDA::IHistogram1D* _hist_zeta;
    AIDA::IHistogram1D* _hist_zphi;
    AIDA::IHistogram1D* _hist_zmass;
    AIDA::IHistogram1D* _hist_zlogmass;
    AIDA::IHistogram1D* _hist_jetcount;
    AIDA::IHistogram1D* _hist_jetpt;
    AIDA::IHistogram1D* _hist_jetlogpt;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_LHC_ZANALYSIS> plugin_MC_LHC_ZANALYSIS;

}
