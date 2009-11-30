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
      setNeedsCrossSection(true);
      _sumWPass = 0.0;
    }
 

    /// @name Analysis methods
    /// @todo change "Weights" to differential cross sections once histos normalised to cross-section.
    //@{

    void init() {
      const ChargedFinalState cfs;
      addProjection(cfs, "CFS");
      /// @todo Handle muon-decay Ws as well
      const WFinder wf(-MAXRAPIDITY, MAXRAPIDITY, 0.0*GeV, ELECTRON, 70.0*GeV, 90.0*GeV, 0.2);
      addProjection(wf, "WF");
      FastJets fastjets(wf.remainingFinalState(), FastJets::KT, 0.7);
      addProjection(fastjets, "Jets");

      _hist_chargemulti = bookHistogram1D("track-n", 50, -0.5, 399.5);
      _hist_chargept = bookHistogram1D("track-pt", 25, 0, 25);
      /// @todo Use profile plots instead:
      _hist_chargemeanpt = bookHistogram1D("track-ptavg", 20, 0, 10);
      /// @todo Use profile plots instead (and this isn't really an "RMS")
      _hist_chargermspt = bookHistogram1D("track-ptrms", 20, 0, 10);
      //
      _hist_wcount = bookHistogram1D("w-n", 6, -0.5, 5.5);
      _hist_wpluscount = bookHistogram1D("wplus-n", 6, -0.5, 5.5);
      _hist_wminuscount = bookHistogram1D("wminus-n", 6, -0.5, 5.5);
      _hist_wpt = bookHistogram1D("w-pt", 25, 0, 100);
      _hist_wpluspt = bookHistogram1D("wplus-pt", 25, 0, 100);
      _hist_wminuspt = bookHistogram1D("wminus-pt", 25, 0, 100);
      _hist_weta = bookHistogram1D("w-eta", 36, -6, 6);
      _hist_wpluseta = bookHistogram1D("wplus-eta", 36, -6, 6);
      _hist_wminuseta = bookHistogram1D("wminus-eta", 36, -6, 6);
      _hist_wphi = bookHistogram1D("w-phi", 25, 0, TWOPI);
      _hist_wplusphi = bookHistogram1D("wplus-phi", 25, 0, TWOPI);
      _hist_wminusphi = bookHistogram1D("wminus-phi", 25, 0, TWOPI);
      _hist_wmass = bookHistogram1D("w-m", 40, 60, 100);
      _hist_wplusmass = bookHistogram1D("wplus-m", 40, 60, 100);
      _hist_wminusmass = bookHistogram1D("wminus-m", 40, 60, 100);
      //_hist_weta_asymm = bookProfile1D("asymm-eta-w", 20, -5.0, 5.0);
      //
      _hist_jetcount = bookHistogram1D("jet-n", 6, -0.5, 5.5);
      _hist_jetpt = bookHistogram1D("jet-pt", 50, 20, 100);
    }
 
 
    void analyze(const Event& event) {
      const WFinder& wf = applyProjection<WFinder>(event, "WF");
      if (wf.size() == 0) {
        getLog() << Log::DEBUG << "No W candidates found: vetoing" << endl;
        vetoEvent;
      }
      const double weight = event.weight();
      _sumWPass += weight;
 
      // Charged particles part
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
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
      uint n_wplus(0), n_wminus(0);
      foreach (const Particle& wp, wf.particles()) {
        const double pT = wp.momentum().pT();
        const double eta = wp.momentum().eta();
        const double phi = wp.momentum().phi();
        const double m = wp.momentum().mass();
        /// @todo When histo handling is easier, build total histos by summing W+/- histos
        _hist_wpt->fill(pT/GeV, weight);
        _hist_weta->fill(eta, weight);
        _hist_wphi->fill(phi, weight);
        _hist_wmass->fill(m/GeV, weight);
        if (wp.pdgId() == WPLUSBOSON) {
          n_wplus += 1;
          _hist_wpluspt->fill(pT/GeV, weight);
          _hist_wpluseta->fill(eta, weight);
          _hist_wplusphi->fill(phi, weight);
          _hist_wplusmass->fill(m/GeV, weight);
        } else if (wp.pdgId() == WMINUSBOSON) {
          n_wminus += 1;
          _hist_wminuspt->fill(pT/GeV, weight);
          _hist_wminuseta->fill(eta, weight);
          _hist_wminusphi->fill(phi, weight);
          _hist_wminusmass->fill(m/GeV, weight);
        } else {
          // Just checking!
          throw Error("There shouldn't be any W candidates without a W PID!");
        }
      }
      _hist_wcount->fill(n_wplus+n_wminus, weight);
      _hist_wpluscount->fill(n_wplus, weight);   
      _hist_wminuscount->fill(n_wminus, weight);

      // Jet part
      const FastJets& fastjets = applyProjection<FastJets>(event, "Jets");
      const Jets jets = fastjets.jetsByPt();
      _hist_jetcount->fill(jets.size(), weight);
      foreach (const Jet& j, jets) {
        const double pT = j.momentum().pT();
        _hist_jetpt->fill(pT/GeV, weight);
      }

    }
 
 
    void finalize() {
      const double xsec_sumw = crossSectionPerEvent()/picobarn;
      /// @todo Actually "measure" separate W+ and W- xsecs 
      const double xsec_sumw_plus = xsec_sumw/2.0;
      const double xsec_sumw_minus = xsec_sumw/2.0;

      scale(_hist_chargemulti, xsec_sumw);
      scale(_hist_chargept, xsec_sumw);
      scale(_hist_chargemeanpt, xsec_sumw);
      scale(_hist_chargermspt, xsec_sumw);
      scale(_hist_wcount, xsec_sumw);
      scale(_hist_wpluscount, xsec_sumw_plus);
      scale(_hist_wminuscount, xsec_sumw_minus);
      scale(_hist_wpt, xsec_sumw);
      scale(_hist_wpluspt, xsec_sumw_plus);
      scale(_hist_wminuspt, xsec_sumw_minus);
      scale(_hist_weta, xsec_sumw);
      scale(_hist_wpluseta, xsec_sumw_plus);
      scale(_hist_wminuseta, xsec_sumw_minus);
      scale(_hist_wphi, xsec_sumw);
      scale(_hist_wplusphi, xsec_sumw_plus);
      scale(_hist_wminusphi, xsec_sumw_minus);
      scale(_hist_wmass, xsec_sumw);
      scale(_hist_wplusmass, xsec_sumw_plus);
      scale(_hist_wminusmass, xsec_sumw_minus);
      scale(_hist_jetcount, xsec_sumw);
      scale(_hist_jetpt, xsec_sumw);
    }
 
    //@}

 
  private:


    /// @name Counters
    //@{
    double _sumWPass;
    //@}


    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _hist_chargemulti;
    AIDA::IHistogram1D* _hist_chargept;
    AIDA::IHistogram1D* _hist_chargemeanpt;
    AIDA::IHistogram1D* _hist_chargermspt;
    AIDA::IHistogram1D *_hist_wcount, *_hist_wpluscount, *_hist_wminuscount;
    AIDA::IHistogram1D *_hist_wpt, *_hist_wpluspt, *_hist_wminuspt;
    AIDA::IHistogram1D *_hist_weta, *_hist_wpluseta, *_hist_wminuseta;
    AIDA::IHistogram1D *_hist_wphi, *_hist_wplusphi, *_hist_wminusphi;
    AIDA::IHistogram1D *_hist_wmass, *_hist_wplusmass, *_hist_wminusmass;
    AIDA::IHistogram1D* _hist_jetcount;
    AIDA::IHistogram1D* _hist_jetpt;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_LHC_WANALYSIS> plugin_MC_LHC_WANALYSIS;

}
