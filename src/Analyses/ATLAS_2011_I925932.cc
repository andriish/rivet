// -*- C++ -*-
// ATLAS W pT analysis

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {


  class ATLAS_2011_I925932 : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_I925932() : Analysis("ATLAS_2011_I925932") {  }


    /// @name Analysis methods
    //@{

    void init() {
      // Set up projections
      WFinder wfinder_dressed_el(-2.4, 2.4, 20.0*GeV, ELECTRON, 0.0*GeV, 1000.0*GeV, 25.0*GeV, 0.2,true,false);
      addProjection(wfinder_dressed_el, "WFinder_dressed_el");
      WFinder wfinder_bare_el(-2.4, 2.4, 20.0*GeV, ELECTRON, 0.0*GeV, 1000.0*GeV, 25.0*GeV, 0.0,true,false);
      addProjection(wfinder_bare_el, "WFinder_bare_el");
      WFinder wfinder_dressed_mu(-2.4, 2.4, 20.0*GeV, MUON, 0.0*GeV, 1000.0*GeV, 25.0*GeV, 0.2,true,false);
      addProjection(wfinder_dressed_mu, "WFinder_dressed_mu");
      WFinder wfinder_bare_mu(-2.4, 2.4, 20.0*GeV, MUON, 0.0*GeV, 1000.0*GeV, 25.0*GeV, 0.0,true,false);
      addProjection(wfinder_bare_mu, "WFinder_bare_mu");

      // Counters
      _sumW_e_dressed = 0;
      _sumW_e_bare = 0;
      _sumW_mu_dressed = 0;
      _sumW_mu_bare = 0;

      // Book histograms
      _hist_wpt_dressed_el  = bookHistogram1D(1, 1, 1);
      _hist_wpt_bare_el     = bookHistogram1D(1, 1, 2);
      _hist_wpt_dressed_mu  = bookHistogram1D(2, 1, 1);
      _hist_wpt_bare_mu     = bookHistogram1D(2, 1, 2);
    }


    inline double calcMT(const FourMomentum& a, const FourMomentum& b) {
      return sqrt(2.0 * a.pT() * b.pT() * (1.0 - cos(a.phi() - b.phi())) );
    }


    /// Do the analysis
    void analyze(const Event& event) {

      const WFinder& wfinder_dressed_el = applyProjection<WFinder>(event, "WFinder_dressed_el");
      const WFinder& wfinder_bare_el    = applyProjection<WFinder>(event, "WFinder_bare_el");
      const WFinder& wfinder_dressed_mu = applyProjection<WFinder>(event, "WFinder_dressed_mu");
      const WFinder& wfinder_bare_mu    = applyProjection<WFinder>(event, "WFinder_bare_mu");

      if (wfinder_dressed_el.particles().empty() && wfinder_bare_el.particles().empty() &&
          wfinder_dressed_mu.particles().empty() && wfinder_bare_mu.particles().empty()) {
        MSG_DEBUG("No W bosons found");
        vetoEvent;
      }

      // "Dressed" electron
      if (!wfinder_dressed_el.particles().empty()) {
	    const FourMomentum el = wfinder_dressed_el.constituentLeptons()[0].momentum();
	    const FourMomentum nu = wfinder_dressed_el.constituentNeutrinos()[0].momentum();
	    if (calcMT(el, nu) > 40.0*GeV && nu.pT() > 25.0*GeV) {
          _sumW_e_dressed += event.weight();
          _hist_wpt_dressed_el->fill(wfinder_dressed_el.bosons()[0].momentum().pT()/GeV, event.weight());
	    }
      }

      // "Bare" electron
      if (!wfinder_bare_el.particles().empty()) {
	    const FourMomentum el = wfinder_bare_el.constituentLeptons()[0].momentum();
	    const FourMomentum nu = wfinder_bare_el.constituentNeutrinos()[0].momentum();
	    if (calcMT(el, nu) > 40.0*GeV && nu.pT() > 25.0*GeV) {
          _sumW_e_bare += event.weight();
          _hist_wpt_bare_el->fill(wfinder_bare_el.bosons()[0].momentum().pT()/GeV, event.weight());
	    }
      }

      // "Dressed" muon
      if (!wfinder_dressed_mu.particles().empty()) {
	    const FourMomentum mu = wfinder_dressed_mu.constituentLeptons()[0].momentum();
	    const FourMomentum nu = wfinder_dressed_mu.constituentNeutrinos()[0].momentum();
	    if (calcMT(mu, nu) > 40.0*GeV && nu.pT() > 25.0*GeV) {
          _sumW_mu_dressed += event.weight();
          _hist_wpt_dressed_mu->fill(wfinder_dressed_mu.bosons()[0].momentum().pT()/GeV, event.weight());
	    }
      }

      // "Bare" muon
      if (!wfinder_bare_mu.particles().empty()) {
	    const FourMomentum mu = wfinder_bare_mu.constituentLeptons()[0].momentum();
	    const FourMomentum nu = wfinder_bare_mu.constituentNeutrinos()[0].momentum();
	    if (calcMT(mu, nu) > 40.0*GeV && nu.pT() > 25.0*GeV) {
          _sumW_mu_bare += event.weight();
          _hist_wpt_bare_mu->fill(wfinder_bare_mu.bosons()[0].momentum().pT()/GeV, event.weight());
	    }
      }

    }


    // Normalize histos
    void finalize() {
      scale(_hist_wpt_dressed_el, 1/_sumW_e_dressed);
      scale(_hist_wpt_bare_el, 1/_sumW_e_bare);
      scale(_hist_wpt_dressed_mu, 1/_sumW_mu_dressed);
      scale(_hist_wpt_bare_mu, 1/_sumW_mu_bare);
    }

    //@}


  private:

	double _sumW_e_dressed, _sumW_e_bare, _sumW_mu_dressed, _sumW_mu_bare;

	AIDA::IHistogram1D* _hist_wpt_dressed_el;
	AIDA::IHistogram1D* _hist_wpt_bare_el;
	AIDA::IHistogram1D* _hist_wpt_dressed_mu;
	AIDA::IHistogram1D* _hist_wpt_bare_mu;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I925932);

}
