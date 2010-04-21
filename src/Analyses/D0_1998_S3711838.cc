// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  /// @brief D0 Run I Z boson \f$ p_\perp \f$ measurement
  class D0_1998_S3711838 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    D0_1998_S3711838()
      : Analysis("D0_1998_S3711838")
    {
      /// @todo Set approriate for your analysis
      setBeams(PROTON, ANTIPROTON);

      /// @todo Set whether your finalize method needs the generator cross section
      setNeedsCrossSection(false);
    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      /// @todo Use separate pT and ETmiss cuts in WFinder
      WFinder wfe(-5, 5, 0.0*GeV, ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      addProjection(wfe, "WFe");

      // Book histogram
      _hist_wpt = bookHistogram1D(1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const WFinder& wf = applyProjection<WFinder>(event, "WFe");
      if (wf.size() == 0) {
        getLog() << Log::DEBUG << "No W candidates found: vetoing" << endl;
        vetoEvent;
      }

      // Require the electron to have ET > 25 GeV  and |eta| < 1.1
      /// @todo Use separate pT and ETmiss cuts in WFinder
      FourMomentum p_e;
      int chg_e = 0;

      foreach (const Particle& l, wf.constituentsFinalState().particles()) {
        const FourMomentum pl = l.momentum();
        if (abs(l.pdgId()) == ELECTRON) {
          chg_e = PID::threeCharge(l.pdgId());
          p_e = pl;
          const double eta_e = fabs(p_e.pseudorapidity());
          if ( (pl.Et()/GeV < 25.0) || (eta_e > 1.1) ) {
            getLog() << Log::DEBUG << l.pdgId() << " ET,pT,eta: " << pl.Et()/GeV << "," << pl.pT()/GeV << "," << eta_e << " fails electron cut" << endl;
            vetoEvent;
          }
        }
        // Require ETmiss to be > 25 GeV
        else if (abs(l.pdgId()) == NU_E) {
          FourMomentum p_nu = l.momentum();
          if (p_nu.Et()/GeV < 25.0) {
            getLog() << Log::DEBUG << l.pdgId() << " ET(miss): " << p_nu.Et() << "fails ETmiss cut" << endl;
            vetoEvent;
          }
        }
      }
      assert(chg_e != 0);

      FourMomentum pW = wf.particles()[0].momentum();
      getLog() << Log::DEBUG << "Dilepton mass = " << pW.mass()/GeV << " GeV"  << endl;
      getLog() << Log::DEBUG << "Dilepton pT   = " << pW.pT()/GeV << " GeV" << endl;
      _hist_wpt->fill(pW.pT()/GeV, weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

    }

    //@}

  private:

    /// @name Histograms
    AIDA::IHistogram1D *_hist_wpt;

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<D0_1998_S3711838> plugin_D0_1998_S3711838;


}
