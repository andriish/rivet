// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Particle.fhh"

namespace Rivet {


  class ATLAS_2014_I1288706 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1288706()
      : Analysis("ATLAS_2014_I1288706")
    {
     }


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

       // Set up projections
       FinalState fs;

       ZFinder zfinder_ext_dressed_mu(fs, Cuts::abseta<2.4 && Cuts::pT>6.0*GeV, PID::MUON, 12.0*GeV, 66.0*GeV, 0.1);
       declare(zfinder_ext_dressed_mu, "ZFinder_ext_dressed_mu");

       ZFinder zfinder_dressed_mu(fs, Cuts::abseta<2.4 && Cuts::pT>12*GeV, PID::MUON, 26.0*GeV, 66.0*GeV, 0.1);
       declare(zfinder_dressed_mu, "ZFinder_dressed_mu");

       ZFinder zfinder_dressed_el(fs, Cuts::abseta<2.4 && Cuts::pT>12*GeV, PID::ELECTRON, 26.0*GeV, 66.0*GeV, 0.1);
       declare(zfinder_dressed_el, "ZFinder_dressed_el");

       book(_hist_ext_mu_dressed ,1, 1, 1); // muon, dressed level, extended phase space
       book(_hist_mu_dressed	    ,2, 1, 1); // muon, dressed level, nominal phase space
       book(_hist_el_dressed	    ,2, 1, 2); // electron, dressed level, nominal phase space
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ZFinder& zfinder_ext_dressed_mu = apply<ZFinder>(event, "ZFinder_ext_dressed_mu");
      const ZFinder& zfinder_dressed_mu     = apply<ZFinder>(event, "ZFinder_dressed_mu"    );
      const ZFinder& zfinder_dressed_el     = apply<ZFinder>(event, "ZFinder_dressed_el"    );

      FillPlots(zfinder_ext_dressed_mu, _hist_ext_mu_dressed, 9.0);
      FillPlots(zfinder_dressed_mu,     _hist_mu_dressed,    15.0);
      FillPlots(zfinder_dressed_el,     _hist_el_dressed,    15.0);
    }


    void FillPlots(const ZFinder& zfinder, Histo1DPtr hist, double leading_pT) {

      if(zfinder.bosons().size() != 1) return;

      const FourMomentum el1 = zfinder.particles()[0].momentum();
      const FourMomentum el2 = zfinder.particles()[1].momentum();

      if (el1.pT() > leading_pT*GeV || el2.pT() > leading_pT*GeV) {
   	    double mass = zfinder.bosons()[0].mass()/GeV;
   	    hist->fill(mass);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_hist_ext_mu_dressed, crossSection()/sumOfWeights());
      scale(_hist_mu_dressed,     crossSection()/sumOfWeights());
      scale(_hist_el_dressed,     crossSection()/sumOfWeights());

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hist_ext_mu_dressed;
    Histo1DPtr _hist_mu_dressed;
    Histo1DPtr _hist_el_dressed;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1288706);

}
