// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetYODA.hh"

namespace Rivet {


  /// @brief MC validation analysis for higgs [-> tau tau] events
  class MC_HINC : public Analysis {
  public:

    /// Default constructor
    MC_HINC()
      : Analysis("MC_HINC")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      ZFinder hfinder(fs, -3.5, 3.5, 25.0*GeV, TAU, 115.0*GeV, 125.0*GeV, 0.0, false, false);
      addProjection(hfinder, "Hfinder");

      _h_H_mass = bookHisto1D("H_mass", 50, 119.7, 120.3);
      _h_H_pT = bookHisto1D("H_pT", logspace(100, 1.0, 0.5*sqrtS()));
      _h_H_pT_peak = bookHisto1D("H_pT_peak", 25, 0.0, 25.0);
      _h_H_y = bookHisto1D("H_y", 40, -4.0, 4.0);
      _h_H_phi = bookHisto1D("H_phi", 25, 0.0, TWOPI);
      _h_lepton_pT = bookHisto1D("lepton_pT", logspace(100, 10.0, 0.25*sqrtS()));
      _h_lepton_eta = bookHisto1D("lepton_eta", 40, -4.0, 4.0);

    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& hfinder = applyProjection<ZFinder>(e, "Hfinder");
      if (hfinder.bosons().size()!=1) {
        vetoEvent;
      }
      const double weight = e.weight();

      FourMomentum hmom(hfinder.bosons()[0].momentum());
      _h_H_mass->fill(hmom.mass(),weight);
      _h_H_pT->fill(hmom.pT(),weight);
      _h_H_pT_peak->fill(hmom.pT(),weight);
      _h_H_y->fill(hmom.rapidity(),weight);
      _h_H_phi->fill(hmom.azimuthalAngle(),weight);
      foreach (const Particle& l, hfinder.constituents()) {
        _h_lepton_pT->fill(l.momentum().pT(), weight);
        _h_lepton_eta->fill(l.momentum().eta(), weight);
      }
    }


    /// Finalize
    void finalize() {
      scale(_h_H_mass, crossSection()/sumOfWeights());
      scale(_h_H_pT, crossSection()/sumOfWeights());
      scale(_h_H_pT_peak, crossSection()/sumOfWeights());
      scale(_h_H_y, crossSection()/sumOfWeights());
      scale(_h_H_phi, crossSection()/sumOfWeights());
      scale(_h_lepton_pT, crossSection()/sumOfWeights());
      scale(_h_lepton_eta, crossSection()/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_H_mass;
    Histo1DPtr _h_H_pT;
    Histo1DPtr _h_H_pT_peak;
    Histo1DPtr _h_H_y;
    Histo1DPtr _h_H_phi;
    Histo1DPtr _h_lepton_pT;
    Histo1DPtr _h_lepton_eta;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_HINC);

}
