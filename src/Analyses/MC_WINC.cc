// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  /// @brief MC validation analysis for inclusive W events
  class MC_WINC : public Analysis {
  public:

    /// Default constructor
    MC_WINC()
      : Analysis("MC_WINC")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      WFinder wfinder(fs, -3.5, 3.5, 25.0*GeV, ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      addProjection(wfinder, "WFinder");

      _h_W_mass = bookHistogram1D("W_mass", 50, 55.0, 105.0);
      _h_W_pT = bookHistogram1D("W_pT", logBinEdges(100, 1.0, 0.5*sqrtS()));
      _h_W_pT_peak = bookHistogram1D("W_pT_peak", 25, 0.0, 125.0);
      _h_W_y = bookHistogram1D("W_y", 40, -4.0, 4.0);
      _h_W_phi = bookHistogram1D("W_phi", 25, 0.0, TWOPI);
      _h_Wplus_pT = bookHistogram1D("Wplus_pT", logBinEdges(25, 10.0, 0.5*sqrtS()));
      _h_Wminus_pT = bookHistogram1D("Wminus_pT", logBinEdges(25, 10.0, 0.5*sqrtS()));
      _h_lepton_pT = bookHistogram1D("lepton_pT", logBinEdges(100, 10.0, 0.25*sqrtS()));
      _h_lepton_eta = bookHistogram1D("lepton_eta", 40, -4.0, 4.0);
      _htmp_dsigminus_deta = bookHistogram1D("lepton_dsigminus_deta", 20, 0.0, 4.0);
      _htmp_dsigplus_deta  = bookHistogram1D("lepton_dsigplus_deta", 20, 0.0, 4.0);
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = e.weight();

      const WFinder& wfinder = applyProjection<WFinder>(e, "WFinder");
      if (wfinder.bosons().size() != 1) {
        vetoEvent;
      }

      double charge3_x_eta = 0;
      int charge3 = 0;
      FourMomentum emom;
      FourMomentum wmom(wfinder.bosons().front().momentum());
      _h_W_mass->fill(wmom.mass(), weight);
      _h_W_pT->fill(wmom.pT(), weight);
      _h_W_pT_peak->fill(wmom.pT(), weight);
      _h_W_y->fill(wmom.rapidity(), weight);
      _h_W_phi->fill(wmom.azimuthalAngle(), weight);
      Particle l=wfinder.constituentLeptons()[0];
      _h_lepton_pT->fill(l.momentum().pT(), weight);
      _h_lepton_eta->fill(l.momentum().eta(), weight);
      if (PID::threeCharge(l.pdgId()) != 0) {
        emom = l.momentum();
        charge3_x_eta = PID::threeCharge(l.pdgId()) * emom.eta();
        charge3 = PID::threeCharge(l.pdgId());
      }
      assert(charge3_x_eta != 0);
      assert(charge3!=0);
      if (emom.Et() > 30/GeV) {
        if (charge3_x_eta < 0) {
          _htmp_dsigminus_deta->fill(emom.eta(), weight);
        } else {
          _htmp_dsigplus_deta->fill(emom.eta(), weight);
        }
      }
      if (charge3 < 0) {
        _h_Wminus_pT->fill(wmom.pT(), weight);
      } else {
        _h_Wplus_pT->fill(wmom.pT(), weight);
      }
    }


    /// Finalize
    void finalize() {
      scale(_h_W_mass, crossSection()/sumOfWeights());
      scale(_h_W_pT, crossSection()/sumOfWeights());
      scale(_h_W_pT_peak, crossSection()/sumOfWeights());
      scale(_h_W_y, crossSection()/sumOfWeights());
      scale(_h_W_phi, crossSection()/sumOfWeights());
      scale(_h_lepton_pT, crossSection()/sumOfWeights());
      scale(_h_lepton_eta, crossSection()/sumOfWeights());

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta) for each Et region
      AIDA::IHistogramFactory& hf = histogramFactory();
      IHistogram1D* numtmp = hf.subtract("/numtmp", *_htmp_dsigplus_deta, *_htmp_dsigminus_deta);
      IHistogram1D* dentmp = hf.add("/dentmp", *_htmp_dsigplus_deta, *_htmp_dsigminus_deta);
      assert(numtmp && dentmp);
      hf.divide(histoDir() + "/W_chargeasymm_eta", *numtmp, *dentmp);
      hf.destroy(numtmp);
      hf.destroy(dentmp);
      hf.destroy(_htmp_dsigminus_deta);
      hf.destroy(_htmp_dsigplus_deta);

      // W charge asymmetry vs. pTW: dsig+/dpT / dsig-/dpT
      hf.divide(histoDir() + "/W_chargeasymm_pT", *_h_Wplus_pT, *_h_Wminus_pT);
      scale(_h_Wplus_pT, crossSection()/sumOfWeights());
      scale(_h_Wminus_pT, crossSection()/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D * _h_W_mass;
    AIDA::IHistogram1D * _h_W_pT;
    AIDA::IHistogram1D * _h_W_pT_peak;
    AIDA::IHistogram1D * _h_W_y;
    AIDA::IHistogram1D * _h_W_phi;
    AIDA::IHistogram1D * _h_Wplus_pT;
    AIDA::IHistogram1D * _h_Wminus_pT;
    AIDA::IHistogram1D * _h_lepton_pT;
    AIDA::IHistogram1D * _h_lepton_eta;

    AIDA::IHistogram1D * _htmp_dsigminus_deta;
    AIDA::IHistogram1D * _htmp_dsigplus_deta;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_WINC);

}
