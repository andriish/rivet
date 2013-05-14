// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetYODA.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/CentralEtHCM.hh"

namespace Rivet {


  /// H1 energy flow in DIS
  ///
  /// @todo Check!
  /// @author Leif Lonnblad
  class H1_1995_S3167097 : public Analysis {
  public:

    /// Constructor
    H1_1995_S3167097()
      : Analysis("H1_1995_S3167097")
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      const DISKinematics& diskin = addProjection(DISKinematics(), "Kinematics");
      const DISFinalState& fshcm = addProjection(DISFinalState(diskin, DISFinalState::HCM), "FS");
      addProjection(CentralEtHCM(fshcm), "Y1HCM");

      const size_t NBINS = 9;
      _hEtFlow = vector<Histo1DPtr>(NBINS);
      _hEtFlowStat = vector<Histo1DPtr>(NBINS);
      _nev = vector<double>(NBINS);
      /// @todo Automate this sort of thing so that the analysis code is more readable.
      for (size_t i = 0; i < NBINS; ++i) {
        string istr(1, char('1' + i));
        _hEtFlow[i] = bookHisto1D(istr, 24, -6, 6);
        _hEtFlowStat[i] = bookHisto1D(istr, 24, -6, 6);
      }
      /// @todo Replace with really temp (non-persistent) histos?
      // _hAvEt = bookScatter2D("21");
      // _hAvX  = bookScatter2D("22");
      // _hAvQ2 = bookScatter2D("23");
      _hAvEt = bookHisto1D("21tmp", NBINS, 1.0, 10.0);
      _hAvX  = bookHisto1D("22tmp", NBINS, 1.0, 10.0);
      _hAvQ2 = bookHisto1D("23tmp", NBINS, 1.0, 10.0);
      _hN    = bookHisto1D("24", NBINS, 1.0, 10.0);
    }


    /// Calculate the bin number from the DISKinematics projection
    size_t _getbin(const DISKinematics& dk) {
      if (inRange(dk.Q2()/GeV2, 5.0, 10.0)) {
        if (inRange(dk.x(), 1e-4, 2e-4)) return 0;
        if (inRange(dk.x(), 2e-4, 5e-4) && dk.Q2() > 6.0*GeV2) return 1;
      } else if (inRange(dk.Q2()/GeV2, 10.0, 20.0)) {
        if (inRange(dk.x(), 2e-4, 5e-4)) return 2;
        if (inRange(dk.x(), 5e-4, 8e-4)) return 3;
        if (inRange(dk.x(), 8e-4, 1.5e-3)) return 4;
        if (inRange(dk.x(), 1.5e-3, 4e-3)) return 5;
      } else if (inRange(dk.Q2()/GeV2, 20.0, 50.0)) {
        if (inRange(dk.x(), 5e-4, 1.4e-3)) return 6;
        if (inRange(dk.x(), 1.4e-3, 3e-3)) return 7;
        if (inRange(dk.x(), 3e-3, 1e-2)) return 8;
      }
      return -1;
    }


    void analyze(const Event& event) {
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      const DISKinematics& dk = applyProjection<DISKinematics>(event, "Kinematics");
      const CentralEtHCM y1 = applyProjection<CentralEtHCM>(event, "Y1HCM");

      const int ibin = _getbin(dk);
      if (ibin < 0) vetoEvent;
      const double weight = event.weight();

      for (size_t i = 0, N = fs.particles().size(); i < N; ++i) {
        const double rap = fs.particles()[i].momentum().rapidity();
        const double et = fs.particles()[i].momentum().Et();
        _hEtFlow[ibin]->fill(rap, weight * et/GeV);
        _hEtFlowStat[ibin]->fill(rap, weight * et/GeV);
      }

      _nev[ibin] += weight;
      _hAvEt->fill(ibin + 1.5, weight * y1.sumEt()/GeV);
      _hAvX->fill(ibin + 1.5, weight * dk.x());
      _hAvQ2->fill(ibin + 1.5, weight * dk.Q2()/GeV2);
      _hN->fill(ibin + 1.5, weight);
    }


    void finalize() {
      for (size_t ibin = 0; ibin < 9; ++ibin) {
        scale(_hEtFlow[ibin], 0.5/_nev[ibin]);
        scale(_hEtFlowStat[ibin], 0.5/_nev[ibin]);
      }

      "/H1_1995_S3167097/21";

      divide(_tmphAvEt, *_hN, );
      h->setTitle(_hAvEt->title());
      histogramFactory().destroy(_hAvEt);

      h = histogramFactory().divide("/H1_1995_S3167097/22", *_hAvX, *_hN);
      h->setTitle(_hAvX->title());
      histogramFactory().destroy(_hAvX);

      h = histogramFactory().divide("/H1_1995_S3167097/23", *_hAvQ2, *_hN);
      h->setTitle(_hAvQ2->title());
      histogramFactory().destroy(_hAvQ2);
    }

    //@}


  private:

    /// Histograms for the \f$ E_T \f$ flows
    vector<Histo1DPtr> _hEtFlow, _hEtFlowStat;

    /// Histograms for averages in different kinematical bins.
    Histo1DPtr _hAvEt, _hAvX, _hAvQ2, _hN;

    /// Helper vector;
    vector<double> _nev;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(H1_1995_S3167097);

}
