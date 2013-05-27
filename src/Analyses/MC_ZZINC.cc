// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {

  /// @brief MC validation analysis for Z[ee]Z[mumu] events
  class MC_ZZINC : public Analysis {
  public:

    /// Default constructor
    MC_ZZINC()
      : Analysis("MC_ZZINC")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      ZFinder zeefinder(fs, -3.5, 3.5, 25.0*GeV, PID::ELECTRON, 65.0*GeV, 115.0*GeV, 0.2, true, true);
      addProjection(zeefinder, "ZeeFinder");

      VetoedFinalState zmminput;
      zmminput.addVetoOnThisFinalState(zeefinder);
      ZFinder zmmfinder(zmminput, -3.5, 3.5, 25.0*GeV, PID::MUON, 65.0*GeV, 115.0*GeV, 0.2, true, true);
      addProjection(zmmfinder, "ZmmFinder");

      // properties of the pair momentum
      _h_ZZ_pT = bookHisto1D("ZZ_pT", logspace(100, 1.0, 0.5*sqrtS()));
      _h_ZZ_pT_peak = bookHisto1D("ZZ_pT_peak", 25, 0.0, 25.0);
      _h_ZZ_eta = bookHisto1D("ZZ_eta", 40, -7.0, 7.0);
      _h_ZZ_phi = bookHisto1D("ZZ_phi", 25, 0.0, TWOPI);
      _h_ZZ_m = bookHisto1D("ZZ_m", logspace(100, 150.0, 180.0+0.25*sqrtS()));

      // correlations between the ZZ
      _h_ZZ_dphi = bookHisto1D("ZZ_dphi", 25, 0.0, PI);  /// @todo non-linear?
      _h_ZZ_deta = bookHisto1D("ZZ_deta", 25, -7.0, 7.0);
      _h_ZZ_dR = bookHisto1D("ZZ_dR", 25, 0.5, 7.0);
      _h_ZZ_dpT = bookHisto1D("ZZ_dpT", logspace(100, 1.0, 0.5*sqrtS()));
      _h_ZZ_costheta_planes = bookHisto1D("ZZ_costheta_planes", 25, -1.0, 1.0);

      // properties of the Z bosons
      _h_Z_pT = bookHisto1D("Z_pT", logspace(100, 10.0, 0.25*sqrtS()));
      _h_Z_eta = bookHisto1D("Z_eta", 70, -7.0, 7.0);

      // properties of the leptons
      _h_Zl_pT = bookHisto1D("Zl_pT", logspace(100, 30.0, 0.1
                                                      *sqrtS()));
      _h_Zl_eta = bookHisto1D("Zl_eta", 40, -3.5, 3.5);

      // correlations between the opposite charge leptons
      _h_ZeZm_dphi = bookHisto1D("ZeZm_dphi", 25, 0.0, PI);
      _h_ZeZm_deta = bookHisto1D("ZeZm_deta", 25, -5.0, 5.0);
      _h_ZeZm_dR = bookHisto1D("ZeZm_dR", 25, 0.5, 5.0);
      _h_ZeZm_m = bookHisto1D("ZeZm_m", 100, 0.0, 300.0);

    }



    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = e.weight();

      const ZFinder& zeefinder = applyProjection<ZFinder>(e, "ZeeFinder");
      if (zeefinder.bosons().size()!=1) {
        vetoEvent;
      }

      const ZFinder& zmmfinder = applyProjection<ZFinder>(e, "ZmmFinder");
      if (zmmfinder.bosons().size()!=1) {
        vetoEvent;
      }

      FourMomentum zee(zeefinder.bosons()[0].momentum());
      FourMomentum zmm(zmmfinder.bosons()[0].momentum());
      FourMomentum zz(zee+zmm);
      // find leptons
      FourMomentum ep(zeefinder.constituents()[0].momentum()),
        em(zeefinder.constituents()[1].momentum()),
        mp(zmmfinder.constituents()[0].momentum()),
        mm(zmmfinder.constituents()[1].momentum());


      _h_ZZ_pT->fill(zz.pT(),weight);
      _h_ZZ_pT_peak->fill(zz.pT(),weight);
      _h_ZZ_eta->fill(zz.eta(),weight);
      _h_ZZ_phi->fill(zz.azimuthalAngle(),weight);
      double mzz2=zz.mass2();
      if (mzz2>0.0) _h_ZZ_m->fill(sqrt(mzz2), weight);

      _h_ZZ_dphi->fill(mapAngle0ToPi(zee.phi()-zmm.phi()), weight);
      _h_ZZ_deta->fill(zee.eta()-zmm.eta(), weight);
      _h_ZZ_dR->fill(deltaR(zee,zmm), weight);
      _h_ZZ_dpT->fill(fabs(zee.pT()-zmm.pT()), weight);

      Vector3 crossZee = ep.vector3().cross(em.vector3());
      Vector3 crossZmm = mp.vector3().cross(mm.vector3());
      double costheta = crossZee.dot(crossZmm)/crossZee.mod()/crossZmm.mod();
      _h_ZZ_costheta_planes->fill(costheta, weight);

      _h_Z_pT->fill(zee.pT(),weight);
      _h_Z_pT->fill(zmm.pT(),weight);
      _h_Z_eta->fill(zee.eta(),weight);
      _h_Z_eta->fill(zmm.eta(),weight);

      _h_Zl_pT->fill(ep.pT(), weight);
      _h_Zl_pT->fill(em.pT(), weight);
      _h_Zl_pT->fill(mp.pT(), weight);
      _h_Zl_pT->fill(mm.pT(), weight);
      _h_Zl_eta->fill(ep.eta(), weight);
      _h_Zl_eta->fill(em.eta(), weight);
      _h_Zl_eta->fill(mp.eta(), weight);
      _h_Zl_eta->fill(mm.eta(), weight);

      _h_ZeZm_dphi->fill(mapAngle0ToPi(ep.phi()-mm.phi()), weight);
      _h_ZeZm_deta->fill(ep.eta()-mm.eta(), weight);
      _h_ZeZm_dR->fill(deltaR(ep,mm), weight);
      double m2=FourMomentum(ep+mm).mass2();
      if (m2 < 0) m2 = 0.0;
      _h_ZeZm_m->fill(sqrt(m2), weight);
    }


    /// Finalize
    void finalize() {
      double norm=crossSection()/sumOfWeights();
      scale(_h_ZZ_pT, norm);
      scale(_h_ZZ_pT_peak, norm);
      scale(_h_ZZ_eta, norm);
      scale(_h_ZZ_phi, norm);
      scale(_h_ZZ_m, norm);
      scale(_h_ZZ_dphi, norm);
      scale(_h_ZZ_deta, norm);
      scale(_h_ZZ_dR, norm);
      scale(_h_ZZ_dpT, norm);
      scale(_h_ZZ_costheta_planes, norm);
      scale(_h_Z_pT, norm);
      scale(_h_Z_eta, norm);
      scale(_h_Zl_pT, norm);
      scale(_h_Zl_eta, norm);
      scale(_h_ZeZm_dphi, norm);
      scale(_h_ZeZm_deta, norm);
      scale(_h_ZeZm_dR, norm);
      scale(_h_ZeZm_m, norm);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_ZZ_pT;
    Histo1DPtr _h_ZZ_pT_peak;
    Histo1DPtr _h_ZZ_eta;
    Histo1DPtr _h_ZZ_phi;
    Histo1DPtr _h_ZZ_m;
    Histo1DPtr _h_ZZ_dphi;
    Histo1DPtr _h_ZZ_deta;
    Histo1DPtr _h_ZZ_dR;
    Histo1DPtr _h_ZZ_dpT;
    Histo1DPtr _h_ZZ_costheta_planes;
    Histo1DPtr _h_Z_pT;
    Histo1DPtr _h_Z_eta;
    Histo1DPtr _h_Zl_pT;
    Histo1DPtr _h_Zl_eta;
    Histo1DPtr _h_ZeZm_dphi;
    Histo1DPtr _h_ZeZm_deta;
    Histo1DPtr _h_ZeZm_dR;
    Histo1DPtr _h_ZeZm_m;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ZZINC);

}
