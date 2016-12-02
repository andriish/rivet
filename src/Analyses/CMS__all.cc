#line 1 "CMS_2010_S8547297.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class CMS_2010_S8547297 : public Analysis {
  public:

    CMS_2010_S8547297() : Analysis("CMS_2010_S8547297") {}


    void init() {
      ChargedFinalState cfs(-2.5, 2.5, 0.0*GeV);
      declare(cfs, "CFS");

      if (fuzzyEquals(sqrtS()/GeV, 900)) {
        for (int d=1; d<=3; d++) {
          for (int y=1; y<=4; y++) {
            _h_dNch_dpT.push_back(bookHisto1D(d, 1, y));
          }
        }
        _h_dNch_dpT_all = bookHisto1D(7, 1, 1);
        _h_dNch_dEta = bookHisto1D(8, 1, 1);
      } else if (fuzzyEquals(sqrtS()/GeV, 2360)) {
        for (int d=4; d<=6; d++) {
          for (int y=1; y<=4; y++) {
            _h_dNch_dpT.push_back(bookHisto1D(d, 1, y));
          }
        }
        _h_dNch_dpT_all = bookHisto1D(7, 1, 2);
        _h_dNch_dEta = bookHisto1D(8, 1, 2);
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      //charged particles
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");

      foreach (const Particle& p, charged.particles()) {
        //selecting only charged hadrons
        if (! PID::isHadron(p.pid())) continue;

        const double pT = p.pT();
        const double eta = p.eta();

        // The data is actually a duplicated folded distribution.  This should mimic it.
        _h_dNch_dEta->fill(eta, 0.5*weight);
        _h_dNch_dEta->fill(-eta, 0.5*weight);
        if (fabs(eta) < 2.4 && pT > 0.1*GeV) {
          if (pT < 4.0*GeV) {
            _h_dNch_dpT_all->fill(pT/GeV, weight/(pT/GeV));
            if (pT < 2.0*GeV) {
              int ietabin = int(fabs(eta)/0.2);
              _h_dNch_dpT[ietabin]->fill(pT/GeV, weight);
            }
          }
        }
      }
    }


    void finalize() {
      const double normfac = 1.0/sumOfWeights(); // Normalizing to unit eta is automatic
      // The pT distributions in bins of eta must be normalized to unit eta.  This is a factor of 2
      // for the |eta| times 0.2 (eta range).
      // The pT distributions over all eta are normalized to unit eta (2.0*2.4) and by 1/2*pi*pT.
      // The 1/pT part is taken care of in the filling.  The 1/2pi is taken care of here.
      const double normpT = normfac/(2.0*0.2);
      const double normpTall = normfac/(2.0*M_PI*2.0*2.4);

      for (size_t ietabin=0; ietabin < _h_dNch_dpT.size(); ietabin++){
        scale(_h_dNch_dpT[ietabin], normpT);
      }
      scale(_h_dNch_dpT_all, normpTall);
      scale(_h_dNch_dEta, normfac);
    }


  private:

    std::vector<Histo1DPtr> _h_dNch_dpT;
    Histo1DPtr _h_dNch_dpT_all;
    Histo1DPtr _h_dNch_dEta;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2010_S8547297);

}
#line 1 "CMS_2010_S8656010.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class CMS_2010_S8656010 : public Analysis {
  public:

    CMS_2010_S8656010() : Analysis("CMS_2010_S8656010") {}


    void init() {
      ChargedFinalState cfs(-2.5, 2.5, 0.0*GeV);
      declare(cfs, "CFS");

      for (int d=1; d<=3; d++) {
        for (int y=1; y<=4; y++) {
          _h_dNch_dpT.push_back(bookHisto1D(d, 1, y));
        }
      }

      _h_dNch_dpT_all = bookHisto1D(4, 1, 1);
      _h_dNch_dEta = bookHisto1D(5, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      //charged particles
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");

      foreach (const Particle& p, charged.particles()) {
        //selecting only charged hadrons
        if (! PID::isHadron(p.pid())) continue;

        const double pT = p.pT();
        const double eta = p.eta();

        // The data is actually a duplicated folded distribution.  This should mimic it.
        _h_dNch_dEta->fill(eta, 0.5*weight);
        _h_dNch_dEta->fill(-eta, 0.5*weight);
        if (fabs(eta) < 2.4 && pT > 0.1*GeV) {
          if (pT < 6.0*GeV) {
            _h_dNch_dpT_all->fill(pT/GeV, weight/(pT/GeV));
            if (pT < 2.0*GeV) {
              int ietabin = int(fabs(eta)/0.2);
              _h_dNch_dpT[ietabin]->fill(pT/GeV, weight);
            }
          }
        }
      }
    }


    void finalize() {
      const double normfac = 1.0/sumOfWeights(); // Normalizing to unit eta is automatic
      // The pT distributions in bins of eta must be normalized to unit eta.  This is a factor of 2
      // for the |eta| times 0.2 (eta range).
      // The pT distributions over all eta are normalized to unit eta (2.0*2.4) and by 1/2*pi*pT.
      // The 1/pT part is taken care of in the filling.  The 1/2pi is taken care of here.
      const double normpT = normfac/(2.0*0.2);
      const double normpTall = normfac/(2.0*M_PI*2.0*2.4);

      for (size_t ietabin=0; ietabin < _h_dNch_dpT.size(); ietabin++){
        scale(_h_dNch_dpT[ietabin], normpT);
      }
      scale(_h_dNch_dpT_all, normpTall);
      scale(_h_dNch_dEta, normfac);
    }


  private:

    std::vector<Histo1DPtr> _h_dNch_dpT;
    Histo1DPtr _h_dNch_dpT_all;
    Histo1DPtr _h_dNch_dEta;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2010_S8656010);

}
#line 1 "CMS_2011_S8884919.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
using namespace std;

namespace Rivet {

  class CMS_2011_S8884919 : public Analysis {
  public:

    CMS_2011_S8884919()
      : Analysis("CMS_2011_S8884919")
    {    }


    void init() {
      ChargedFinalState cfs(-2.4, 2.4, 0.0*GeV);
      declare(cfs, "CFS");

      // eta bins
      _etabins.push_back(0.5);
      _etabins.push_back(1.0);
      _etabins.push_back(1.5);
      _etabins.push_back(2.0);
      _etabins.push_back(2.4) ;

      if (fuzzyEquals(sqrtS()/GeV, 900)) {
        for (size_t ietabin=0; ietabin < _etabins.size(); ietabin++) {
          _h_dNch_dn.push_back( bookHisto1D( 2 + ietabin, 1, 1) );
        }
        _h_dNch_dn_pt500_eta24 = bookHisto1D(20, 1, 1);
        _h_dmpt_dNch_eta24 = bookProfile1D(23, 1, 1);
      }

      if (fuzzyEquals(sqrtS()/GeV, 2360)) {
        for (size_t ietabin=0; ietabin < _etabins.size(); ietabin++) {
          _h_dNch_dn.push_back( bookHisto1D(7 + ietabin, 1, 1) );
        }
        _h_dNch_dn_pt500_eta24 = bookHisto1D(21, 1, 1);
        _h_dmpt_dNch_eta24 = bookProfile1D(24, 1, 1);
      }

      if (fuzzyEquals(sqrtS()/GeV, 7000)) {
        for (size_t ietabin=0; ietabin < _etabins.size(); ietabin++) {
          _h_dNch_dn.push_back( bookHisto1D(12 + ietabin, 1, 1) );
        }
        _h_dNch_dn_pt500_eta24 = bookHisto1D(22, 1, 1);
        _h_dmpt_dNch_eta24 = bookProfile1D(25, 1, 1);
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get the charged particles
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");

      // Resetting the multiplicity for the event to 0;
      vector<int> _nch_in_Evt;
      vector<int> _nch_in_Evt_pt500;
      _nch_in_Evt.assign(_etabins.size(), 0);
      _nch_in_Evt_pt500.assign(_etabins.size(), 0);
      double sumpt = 0;

      // Loop over particles in event
      foreach (const Particle& p, charged.particles()) {
        // Selecting only charged hadrons
        if (! PID::isHadron(p.pid())) continue;

        double pT = p.pT();
        double eta = p.eta();
        sumpt += pT;
        for (size_t ietabin = _etabins.size(); ietabin > 0; --ietabin) {
          if (fabs(eta) > _etabins[ietabin-1]) break;
          ++_nch_in_Evt[ietabin-1];
          if (pT > 0.5/GeV) ++_nch_in_Evt_pt500[ietabin-1];
        }
      }

      // Filling multiplicity-dependent histogramms
      for (size_t ietabin = 0; ietabin < _etabins.size(); ietabin++) {
        _h_dNch_dn[ietabin]->fill(_nch_in_Evt[ietabin], weight);
      }

      // Do only if eta bins are the needed ones
      if (_etabins[4] == 2.4 && _etabins[0] == 0.5) {
        if (_nch_in_Evt[4] != 0) {
          _h_dmpt_dNch_eta24->fill(_nch_in_Evt[4], sumpt/GeV / _nch_in_Evt[4], weight);
        }
        _h_dNch_dn_pt500_eta24->fill(_nch_in_Evt_pt500[4], weight);
      } else {
        MSG_WARNING("You changed the number of eta bins, but forgot to propagate it everywhere !!");
      }
    }


    void finalize() {
      for (size_t ietabin = 0; ietabin < _etabins.size(); ietabin++){
        normalize(_h_dNch_dn[ietabin]);
      }
      normalize(_h_dNch_dn_pt500_eta24);
    }


  private:

    vector<Histo1DPtr> _h_dNch_dn;
    Histo1DPtr _h_dNch_dn_pt500_eta24;
    Profile1DPtr _h_dmpt_dNch_eta24;

    vector<double> _etabins;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S8884919);

}
#line 1 "CMS_2011_S8941262.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  class CMS_2011_S8941262 : public Analysis {
  public:

    /// Constructor
    CMS_2011_S8941262() : Analysis("CMS_2011_S8941262") {  }


    /// Book histograms and initialise projections before the run
    void init() {
      _h_total = bookHisto1D(1, 1, 1);
      _h_mupt  = bookHisto1D(2, 1, 1);
      _h_mueta = bookHisto1D(3, 1, 1);
      nbtot=0.;   nbmutot=0.;

      IdentifiedFinalState ifs(Cuts::abseta < 2.1 && Cuts::pT > 6*GeV);
      ifs.acceptIdPair(PID::MUON);
      declare(ifs, "IFS");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // a b-quark must have been produced
      /// @todo Ouch. Use hadron tagging...
      int nb = 0;
      foreach (const GenParticle* p, particles(event.genEvent())) {
        if (abs(p->pdg_id()) == PID::BQUARK) nb += 1;
      }
      if (nb == 0) vetoEvent;
      nbtot += weight;

      // Event must contain a muon
      Particles muons = apply<IdentifiedFinalState>(event, "IFS").particlesByPt();
      if (muons.size() < 1) vetoEvent;
      nbmutot += weight;

      FourMomentum pmu = muons[0].momentum();
      _h_total->fill(      7000/GeV, weight);
      _h_mupt->fill(   pmu.pT()/GeV, weight);
      _h_mueta->fill( pmu.eta()/GeV, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_total, crossSection()/microbarn/sumOfWeights());
      scale(_h_mupt,  crossSection()/nanobarn/sumOfWeights());
      scale(_h_mueta, crossSection()/nanobarn/sumOfWeights());
    }


  private:

    double nbtot, nbmutot;

    Histo1DPtr _h_total;
    Histo1DPtr _h_mupt;
    Histo1DPtr _h_mueta;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S8941262);

}
#line 1 "CMS_2011_S8950903.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  // CMS azimuthal decorrelations
  class CMS_2011_S8950903 : public Analysis {
  public:

    CMS_2011_S8950903() : Analysis("CMS_2011_S8950903") {}


    void init() {
      FinalState fs;
      FastJets akt(fs, FastJets::ANTIKT, 0.5);
      declare(akt, "antikT");

      _h_deltaPhi.addHistogram( 80.,  110., bookHisto1D(1, 1, 1));
      _h_deltaPhi.addHistogram(110.,  140., bookHisto1D(2, 1, 1));
      _h_deltaPhi.addHistogram(140.,  200., bookHisto1D(3, 1, 1));
      _h_deltaPhi.addHistogram(200.,  300., bookHisto1D(4, 1, 1));
      _h_deltaPhi.addHistogram(300., 7000., bookHisto1D(5, 1, 1));
    }


    void analyze(const Event & event) {
      const double weight = event.weight();

      const Jets& jets = apply<JetAlg>(event, "antikT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      if (fabs(jets[0].eta()) > 1.1 || jets[0].pT() < 80.) vetoEvent;
      if (fabs(jets[1].eta()) > 1.1 || jets[1].pT() < 30.) vetoEvent;

      double dphi = deltaPhi(jets[0].momentum(), jets[1].phi());

      _h_deltaPhi.fill(jets[0].pT(), dphi, weight);
    }


    void finalize() {
      foreach (Histo1DPtr histo, _h_deltaPhi.getHistograms()) {
        normalize(histo, 1.);
      }
    }

  private:

    BinnedHistogram<double> _h_deltaPhi;

  };

  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S8950903);

}

#line 1 "CMS_2011_S8957746.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {

  /// Rivet analysis class for CMS_2011_S8957746 dataset
  class CMS_2011_S8957746 : public Analysis {
  public:

    /// Constructor
    CMS_2011_S8957746()
      : Analysis("CMS_2011_S8957746") {  }


    /// Initialization, called once before running
    void init() {
      // Projections
      const FastJets jets(FinalState(-5.0, 5.0, 0.0*GeV), FastJets::ANTIKT, 0.5);
      declare(jets, "Jets");

      // Book histograms
      _hist_T_90  = bookHisto1D(1, 1, 1);
      _hist_m_90  = bookHisto1D(2, 1, 1);
      _hist_T_125 = bookHisto1D(3, 1, 1);
      _hist_m_125 = bookHisto1D(4, 1, 1);
      _hist_T_200 = bookHisto1D(5, 1, 1);
      _hist_m_200 = bookHisto1D(6, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(30.0*GeV);
      if (jets.size() < 2 ||
          fabs(jets[0].eta()) >= 1.3 ||
          fabs(jets[1].eta()) >= 1.3 ||
          jets[0].pT() < 90*GeV) {
        vetoEvent;
      }
      std::vector<Vector3> momenta;
      foreach (const Jet& j, jets) {
        if (j.abseta() < 1.3) {
          Vector3 mom = j.p3();
          mom.setZ(0.0);
          momenta.push_back(mom);
        }
      }
      if (momenta.size() == 2) {
        // We need to use a ghost so that Thrust.calc() doesn't return 1.
        momenta.push_back(Vector3(1e-10*MeV, 0., 0.));
      }
      Thrust thrust;
      thrust.calc(momenta);

      // The lowest bin also includes the underflow:
      const double T = max(log(1-thrust.thrust()), -12.0);
      const double M = max(log(thrust.thrustMajor()), -6.0);
      if (jets[0].pT()/GeV > 200) {
        _hist_T_200->fill(T, weight);
        _hist_m_200->fill(M, weight);
      } else if (jets[0].pT()/GeV > 125) {
        _hist_T_125->fill(T, weight);
        _hist_m_125->fill(M, weight);
      } else if (jets[0].pT()/GeV > 90) {
        _hist_T_90->fill(T, weight);
        _hist_m_90->fill(M, weight);
      }
    }


    void finalize() {
      normalize(_hist_T_90);
      normalize(_hist_m_90);
      normalize(_hist_T_125);
      normalize(_hist_m_125);
      normalize(_hist_T_200);
      normalize(_hist_m_200);
    }


  private:

    Histo1DPtr _hist_T_90;
    Histo1DPtr _hist_m_90;
    Histo1DPtr _hist_T_125;
    Histo1DPtr _hist_m_125;
    Histo1DPtr _hist_T_200;
    Histo1DPtr _hist_m_200;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S8957746);

}
#line 1 "CMS_2011_S8968497.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CMS_2011_S8968497 : public Analysis {
  public:

    CMS_2011_S8968497()
      : Analysis("CMS_2011_S8968497")
    { }


    void init() {
      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");
      _h_chi_dijet.addHistogram(2200., 7000., bookHisto1D(1, 1, 1));
      _h_chi_dijet.addHistogram(1800., 2200., bookHisto1D(2, 1, 1));
      _h_chi_dijet.addHistogram(1400., 1800., bookHisto1D(3, 1, 1));
      _h_chi_dijet.addHistogram(1100., 1400., bookHisto1D(4, 1, 1));
      _h_chi_dijet.addHistogram( 850., 1100., bookHisto1D(5, 1, 1));
      _h_chi_dijet.addHistogram( 650.,  850., bookHisto1D(6, 1, 1));
      _h_chi_dijet.addHistogram( 500.,  650., bookHisto1D(7, 1, 1));
      _h_chi_dijet.addHistogram( 350.,  500., bookHisto1D(8, 1, 1));
      _h_chi_dijet.addHistogram( 250.,  350., bookHisto1D(9, 1, 1));
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;
      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      double y0 = j0.rapidity();
      double y1 = j1.rapidity();
      if (fabs(y0+y1)/2. > 1.11) vetoEvent;
      double mjj = FourMomentum(j0+j1).mass();
      double chi = exp(fabs(y0-y1));
      if(chi<16.)  _h_chi_dijet.fill(mjj, chi, weight);
    }


    void finalize() {
      foreach (Histo1DPtr hist, _h_chi_dijet.getHistograms()) {
        normalize(hist);
      }
    }


  private:

    BinnedHistogram<double> _h_chi_dijet;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S8968497);

}
#line 1 "CMS_2011_S8973270.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CMS_2011_S8973270 : public Analysis {
  public:

    /// Constructor
    CMS_2011_S8973270() : Analysis("CMS_2011_S8973270") {  }


    void init() {
      FinalState fs;
      FastJets jetproj(fs, FastJets::ANTIKT, 0.5);
      jetproj.useInvisibles();
      declare(jetproj, "Jets");

      UnstableFinalState ufs;
      declare(ufs, "UFS");

      // Book histograms
      _h_dsigma_dR_56GeV = bookHisto1D(1,1,1);
      _h_dsigma_dR_84GeV = bookHisto1D(2,1,1);
      _h_dsigma_dR_120GeV = bookHisto1D(3,1,1);
      _h_dsigma_dPhi_56GeV = bookHisto1D(4,1,1);
      _h_dsigma_dPhi_84GeV = bookHisto1D(5,1,1);
      _h_dsigma_dPhi_120GeV = bookHisto1D(6,1,1);

      _countMCDR56 = 0;
      _countMCDR84 = 0;
      _countMCDR120 = 0;
      _countMCDPhi56 = 0;
      _countMCDPhi84 = 0;
      _countMCDPhi120 = 0;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const Jets& jets = apply<FastJets>(event,"Jets").jetsByPt();
      const UnstableFinalState& ufs = apply<UnstableFinalState>(event, "UFS");

      // Find the leading jet pT and eta
      if (jets.size() == 0) vetoEvent;
      const double ljpT = jets[0].pT();
      const double ljeta = jets[0].eta();
      MSG_DEBUG("Leading jet pT / eta: " << ljpT << " / " << ljeta);

      // Minimum requirement for event
      if (ljpT > 56*GeV && fabs(ljeta) < 3.0) {
        // Find B hadrons in event
        int nab = 0, nb = 0; //counters for all B and independent B hadrons
        double etaB1 = 7.7, etaB2 = 7.7;
        double phiB1 = 7.7, phiB2 = 7.7;
        double pTB1 = 7.7, pTB2 = 7.7;

        foreach (const Particle& p, ufs.particles()) {
          int aid = p.abspid();
          if (aid/100 == 5 || aid/1000==5) {
            nab++;
            // 2J+1 == 1 (mesons) or 2 (baryons)
            if (aid%10 == 1 || aid%10 == 2) {
              // No B decaying to B
              if (aid != 5222 && aid != 5112 && aid != 5212 && aid != 5322) {
                if (nb==0) {
                  etaB1 = p.eta();
                  phiB1 = p.phi();
                  pTB1 = p.pT();
                } else if (nb==1) {
                  etaB2 = p.eta();
                  phiB2 = p.phi();
                  pTB2 = p.pT();
                }
                nb++;
              }
            }
            MSG_DEBUG("ID " << aid <<  " B hadron");
          }
        }

        if (nb==2 && pTB1 > 15*GeV && pTB2 > 15*GeV && fabs(etaB1) < 2.0 && fabs(etaB2) < 2.0) {
          double dPhi = deltaPhi(phiB1, phiB2);
          double dR = deltaR(etaB1, phiB1, etaB2, phiB2);
          MSG_DEBUG("DR/DPhi " << dR << " " << dPhi);

          // MC counters
          if (dR > 2.4) _countMCDR56 += weight;
          if (dR > 2.4 && ljpT > 84*GeV) _countMCDR84 += weight;
          if (dR > 2.4 && ljpT > 120*GeV) _countMCDR120 += weight;
          if (dPhi > 3.*PI/4.) _countMCDPhi56 += weight;
          if (dPhi > 3.*PI/4. && ljpT > 84*GeV) _countMCDPhi84 += weight;
          if (dPhi > 3.*PI/4. && ljpT > 120*GeV) _countMCDPhi120 += weight;

          _h_dsigma_dR_56GeV->fill(dR, weight);
          if (ljpT > 84*GeV) _h_dsigma_dR_84GeV->fill(dR, weight);
          if (ljpT > 120*GeV) _h_dsigma_dR_120GeV->fill(dR, weight);
          _h_dsigma_dPhi_56GeV->fill(dPhi, weight);
          if (ljpT > 84*GeV) _h_dsigma_dPhi_84GeV->fill(dPhi, weight);
          if (ljpT > 120*GeV) _h_dsigma_dPhi_120GeV->fill(dPhi, weight);
          //MSG_DEBUG("nb " << nb << " " << nab);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      MSG_DEBUG("crossSection " << crossSection() << " sumOfWeights " << sumOfWeights());

      // Hardcoded bin widths
      double DRbin = 0.4;
      double DPhibin = PI/8.0;
      // Find out the correct numbers
      double nDataDR56 = 25862.20;
      double nDataDR84 = 5675.55;
      double nDataDR120 = 1042.72;
      double nDataDPhi56 = 24220.00;
      double nDataDPhi84 = 4964.00;
      double nDataDPhi120 = 919.10;
      double normDR56 = (_countMCDR56 > 0.) ? nDataDR56/_countMCDR56 : crossSection()/sumOfWeights();
      double normDR84 = (_countMCDR84 > 0.) ? nDataDR84/_countMCDR84 : crossSection()/sumOfWeights();
      double normDR120 = (_countMCDR120 > 0.) ? nDataDR120/_countMCDR120 : crossSection()/sumOfWeights();
      double normDPhi56 = (_countMCDPhi56 > 0.) ? nDataDPhi56/_countMCDPhi56 : crossSection()/sumOfWeights();
      double normDPhi84 = (_countMCDPhi84 > 0.) ? nDataDPhi84/_countMCDPhi84 : crossSection()/sumOfWeights();
      double normDPhi120 = (_countMCDPhi120 > 0.) ? nDataDPhi120/_countMCDPhi120 : crossSection()/sumOfWeights();
      scale(_h_dsigma_dR_56GeV, normDR56*DRbin);
      scale(_h_dsigma_dR_84GeV, normDR84*DRbin);
      scale(_h_dsigma_dR_120GeV, normDR120*DRbin);
      scale(_h_dsigma_dPhi_56GeV, normDPhi56*DPhibin);
      scale(_h_dsigma_dPhi_84GeV, normDPhi84*DPhibin);
      scale(_h_dsigma_dPhi_120GeV, normDPhi120*DPhibin);
    }

    //@}


  private:

    /// @name Counters
    //@{
    double _countMCDR56, _countMCDR84, _countMCDR120;
    double _countMCDPhi56, _countMCDPhi84, _countMCDPhi120;
    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _h_dsigma_dR_56GeV, _h_dsigma_dR_84GeV, _h_dsigma_dR_120GeV;
    Histo1DPtr _h_dsigma_dPhi_56GeV, _h_dsigma_dPhi_84GeV, _h_dsigma_dPhi_120GeV;
    //@}

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S8973270);

}
#line 1 "CMS_2011_S8978280.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief CMS strange particle spectra (Ks, Lambda, Cascade) in pp at 900 and 7000 GeV
  /// @author Kevin Stenson
  class CMS_2011_S8978280 : public Analysis {
  public:

    /// Constructor
    CMS_2011_S8978280()
      : Analysis("CMS_2011_S8978280")
    { }


    void init() {
      UnstableFinalState ufs(Cuts::absrap < 2);
      declare(ufs, "UFS");

      // Particle distributions versus rapidity and transverse momentum
      if (fuzzyEquals(sqrtS()/GeV, 900*GeV)){
        _h_dNKshort_dy  = bookHisto1D(1, 1, 1);
        _h_dNKshort_dpT = bookHisto1D(2, 1, 1);
        _h_dNLambda_dy  = bookHisto1D(3, 1, 1);
        _h_dNLambda_dpT = bookHisto1D(4, 1, 1);
        _h_dNXi_dy      = bookHisto1D(5, 1, 1);
        _h_dNXi_dpT     = bookHisto1D(6, 1, 1);
        //
        _h_LampT_KpT    = bookScatter2D(7, 1, 1);
        _h_XipT_LampT   = bookScatter2D(8, 1, 1);
        _h_Lamy_Ky      = bookScatter2D(9, 1, 1);
        _h_Xiy_Lamy     = bookScatter2D(10, 1, 1);

      } else if (fuzzyEquals(sqrtS()/GeV, 7000*GeV)){
        _h_dNKshort_dy  = bookHisto1D(1, 1, 2);
        _h_dNKshort_dpT = bookHisto1D(2, 1, 2);
        _h_dNLambda_dy  = bookHisto1D(3, 1, 2);
        _h_dNLambda_dpT = bookHisto1D(4, 1, 2);
        _h_dNXi_dy      = bookHisto1D(5, 1, 2);
        _h_dNXi_dpT     = bookHisto1D(6, 1, 2);
        //
        _h_LampT_KpT    = bookScatter2D(7, 1, 2);
        _h_XipT_LampT   = bookScatter2D(8, 1, 2);
        _h_Lamy_Ky      = bookScatter2D(9, 1, 2);
        _h_Xiy_Lamy     = bookScatter2D(10, 1, 2);
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      const UnstableFinalState& parts = apply<UnstableFinalState>(event, "UFS");
      foreach (const Particle& p, parts.particles()) {
        switch (p.abspid()) {
        case PID::K0S:
          _h_dNKshort_dy->fill(p.absrap(), weight);
          _h_dNKshort_dpT->fill(p.pT(), weight);
          break;

        case PID::LAMBDA:
          // Lambda should not have Cascade or Omega ancestors since they should not decay. But just in case...
          if ( !( p.hasAncestor(3322) || p.hasAncestor(-3322) || p.hasAncestor(3312) || p.hasAncestor(-3312) || p.hasAncestor(3334) || p.hasAncestor(-3334) ) ) {
            _h_dNLambda_dy->fill(p.absrap(), weight);
            _h_dNLambda_dpT->fill(p.pT(), weight);
          }
          break;

        case PID::XIMINUS:
          // Cascade should not have Omega ancestors since it should not decay.  But just in case...
          if ( !( p.hasAncestor(3334) || p.hasAncestor(-3334) ) ) {
            _h_dNXi_dy->fill(p.absrap(), weight);
            _h_dNXi_dpT->fill(p.pT(), weight);
          }
          break;
        }

      }
    }


    void finalize() {
      divide(_h_dNLambda_dpT,_h_dNKshort_dpT, _h_LampT_KpT);
      divide(_h_dNXi_dpT,_h_dNLambda_dpT, _h_XipT_LampT);
      divide(_h_dNLambda_dy,_h_dNKshort_dy, _h_Lamy_Ky);
      divide(_h_dNXi_dy,_h_dNLambda_dy, _h_Xiy_Lamy);
      const double normpT = 1.0/sumOfWeights();
      const double normy = 0.5*normpT; // Accounts for using |y| instead of y
      scale(_h_dNKshort_dy, normy);
      scale(_h_dNKshort_dpT, normpT);
      scale(_h_dNLambda_dy, normy);
      scale(_h_dNLambda_dpT, normpT);
      scale(_h_dNXi_dy, normy);
      scale(_h_dNXi_dpT, normpT);
    }


  private:

    // Particle distributions versus rapidity and transverse momentum
    Histo1DPtr _h_dNKshort_dy, _h_dNKshort_dpT, _h_dNLambda_dy, _h_dNLambda_dpT, _h_dNXi_dy, _h_dNXi_dpT;
    Scatter2DPtr _h_LampT_KpT, _h_XipT_LampT, _h_Lamy_Ky, _h_Xiy_Lamy;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S8978280);

}
#line 1 "CMS_2011_S9086218.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  // Inclusive jet pT
  class CMS_2011_S9086218 : public Analysis {
  public:

    // Constructor
    CMS_2011_S9086218() : Analysis("CMS_2011_S9086218") {}


    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;

      // Initialize the projectors:
      declare(FastJets(fs, FastJets::ANTIKT, 0.5),"Jets");

      // Book histograms:
      _hist_sigma.addHistogram(0.0, 0.5, bookHisto1D(1, 1, 1));
      _hist_sigma.addHistogram(0.5, 1.0, bookHisto1D(2, 1, 1));
      _hist_sigma.addHistogram(1.0, 1.5, bookHisto1D(3, 1, 1));
      _hist_sigma.addHistogram(1.5, 2.0, bookHisto1D(4, 1, 1));
      _hist_sigma.addHistogram(2.0, 2.5, bookHisto1D(5, 1, 1));
      _hist_sigma.addHistogram(2.5, 3.0, bookHisto1D(6, 1, 1));
    }

    // Analysis
    void analyze(const Event &event) {
      const double weight = event.weight();
      const FastJets& fj = apply<FastJets>(event,"Jets");
      const Jets& jets = fj.jets(Cuts::ptIn(18*GeV, 1100.0*GeV) && Cuts::absrap < 4.7);

      // Fill the relevant histograms:
      foreach(const Jet& j, jets) {
        _hist_sigma.fill(j.absrap(), j.pT(), weight);
      }
    }

    // Finalize
    void finalize() {
      _hist_sigma.scale(crossSection()/sumOfWeights()/2.0, this);
    }

  private:
    BinnedHistogram<double> _hist_sigma;
  };

  // This global object acts as a hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_2011_S9086218);

}
#line 1 "CMS_2011_S9088458.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


   /// CMS ratio of 3-jet to 2-jet cross-sections
   class CMS_2011_S9088458 : public Analysis {
   public:


     CMS_2011_S9088458()
       : Analysis("CMS_2011_S9088458") {  }


     void init() {
       FinalState fs;
       FastJets akt(fs, FastJets::ANTIKT, 0.5);
       declare(akt, "antikT");

       _h_tmp_dijet = Histo1D(refData(1, 1, 1));
       _h_tmp_trijet = Histo1D(refData(1, 1, 1));
       _h_r32 = bookScatter2D(1, 1, 1);
     }


     void analyze(const Event & event) {
       const double weight = event.weight();

       Jets highpT_jets;
       double HT = 0;
       foreach(const Jet & jet, apply<JetAlg>(event, "antikT").jetsByPt(50.0*GeV)) {
         if (jet.abseta() < 2.5) {
           highpT_jets.push_back(jet);
           HT += jet.pT();
         }
       }
       if (highpT_jets.size() < 2) vetoEvent;
       if (highpT_jets.size() >= 2) _h_tmp_dijet.fill(HT/TeV, weight);
       if (highpT_jets.size() >= 3) _h_tmp_trijet.fill(HT/TeV, weight);
     }


     void finalize() {
       divide(_h_tmp_trijet, _h_tmp_dijet, _h_r32);
     }


   private:

     Histo1D _h_tmp_dijet, _h_tmp_trijet;
     Scatter2DPtr _h_r32;

  };


  // A hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S9088458);

}
#line 1 "CMS_2011_S9120041.cc"
// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

using namespace std;

namespace Rivet {


  // UE charged particles vs. leading jet
  class CMS_2011_S9120041 : public Analysis {
  public:

    /// Constructor
    CMS_2011_S9120041() : Analysis("CMS_2011_S9120041") {}


    void init() {
      const ChargedFinalState cfs(-2.0, 2.0, 500*MeV);
      declare(cfs, "CFS");

      const ChargedFinalState cfsforjet(-2.5, 2.5, 500*MeV);
      const FastJets jetpro(cfsforjet, FastJets::SISCONE, 0.5);
      declare(jetpro, "Jets");

      if (fuzzyEquals(sqrtS(), 7.0*TeV)) {
        _h_Nch_vs_pT = bookProfile1D(1, 1, 1); // Nch vs. pT_max
        _h_Sum_vs_pT = bookProfile1D(2, 1, 1); // sum(pT) vs. pT_max
        _h_pT3_Nch   = bookHisto1D(5, 1, 1);   // transverse Nch,     pT_max > 3GeV
        _h_pT3_Sum   = bookHisto1D(6, 1, 1);   // transverse sum(pT), pT_max > 3GeV
        _h_pT3_pT    = bookHisto1D(7, 1, 1);   // transverse pT,      pT_max > 3GeV
        _h_pT20_Nch  = bookHisto1D(8, 1, 1);   // transverse Nch,     pT_max > 20GeV
        _h_pT20_Sum  = bookHisto1D(9, 1, 1);   // transverse sum(pT), pT_max > 20GeV
        _h_pT20_pT   = bookHisto1D(10, 1, 1);  // transverse pT,      pT_max > 20GeV
      }

      if (fuzzyEquals(sqrtS(), 0.9*TeV)) {
        _h_Nch_vs_pT = bookProfile1D(3, 1, 1); // Nch vs. pT_max
        _h_Sum_vs_pT = bookProfile1D(4, 1, 1); // sum(pT) vs. pT_max
        _h_pT3_Nch   = bookHisto1D(11, 1, 1);  // transverse Nch,     pT_max > 3GeV
        _h_pT3_Sum   = bookHisto1D(12, 1, 1);  // transverse sum(pT), pT_max > 3GeV
        _h_pT3_pT    = bookHisto1D(13, 1, 1);  // transverse pT,      pT_max > 3GeV
      }

      sumOfWeights3  = 0.0;
      sumOfWeights20 = 0.0;

      _nch_tot_pT3  = 0.0;
      _nch_tot_pT20 = 0.0;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Find the lead jet, applying a restriction that the jets must be within |eta| < 2.
      FourMomentum p_lead;
      foreach (const Jet& j, apply<FastJets>(event, "Jets").jetsByPt(1.0*GeV)) {
        if (j.abseta() < 2.0) {
          p_lead = j.momentum();
          break;
        }
      }
      if (p_lead.isZero()) vetoEvent;
      const double philead = p_lead.phi();
      const double pTlead  = p_lead.pT();

      Particles particles = apply<ChargedFinalState>(event, "CFS").particlesByPt();

      int nTransverse = 0;
      double ptSumTransverse = 0.;
      foreach (const Particle& p, particles) {
        double dphi = fabs(deltaPhi(philead, p.phi()));
        if (dphi>PI/3. && dphi<PI*2./3.) {   // Transverse region
          nTransverse++;

          const double pT = p.pT()/GeV;
          ptSumTransverse += pT;

          if (pTlead > 3.0*GeV) _h_pT3_pT->fill(pT, weight);
          if (fuzzyEquals(sqrtS(), 7.0*TeV) && pTlead > 20.0*GeV) _h_pT20_pT->fill(pT, weight);
        }
      }

      const double area = 8./3. * PI;
      _h_Nch_vs_pT->fill(pTlead/GeV, 1./area*nTransverse, weight);
      _h_Sum_vs_pT->fill(pTlead/GeV, 1./area*ptSumTransverse, weight);
      if(pTlead > 3.0*GeV) {
        _h_pT3_Nch->fill(nTransverse, weight);
        _h_pT3_Sum->fill(ptSumTransverse, weight);
        sumOfWeights3 += weight;
        _nch_tot_pT3  += weight*nTransverse;
      }
      if (fuzzyEquals(sqrtS(), 7.0*TeV) && pTlead > 20.0*GeV) {
        _h_pT20_Nch->fill(nTransverse, weight);
        _h_pT20_Sum->fill(ptSumTransverse, weight);
        sumOfWeights20 += weight;
        _nch_tot_pT20  += weight*nTransverse;
      }
    }



    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pT3_Nch);
      normalize(_h_pT3_Sum);
      if (sumOfWeights3 != 0.0) normalize(_h_pT3_pT, _nch_tot_pT3 / sumOfWeights3);

      if (fuzzyEquals(sqrtS(), 7.0*TeV)) {
        normalize(_h_pT20_Nch);
        normalize(_h_pT20_Sum);
        if (sumOfWeights20 != 0.0) normalize(_h_pT20_pT, _nch_tot_pT20 / sumOfWeights20);
      }
    }



  private:

    double sumOfWeights3;
    double sumOfWeights20;

    double _nch_tot_pT3;
    double _nch_tot_pT20;

    Profile1DPtr _h_Nch_vs_pT;
    Profile1DPtr _h_Sum_vs_pT;
    Histo1DPtr _h_pT3_Nch;
    Histo1DPtr _h_pT3_Sum;
    Histo1DPtr _h_pT3_pT;
    Histo1DPtr _h_pT20_Nch;
    Histo1DPtr _h_pT20_Sum;
    Histo1DPtr _h_pT20_pT;

  };


  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S9120041);
}

#line 1 "CMS_2011_S9215166.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  class CMS_2011_S9215166 : public Analysis {
  public:

    /// Constructor
    CMS_2011_S9215166() : Analysis("CMS_2011_S9215166"), _weightMB(0.), _weightDiJet(0.) {  }


    void init() {
      const FinalState fs(-6.0, 6.0, 0.0*GeV);
      declare(fs, "FS");
      declare(FastJets(fs, FastJets::ANTIKT, 0.5), "Jets");

      VetoedFinalState fsv(fs);
      fsv.vetoNeutrinos();
      fsv.addVetoPairDetail(PID::MUON, 0.0*GeV, 99999.9*GeV);
      declare(fsv, "fsv");

      // For the MB ND selection
      const ChargedFinalState fschrgd(-6.0,6.0,0.0*GeV);
      declare(fschrgd, "fschrgd");
      VetoedFinalState fschrgdv(fschrgd);
      fschrgdv.vetoNeutrinos();
      declare(fschrgdv, "fschrgdv");

      if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3)) {
        _hist_mb      = bookHisto1D(1, 1, 1); // energy flow in MB, 0.9 TeV
        _hist_dijet = bookHisto1D(2, 1, 1); // energy flow in dijet events, 0.9 TeV
      } else if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) {
        _hist_mb      = bookHisto1D(3, 1, 1); // energy flow in MB, 7 TeV
        _hist_dijet = bookHisto1D(4, 1, 1); // energy flow in dijet events, 7 TeV
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      // Skip if the event is empty
      const FinalState& fsv = apply<FinalState>(event, "fsv");
      if (fsv.empty()) vetoEvent;

      // Veto diffractive topologies according to defined hadron level
      double count_chrg_forward = 0;
      double count_chrg_backward = 0;
      const FinalState& fschrgdv = apply<FinalState>(event, "fschrgdv");
      foreach (const Particle& p, fschrgdv.particles()) {
        if (3.9 < p.eta() && p.eta() < 4.4) count_chrg_forward++;
        if (-4.4 < p.eta() && p.eta() < -3.9) count_chrg_backward++;
      }
      if (count_chrg_forward == 0 || count_chrg_backward == 0) vetoEvent;
      /// @todo "Diffractive" veto should really also veto dijet events?


      // MINIMUM BIAS EVENTS
      _weightMB += weight;
      foreach (const Particle& p, fsv.particles()) {
        _hist_mb->fill(p.abseta(), weight*p.E()/GeV);
      }


      // DIJET EVENTS
      double PTCUT = -1.0;
      if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3)) PTCUT = 8.0*GeV;
      else if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) PTCUT = 20.0*GeV;
      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets jets = jetpro.jetsByPt(PTCUT);
      if (jets.size() >= 2) {
        // eta cut for the central jets
        if (fabs(jets[0].eta()) < 2.5 && fabs(jets[1].eta()) < 2.5) {
          // Back to back condition of the jets
          const double diffphi = deltaPhi(jets[1].phi(), jets[0].phi());
          if (diffphi-PI < 1.0) {
	    _weightDiJet += weight;
            foreach (const Particle& p, fsv.particles()) {
              _hist_dijet->fill(p.abseta(), weight*p.E()/GeV);
            }
          }
        }
      }

    }


    void finalize() {
      scale(_hist_mb   , 0.5/_weightMB   );
      scale(_hist_dijet, 0.5/_weightDiJet);
    }


  private:

    Histo1DPtr _hist_mb, _hist_dijet;
    double _weightMB,_weightDiJet;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_S9215166);

}
#line 1 "CMS_2012_I941555.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief CMS Z pT and rapidity in Drell-Yan events at 7 TeV
  /// @author Justin Hugon, Luca Perrozzi
  class CMS_2012_I941555 : public Analysis {
  public:

    /// Constructor
    CMS_2012_I941555()
      : Analysis("CMS_2012_I941555")
    {
      _sumw_mu_dressed_pt  = 0;
      _sumwpeak_mu_dressed = 0;
      _sumw_el_dressed_rap = 0;
      _sumw_el_dressed_pt  = 0;
      _sumwpeak_el_dressed = 0;
    }


    /// @name Analysis methods
    //@{

    void init() {

      // Set up projections
      /// @todo Really?: ZFinder zfinder_dressed_mu_pt(-2.1, 2.1, 20, PID::MUON, 60*GeV, 120*GeV, 0.2, false, true);
      FinalState fs;
      Cut cuts = Cuts::abseta < 2.1 && Cuts::pT > 20*GeV;
      ZFinder zfinder_dressed_mu_pt(fs, cuts, PID::MUON, 60*GeV, 120*GeV, 0.2);
      declare(zfinder_dressed_mu_pt, "ZFinder_dressed_mu_pt");
      ZFinder zfinder_dressed_el_pt(fs, cuts, PID::ELECTRON, 60*GeV, 120*GeV, 0.1);
      declare(zfinder_dressed_el_pt, "ZFinder_dressed_el_pt");

      ZFinder zfinder_dressed_mu_rap(fs, Cuts::open(), PID::MUON, 60*GeV, 120*GeV, 0.1);
      declare(zfinder_dressed_mu_rap, "ZFinder_dressed_mu_rap");
      ZFinder zfinder_dressed_el_rap(fs, Cuts::open(), PID::ELECTRON, 60*GeV, 120*GeV, 0.1);
      declare(zfinder_dressed_el_rap, "ZFinder_dressed_el_rap");

      // Book histograms
      _hist_zrap_mu_dressed      = bookHisto1D(1, 1, 1);  // muon "dressed" rapidity
      _hist_zrap_el_dressed      = bookHisto1D(1, 1, 2);  // electron "dressed" rapidity
      _hist_zrap_comb_dressed    = bookHisto1D(1, 1, 3);  // electron "dressed" rapidity

      _hist_zpt_mu_dressed       = bookHisto1D(2, 1, 1);  // muon "dressed" pt
      _hist_zpt_el_dressed       = bookHisto1D(2, 1, 2);  // electron "dressed" pt
      _hist_zpt_comb_dressed     = bookHisto1D(2, 1, 3);  // electron "dressed" pt

      _hist_zptpeak_mu_dressed   = bookHisto1D(3, 1, 1);  // muon "dressed" pt peak
      _hist_zptpeak_el_dressed   = bookHisto1D(3, 1, 2);  // electron "dressed" pt peak
      _hist_zptpeak_comb_dressed = bookHisto1D(3, 1, 3);  // electron "dressed" pt peak
    }


    /// Do the analysis
    void analyze(const Event& evt) {
      const double weight = evt.weight();

      const ZFinder& zfinder_dressed_mu_rap = apply<ZFinder>(evt, "ZFinder_dressed_mu_rap");
      if (!zfinder_dressed_mu_rap.bosons().empty()) {
        _sumw_mu_dressed_rap += weight;
        const FourMomentum pZ = zfinder_dressed_mu_rap.bosons()[0].momentum();
        _hist_zrap_mu_dressed->fill(pZ.rapidity()/GeV, weight);
        _hist_zrap_comb_dressed->fill(pZ.rapidity()/GeV, weight);
      }

      const ZFinder& zfinder_dressed_mu_pt = apply<ZFinder>(evt, "ZFinder_dressed_mu_pt");
      if (!zfinder_dressed_mu_pt.bosons().empty()) {
        _sumw_mu_dressed_pt += weight;
        const FourMomentum pZ = zfinder_dressed_mu_pt.bosons()[0].momentum();
        _hist_zpt_mu_dressed->fill(pZ.pT()/GeV, weight);
        _hist_zpt_comb_dressed->fill(pZ.pT()/GeV, weight);
        if (pZ.pT() < 30*GeV) {
          _sumwpeak_mu_dressed += weight;
          _hist_zptpeak_mu_dressed->fill(pZ.pT()/GeV, weight);
          _hist_zptpeak_comb_dressed->fill(pZ.pT()/GeV, weight);
        }
      }

      const ZFinder& zfinder_dressed_el_rap = apply<ZFinder>(evt, "ZFinder_dressed_el_rap");
      if (!zfinder_dressed_el_rap.bosons().empty()) {
        _sumw_el_dressed_rap += weight;
        const FourMomentum pZ = zfinder_dressed_el_rap.bosons()[0].momentum();
        _hist_zrap_el_dressed->fill(pZ.rapidity()/GeV, weight);
        _hist_zrap_comb_dressed->fill(pZ.rapidity()/GeV, weight);
      }

      const ZFinder& zfinder_dressed_el_pt = apply<ZFinder>(evt, "ZFinder_dressed_el_pt");
      if (!zfinder_dressed_el_pt.bosons().empty()) {
        _sumw_el_dressed_pt += weight;
        const FourMomentum pZ = zfinder_dressed_el_pt.bosons()[0].momentum();
        _hist_zpt_el_dressed->fill(pZ.pT()/GeV, weight);
        _hist_zpt_comb_dressed->fill(pZ.pT()/GeV, weight);
        if (pZ.pT() < 30*GeV) {
          _sumwpeak_el_dressed += weight;
          _hist_zptpeak_el_dressed->fill(pZ.pT()/GeV, weight);
          _hist_zptpeak_comb_dressed->fill(pZ.pT()/GeV, weight);
        }
      }

    }


    void finalize() {
      scale(_hist_zrap_mu_dressed, safediv(1, _sumw_mu_dressed_rap, 1));
      scale(_hist_zpt_mu_dressed, safediv(1, _sumw_mu_dressed_pt, 1));
      scale(_hist_zptpeak_mu_dressed, safediv(1, _sumwpeak_mu_dressed, 1));

      scale(_hist_zrap_el_dressed, safediv(1, _sumw_el_dressed_rap, 1));
      scale(_hist_zpt_el_dressed, safediv(1, _sumw_el_dressed_pt, 1));
      scale(_hist_zptpeak_el_dressed, safediv(1, _sumwpeak_el_dressed, 1));

      scale(_hist_zrap_comb_dressed, safediv(1, _sumw_el_dressed_rap+_sumw_mu_dressed_rap, 1));
      scale(_hist_zpt_comb_dressed, safediv(1, _sumw_el_dressed_pt+_sumw_mu_dressed_pt, 1));
      scale(_hist_zptpeak_comb_dressed, safediv(1, _sumwpeak_el_dressed+_sumwpeak_mu_dressed, 1));
    }

    //@}


  private:

    double _sumw_mu_dressed_rap;
    double _sumw_mu_dressed_pt;
    double _sumwpeak_mu_dressed;

    double _sumw_el_dressed_rap;
    double _sumw_el_dressed_pt;
    double _sumwpeak_el_dressed;

    Histo1DPtr _hist_zrap_mu_dressed;
    Histo1DPtr _hist_zpt_mu_dressed;
    Histo1DPtr _hist_zptpeak_mu_dressed;

    Histo1DPtr _hist_zrap_el_dressed;
    Histo1DPtr _hist_zpt_el_dressed;
    Histo1DPtr _hist_zptpeak_el_dressed;

    Histo1DPtr _hist_zrap_comb_dressed;
    Histo1DPtr _hist_zpt_comb_dressed;
    Histo1DPtr _hist_zptpeak_comb_dressed;

  };


  DECLARE_RIVET_PLUGIN(CMS_2012_I941555);

}
#line 1 "CMS_2011_I954992.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  class CMS_2011_I954992 : public Analysis {
  public:

    CMS_2011_I954992()
      : Analysis("CMS_2011_I954992")
    {    }


  public:

    void init() {
      ChargedFinalState cfs(Cuts::abseta < 2.4);
      declare(cfs,"CFS");

      /// Get muons which pass the initial kinematic cuts
      IdentifiedFinalState muon_fs(Cuts::abseta < 2.1 && Cuts::pT > 4*GeV);
      muon_fs.acceptIdPair(PID::MUON);
      declare(muon_fs, "MUON_FS");

      _h_sigma = bookHisto1D(1,1,1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() != 2) vetoEvent; // no other charged particles in 2.4

      const Particles& muonFS = apply<IdentifiedFinalState>(event, "MUON_FS").particles();
      if (muonFS.size() != 2) vetoEvent;

      if (charge(muonFS[0]) != charge(muonFS[1])) {
         const double dimuon_mass = (muonFS[0].momentum() + muonFS[1].momentum()).mass();
         const double v_angle     = muonFS[0].momentum().angle(muonFS[1].momentum());
         const double dPhi        = deltaPhi(muonFS[0], muonFS[1]);
         const double deltaPt     = fabs(muonFS[0].pT() - muonFS[1].pT());

         if (dimuon_mass >= 11.5*GeV &&
             v_angle < 0.95*PI       &&
             dPhi    > 0.9*PI        &&
             deltaPt < 1.*GeV        ) {
           _h_sigma->fill(sqrtS()/GeV, weight);
         }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_sigma, crossSection()/picobarn/sumOfWeights());
    }

  private:

    Histo1DPtr _h_sigma;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2011_I954992);

}
#line 1 "CMS_2012_I1087342.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"


namespace Rivet {

  // This analysis is a derived from the class Analysis:
  class CMS_2012_I1087342 : public Analysis {

  public:

    // Constructor
    CMS_2012_I1087342() : Analysis("CMS_2012_I1087342") {
    }

    void init() {
      const FinalState fs;
      declare(FastJets(fs, FastJets::ANTIKT, 0.5),"Jets");

      _hist_jetpt_fwdincl = bookHisto1D(1, 1, 1);
      _hist_jetpt_forward = bookHisto1D(2, 1, 1);
      _hist_jetpt_central = bookHisto1D(3, 1, 1);
    }

    void analyze(const Event& event) {
      const double weight = event.weight();

      const FastJets& fj = apply<FastJets>(event,"Jets");
      const Jets jets = fj.jets(Cuts::ptIn(35*GeV, 150*GeV) && Cuts::abseta < 4.7);

      double cjet_pt = 0.0;
      double fjet_pt = 0.0;

      foreach(const Jet& j, jets) {
        double pT = j.pT();
        if (j.abseta() > 3.2) {
          _hist_jetpt_fwdincl->fill(j.pT()/GeV, weight);
        }
        if (j.abseta() < 2.8) {
          if (cjet_pt < pT) cjet_pt = pT;
        }
        if (inRange(j.abseta(), 3.2, 4.7)) {
          if (fjet_pt < pT) fjet_pt = pT;
        }
      }

      if (cjet_pt > 35*GeV && fjet_pt > 35*GeV) {
        _hist_jetpt_forward->fill(fjet_pt/GeV, weight);
        _hist_jetpt_central->fill(cjet_pt/GeV, weight);
      }

    }


    void finalize() {
      scale(_hist_jetpt_fwdincl, crossSection() / picobarn / sumOfWeights() / 3.0);
      scale(_hist_jetpt_forward, crossSection() / picobarn / sumOfWeights() / 3.0);
      scale(_hist_jetpt_central, crossSection() / picobarn / sumOfWeights() / 5.6);
    }


  private:

    Histo1DPtr _hist_jetpt_fwdincl;
    Histo1DPtr _hist_jetpt_forward;
    Histo1DPtr _hist_jetpt_central;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1087342);

}
#line 1 "CMS_2012_I1090423.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CMS_2012_I1090423 : public Analysis {
  public:

    CMS_2012_I1090423()
      : Analysis("CMS_2012_I1090423")
    { }


    void init() {
      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");
      _h_chi_dijet.addHistogram(3000, 7000, bookHisto1D(1, 1, 1));
      _h_chi_dijet.addHistogram(2400, 3000, bookHisto1D(2, 1, 1));
      _h_chi_dijet.addHistogram(1900, 2400, bookHisto1D(3, 1, 1));
      _h_chi_dijet.addHistogram(1500, 1900, bookHisto1D(4, 1, 1));
      _h_chi_dijet.addHistogram(1200, 1500, bookHisto1D(5, 1, 1));
      _h_chi_dijet.addHistogram(1000, 1200, bookHisto1D(6, 1, 1));
      _h_chi_dijet.addHistogram( 800, 1000, bookHisto1D(7, 1, 1));
      _h_chi_dijet.addHistogram( 600,  800, bookHisto1D(8, 1, 1));
      _h_chi_dijet.addHistogram( 400,  600, bookHisto1D(9, 1, 1));
    }


    void analyze(const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      const double y0 = jets[0].rapidity();
      const double y1 = jets[1].rapidity();
      if (fabs(y0+y1)/2 > 1.11) vetoEvent;

      const double chi = exp(fabs(y0-y1));
      if (chi > 16) vetoEvent;

      const FourMomentum jj = jets[0].momentum() + jets[1].momentum();
       _h_chi_dijet.fill(jj.mass(), chi, event.weight());
    }


    void finalize() {
      foreach (Histo1DPtr hist, _h_chi_dijet.getHistograms()) {
        normalize(hist);
      }
    }


  private:

    BinnedHistogram<double> _h_chi_dijet;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1090423);

}
#line 1 "CMS_2012_I1102908.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <sstream>

namespace Rivet {

  /// @cond
  inline double _invert(double x) { return (x > 0) ? 1/x : 0; }
  /// @endcond


  /// @brief CMS inclusive and exclusive dijet production ratio at large rapidity intervals
  class CMS_2012_I1102908 : public Analysis {
  public:

    CMS_2012_I1102908()
      : Analysis("CMS_2012_I1102908")
    {  }


  void init() {
    // Projections
    declare(FastJets(FinalState(), FastJets::ANTIKT, 0.5), "antikT");

    // Histograms
    /// @todo Can we manage to only register these as they are "really" created in the finalize()?
    _h_dijet_ratio = bookScatter2D(1, 1, 1);
    _h_MN_dijet_ratio = bookScatter2D(2, 1, 1);

    // Temporary histograms (directly instantiated)
    _h_DeltaY_exclusive = bookHisto1D("TMP/excl",refData(1, 1, 1));
    _h_DeltaY_inclusive = bookHisto1D("TMP/incl",refData(1, 1, 1));
    _h_DeltaY_MN = bookHisto1D("TMP/YMN",refData(1, 1, 1));
  }


  void analyze(const Event & event) {
    const double weight = event.weight();

    // Jets with  pT > 35.0, -4.7 < y < 4.7
    const JetAlg& jet_alg = apply<JetAlg>(event, "antikT");
    const Jets& jets = jet_alg.jets(Cuts::pT > 35*GeV && Cuts::absrap < 4.7);

    // Veto event if number of jets less than 2
    if (jets.size() < 2) return;

    // Loop over jet pairs
    double deltaY_MN = 0.0;
    for (size_t ij1 = 0; ij1 < jets.size(); ++ij1) {
      for (size_t ij2 = ij1 + 1; ij2 < jets.size(); ++ij2) {
        const double deltaY = fabs(jets[ij1].rapidity() - jets[ij2].rapidity());
        // Exclusive dijet case:
        if (jets.size() == 2) _h_DeltaY_exclusive->fill(deltaY, weight);
        // Inclusive jets case:
        _h_DeltaY_inclusive->fill(deltaY, weight);
        // Mueller-Navelet:
        if (deltaY > deltaY_MN) deltaY_MN = deltaY;
      }
    }
    _h_DeltaY_MN->fill(deltaY_MN, weight);
  }



  void finalize() {
    *_h_dijet_ratio    = YODA::efficiency(*_h_DeltaY_exclusive, *_h_DeltaY_inclusive);
    *_h_MN_dijet_ratio = YODA::efficiency(*_h_DeltaY_exclusive, *_h_DeltaY_MN);
    transformY(*_h_dijet_ratio, _invert);
    transformY(*_h_MN_dijet_ratio, _invert);
  }


  private:

    /// @name Histograms
    //@{
    Scatter2DPtr _h_dijet_ratio, _h_MN_dijet_ratio;
    Histo1DPtr _h_DeltaY_inclusive, _h_DeltaY_exclusive, _h_DeltaY_MN;
    //@}

  };


  DECLARE_RIVET_PLUGIN(CMS_2012_I1102908);

}
#line 1 "CMS_2012_I1107658.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Underlying event activity in the Drell-Yan process at 7 TeV
  class CMS_2012_I1107658 : public Analysis {
  public:

    /// Constructor
    CMS_2012_I1107658()
      : Analysis("CMS_2012_I1107658")
    {   }


    /// Initialization
    void init() {

      /// @note Using a bare muon Z (but with a clustering radius!?)
      Cut cut = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      ZFinder zfinder(FinalState(), cut, PID::MUON, 4*GeV, 140*GeV, 0.2, ZFinder::NOCLUSTER);
      declare(zfinder, "ZFinder");

      ChargedFinalState cfs(-2, 2, 500*MeV);
      VetoedFinalState nonmuons(cfs);
      nonmuons.addVetoPairId(PID::MUON);
      declare(nonmuons, "nonmuons");

      _h_Nchg_towards_pTmumu                 = bookProfile1D(1, 1, 1);
      _h_Nchg_transverse_pTmumu              = bookProfile1D(2, 1, 1);
      _h_Nchg_away_pTmumu                    = bookProfile1D(3, 1, 1);
      _h_pTsum_towards_pTmumu                = bookProfile1D(4, 1, 1);
      _h_pTsum_transverse_pTmumu             = bookProfile1D(5, 1, 1);
      _h_pTsum_away_pTmumu                   = bookProfile1D(6, 1, 1);
      _h_avgpT_towards_pTmumu                = bookProfile1D(7, 1, 1);
      _h_avgpT_transverse_pTmumu             = bookProfile1D(8, 1, 1);
      _h_avgpT_away_pTmumu                   = bookProfile1D(9, 1, 1);
      _h_Nchg_towards_plus_transverse_Mmumu  = bookProfile1D(10, 1, 1);
      _h_pTsum_towards_plus_transverse_Mmumu = bookProfile1D(11, 1, 1);
      _h_avgpT_towards_plus_transverse_Mmumu = bookProfile1D(12, 1, 1);
      _h_Nchg_towards_zmass_81_101           = bookHisto1D(13, 1, 1);
      _h_Nchg_transverse_zmass_81_101        = bookHisto1D(14, 1, 1);
      _h_Nchg_away_zmass_81_101              = bookHisto1D(15, 1, 1);
      _h_pT_towards_zmass_81_101             = bookHisto1D(16, 1, 1);
      _h_pT_transverse_zmass_81_101          = bookHisto1D(17, 1, 1);
      _h_pT_away_zmass_81_101                = bookHisto1D(18, 1, 1);
      _h_Nchg_transverse_zpt_5               = bookHisto1D(19, 1, 1);
      _h_pT_transverse_zpt_5                 = bookHisto1D(20, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");

      if (zfinder.bosons().size() != 1) vetoEvent;

      double Zpt = zfinder.bosons()[0].pT()/GeV;
      double Zphi = zfinder.bosons()[0].phi();
      double Zmass = zfinder.bosons()[0].mass()/GeV;

      Particles particles = apply<VetoedFinalState>(event, "nonmuons").particles();

      int nTowards = 0;
      int nTransverse = 0;
      int nAway = 0;
      double ptSumTowards = 0;
      double ptSumTransverse = 0;
      double ptSumAway = 0;

      foreach (const Particle& p, particles) {
        double dphi = fabs(deltaPhi(Zphi, p.phi()));
        double pT = p.pT();

        if ( dphi < M_PI/3 ) {
          nTowards++;
          ptSumTowards += pT;
          if (Zmass > 81. && Zmass < 101.) _h_pT_towards_zmass_81_101->fill(pT, weight);
        } else if ( dphi < 2.*M_PI/3 ) {
          nTransverse++;
          ptSumTransverse += pT;
          if (Zmass > 81. && Zmass < 101.) _h_pT_transverse_zmass_81_101->fill(pT, weight);
          if (Zpt < 5.) _h_pT_transverse_zpt_5->fill(pT, weight);
        } else {
          nAway++;
          ptSumAway += pT;
          if (Zmass > 81. && Zmass < 101.) _h_pT_away_zmass_81_101->fill(pT, weight);
        }

      } // Loop over particles


      const double area = 8./3.*M_PI;
      if (Zmass > 81. && Zmass < 101.) {
        _h_Nchg_towards_pTmumu->         fill(Zpt, 1./area * nTowards, weight);
        _h_Nchg_transverse_pTmumu->      fill(Zpt, 1./area * nTransverse, weight);
        _h_Nchg_away_pTmumu->            fill(Zpt, 1./area * nAway, weight);
        _h_pTsum_towards_pTmumu->        fill(Zpt, 1./area * ptSumTowards, weight);
        _h_pTsum_transverse_pTmumu->     fill(Zpt, 1./area * ptSumTransverse, weight);
        _h_pTsum_away_pTmumu->           fill(Zpt, 1./area * ptSumAway, weight);
        if (nTowards > 0)    _h_avgpT_towards_pTmumu->    fill(Zpt, ptSumTowards/nTowards, weight);
        if (nTransverse > 0) _h_avgpT_transverse_pTmumu-> fill(Zpt, ptSumTransverse/nTransverse, weight);
        if (nAway > 0)       _h_avgpT_away_pTmumu->       fill(Zpt, ptSumAway/nAway, weight);
        _h_Nchg_towards_zmass_81_101->   fill(nTowards, weight);
        _h_Nchg_transverse_zmass_81_101->fill(nTransverse, weight);
        _h_Nchg_away_zmass_81_101->      fill(nAway, weight);
      }

      if (Zpt < 5.) {
        _h_Nchg_towards_plus_transverse_Mmumu->fill(Zmass, (nTowards + nTransverse)/(2.*area), weight);
        _h_pTsum_towards_plus_transverse_Mmumu->fill(Zmass, (ptSumTowards + ptSumTransverse)/(2.*area), weight);
        if ((nTowards + nTransverse) > 0) _h_avgpT_towards_plus_transverse_Mmumu->fill(Zmass, (ptSumTowards + ptSumTransverse)/(nTowards + nTransverse), weight);
        _h_Nchg_transverse_zpt_5->fill(nTransverse, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pT_towards_zmass_81_101,    safediv(1, _h_Nchg_towards_zmass_81_101->integral(), 0));
      scale(_h_pT_transverse_zmass_81_101, safediv(1, _h_Nchg_transverse_zmass_81_101->integral(), 0));
      scale(_h_pT_away_zmass_81_101,       safediv(1, _h_Nchg_away_zmass_81_101->integral(), 0));
      scale(_h_pT_transverse_zpt_5,        safediv(1, _h_Nchg_transverse_zpt_5->integral(), 0));
      normalize(_h_Nchg_towards_zmass_81_101);
      normalize(_h_Nchg_transverse_zmass_81_101);
      normalize(_h_Nchg_away_zmass_81_101);
      normalize(_h_Nchg_transverse_zpt_5);
    }


  private:


    /// @name Histogram objects
    //@{
    Profile1DPtr _h_Nchg_towards_pTmumu;
    Profile1DPtr _h_Nchg_transverse_pTmumu;
    Profile1DPtr _h_Nchg_away_pTmumu;
    Profile1DPtr _h_pTsum_towards_pTmumu;
    Profile1DPtr _h_pTsum_transverse_pTmumu;
    Profile1DPtr _h_pTsum_away_pTmumu;
    Profile1DPtr _h_avgpT_towards_pTmumu;
    Profile1DPtr _h_avgpT_transverse_pTmumu;
    Profile1DPtr _h_avgpT_away_pTmumu;
    Profile1DPtr _h_Nchg_towards_plus_transverse_Mmumu;
    Profile1DPtr _h_pTsum_towards_plus_transverse_Mmumu;
    Profile1DPtr _h_avgpT_towards_plus_transverse_Mmumu;
    Histo1DPtr _h_Nchg_towards_zmass_81_101;
    Histo1DPtr _h_Nchg_transverse_zmass_81_101;
    Histo1DPtr _h_Nchg_away_zmass_81_101;
    Histo1DPtr _h_pT_towards_zmass_81_101;
    Histo1DPtr _h_pT_transverse_zmass_81_101;
    Histo1DPtr _h_pT_away_zmass_81_101;
    Histo1DPtr _h_Nchg_transverse_zpt_5;
    Histo1DPtr _h_pT_transverse_zpt_5;
    //@}

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1107658);

}
#line 1 "CMS_2012_I1184941.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class CMS_2012_I1184941 : public Analysis {
  public:

    CMS_2012_I1184941()
      : Analysis("CMS_2012_I1184941")
    {   }


    void init() {
      FinalState fs;
      declare(fs, "FS");

      const FastJets jets(FinalState(-4.9, 4.9, 0.0*GeV), FastJets::ANTIKT, 0.5);
      declare(jets, "AntiKtJets05");

      _h_xi = bookHisto1D(1, 1, 1);
    }


    void analyze(const Event& event) {
      double xiM = 0.;
      double xiP = 0.;

      const Jets jets = apply<FastJets>(event, "AntiKtJets05").jetsByPt(20.*GeV);
      if (jets.size() < 2) vetoEvent;  // require a dijet system with a 20 GeV cut on both jets
      if (fabs(jets[0].eta()) > 4.4 || fabs(jets[1].eta()) > 4.4) vetoEvent;

      const FinalState& fsp = apply<FinalState>(event, "FS");

      foreach (const Particle& p, fsp.particles(cmpMomByEta)) {
        const double eta = p.eta();
        const double energy = p.E();
        const double costheta = cos(p.theta());
        // Yes, they really correct to +/- infinity, using Pythia 8 ...
        if (eta < 4.9)  xiP += (energy + energy*costheta);
        if (eta > -4.9 ) xiM += (energy - energy*costheta);
      }

      xiP = xiP / (sqrtS()/GeV);
      xiM = xiM / (sqrtS()/GeV);

      const double weight = event.weight();
      _h_xi->fill( xiM, weight ); // Fill the histogram both with xiP and xiM, and get the average in the endjob.
      _h_xi->fill( xiP, weight );
    }


    void finalize() {
      scale( _h_xi, crossSection()/microbarn/sumOfWeights() / 2.);
    }


  private:

    Histo1DPtr _h_xi;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1184941);

}
#line 1 "CMS_2012_I1193338.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class CMS_2012_I1193338 : public Analysis {
  public:

    CMS_2012_I1193338()
      : Analysis("CMS_2012_I1193338")
    {    }


    void init() {
      declare(ChargedFinalState(-2.4, 2.4, 0.2*GeV), "CFS");
      declare(FinalState(), "FS");

      _h_sigma = bookHisto1D(1, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() > 1) {_h_sigma->fill(1.5, weight);}
      if (cfs.size() > 2) {_h_sigma->fill(2.5, weight);}
      if (cfs.size() > 3) {_h_sigma->fill(3.5, weight);}

      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.size() < 2) vetoEvent; // need at least two particles to calculate gaps

      double gapcenter = 0.;
      double LRG = 0.;
      double etapre = 0.;
      bool first = true;

      foreach(const Particle& p, fs.particles(cmpMomByEta)) { // sorted from minus to plus
        if (first) { // First particle
          first = false;
          etapre = p.eta();
        } else {
          double gap = fabs(p.eta()-etapre);
          if (gap > LRG) {
            LRG = gap; // largest gap
            gapcenter = (p.eta()+etapre)/2.; // find the center of the gap to separate the X and Y systems.
          }
          etapre = p.eta();
        }
      }


      FourMomentum mxFourVector, myFourVector;
      foreach(const Particle& p, fs.particles(cmpMomByEta)) {
        ((p.eta() > gapcenter) ? mxFourVector : myFourVector) += p.momentum();
      }
      const double M2 = max(mxFourVector.mass2(), myFourVector.mass2());
      const double xi = M2/sqr(sqrtS()); // sqrt(s)=7000 GeV, note that units cancel
      if (xi < 5e-6) vetoEvent;

      _h_sigma->fill(0.5, weight);
    }


    void finalize() {
      scale(_h_sigma, crossSection()/millibarn/sumOfWeights());
    }

  private:

    Histo1DPtr _h_sigma;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1193338);

}
#line 1 "CMS_2013_I1122847.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CMS_2013_I1122847 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1122847()
      : Analysis("CMS_2013_I1122847")  {}


    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;

      Cut cuts_mu = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      ZFinder zfinder_mu(fs, cuts_mu, PID::MUON, 40.0*GeV, MAXDOUBLE,
                         0.0, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfinder_mu, "zfinder_mu");

      Cut cuts_el = (Cuts::pT >= 20*GeV && Cuts::abseta < 2.4 && !Cuts::absetaIn(1.447, 1.57));
      ZFinder zfinder_el(fs, cuts_el, PID::ELECTRON, 40.0*GeV, MAXDOUBLE,
                         0.0, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfinder_el, "zfinder_el");


      /// Histograms
      // dimuon
      _hist_mm_100_num = Histo1D(refData(1, 1, 1));
      _hist_mm_125_num = Histo1D(refData(1, 1, 2));
      _hist_mm_150_num = Histo1D(refData(1, 1, 3));
      _hist_mm_240_num = Histo1D(refData(1, 1, 4));

      _hist_mm_100_den = Histo1D(refData(1, 1, 1));
      _hist_mm_125_den = Histo1D(refData(1, 1, 2));
      _hist_mm_150_den = Histo1D(refData(1, 1, 3));
      _hist_mm_240_den = Histo1D(refData(1, 1, 4));

      // Dielectron
      _hist_ee_100_num = Histo1D(refData(2, 1, 1));
      _hist_ee_125_num = Histo1D(refData(2, 1, 2));
      _hist_ee_150_num = Histo1D(refData(2, 1, 3));
      _hist_ee_240_num = Histo1D(refData(2, 1, 4));

      _hist_ee_100_den = Histo1D(refData(2, 1, 1));
      _hist_ee_125_den = Histo1D(refData(2, 1, 2));
      _hist_ee_150_den = Histo1D(refData(2, 1, 3));
      _hist_ee_240_den = Histo1D(refData(2, 1, 4));

      // Dilepton
      _hist_ll_100_num = Histo1D(refData(3, 1, 1));
      _hist_ll_125_num = Histo1D(refData(3, 1, 2));
      _hist_ll_150_num = Histo1D(refData(3, 1, 3));
      _hist_ll_240_num = Histo1D(refData(3, 1, 4));

      _hist_ll_100_den = Histo1D(refData(3, 1, 1));
      _hist_ll_125_den = Histo1D(refData(3, 1, 2));
      _hist_ll_150_den = Histo1D(refData(3, 1, 3));
      _hist_ll_240_den = Histo1D(refData(3, 1, 4));
    }


    double cosThetaCS(const Particle& l1, const Particle& l2) {
      const FourMomentum mom1 = l1.mom();
      const FourMomentum mom2 = l2.mom();
      const FourMomentum mom12 = mom1 + mom2;
      const double Q = mom12.mass();
      const double QT = mom12.pT();
      const double QZ = mom12.pz();

      /// @todo Why include factors of sqrt2 which then get immediately multiplied then divided out?
      const double sqrt2 = sqrt(2.0);
      /// @todo Can be done more nicely via PID-ordered references to mom1, mom2
      const double P1p = ((l1.pid() > 0) ? (mom1.E() + mom1.pz()) : (mom2.E() + mom2.pz())) / sqrt2;
      const double P1m = ((l1.pid() > 0) ? (mom1.E() - mom1.pz()) : (mom2.E() - mom2.pz())) / sqrt2;
      const double P2p = ((l1.pid() > 0) ? (mom2.E() + mom2.pz()) : (mom1.E() + mom1.pz())) / sqrt2;
      const double P2m = ((l1.pid() > 0) ? (mom2.E() - mom2.pz()) : (mom1.E() - mom1.pz())) / sqrt2;

      const double cosThetaCS = sign(QZ) * (2 / (Q * add_quad(Q, QT))) * (P1p*P2m - P1m*P2p);
      return cosThetaCS;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const ZFinder& zfinder_el = apply<ZFinder>(event, "zfinder_el");
      if (zfinder_el.bosons().size() > 0) {
        const Particle& z  = zfinder_el.bosons()[0];
        const Particle& l1 = zfinder_el.constituents()[0];
        const Particle& l2 = zfinder_el.constituents()[1];

        // Prepare variables for filling
        const double rap = z.absrap();
        const double costhetacs = cosThetaCS(l1, l2);
        const double sgn = sign(costhetacs);

        // Fill the histograms
        if (rap < 1.0) {
          _hist_ee_100_num.fill(z.mass(), weight * sgn);
          _hist_ll_100_num.fill(z.mass(), weight * sgn);
          _hist_ee_100_den.fill(z.mass(), weight);
          _hist_ll_100_den.fill(z.mass(), weight);
        } else if (rap < 1.25) {
          _hist_ee_125_num.fill(z.mass(), weight * sgn);
          _hist_ll_125_num.fill(z.mass(), weight * sgn);
          _hist_ee_125_den.fill(z.mass(), weight);
          _hist_ll_125_den.fill(z.mass(), weight);
        } else if (rap < 1.50) {
          _hist_ee_150_num.fill(z.mass(), weight * sgn);
          _hist_ll_150_num.fill(z.mass(), weight * sgn);
          _hist_ee_150_den.fill(z.mass(), weight);
          _hist_ll_150_den.fill(z.mass(), weight);
        } else if (rap < 2.40) {
          _hist_ee_240_num.fill(z.mass(), weight * sgn);
          _hist_ll_240_num.fill(z.mass(), weight * sgn);
          _hist_ee_240_den.fill(z.mass(), weight);
          _hist_ll_240_den.fill(z.mass(), weight);
        }
      }

      const ZFinder& zfinder_mu = apply<ZFinder>(event, "zfinder_mu");
      if (zfinder_mu.bosons().size() > 0) {
        const Particle& z  = zfinder_mu.bosons()[0];
        const Particle& l1 = zfinder_mu.constituents()[0];
        const Particle& l2 = zfinder_mu.constituents()[1];

        // Prepare variables for filling
        const double rap = z.absrap();
        const double costhetacs = cosThetaCS(l1, l2);
        const double sgn = sign(costhetacs);

        // Fill the histograms
        if (rap < 1.0) {
          _hist_mm_100_num.fill(z.mass(), weight * sgn);
          _hist_ll_100_num.fill(z.mass(), weight * sgn);
          _hist_mm_100_den.fill(z.mass(), weight);
          _hist_ll_100_den.fill(z.mass(), weight);
        } else if (rap < 1.25) {
          _hist_mm_125_num.fill(z.mass(), weight * sgn);
          _hist_ll_125_num.fill(z.mass(), weight * sgn);
          _hist_mm_125_den.fill(z.mass(), weight);
          _hist_ll_125_den.fill(z.mass(), weight);
        } else if (rap < 1.50) {
          _hist_mm_150_num.fill(z.mass(), weight * sgn);
          _hist_ll_150_num.fill(z.mass(), weight * sgn);
          _hist_mm_150_den.fill(z.mass(), weight);
          _hist_ll_150_den.fill(z.mass(), weight);
        } else if (rap < 2.40) {
          _hist_mm_240_num.fill(z.mass(), weight * sgn);
          _hist_ll_240_num.fill(z.mass(), weight * sgn);
          _hist_mm_240_den.fill(z.mass(), weight);
          _hist_ll_240_den.fill(z.mass(), weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      divide(_hist_mm_100_num, _hist_mm_100_den, bookScatter2D(1, 1, 1));
      divide(_hist_mm_125_num, _hist_mm_125_den, bookScatter2D(1, 1, 2));
      divide(_hist_mm_150_num, _hist_mm_150_den, bookScatter2D(1, 1, 3));
      divide(_hist_mm_240_num, _hist_mm_240_den, bookScatter2D(1, 1, 4));

      divide(_hist_ee_100_num, _hist_ee_100_den, bookScatter2D(2, 1, 1));
      divide(_hist_ee_125_num, _hist_ee_125_den, bookScatter2D(2, 1, 2));
      divide(_hist_ee_150_num, _hist_ee_150_den, bookScatter2D(2, 1, 3));
      divide(_hist_ee_240_num, _hist_ee_240_den, bookScatter2D(2, 1, 4));

      divide(_hist_ll_100_num, _hist_ll_100_den, bookScatter2D(3, 1, 1));
      divide(_hist_ll_125_num, _hist_ll_125_den, bookScatter2D(3, 1, 2));
      divide(_hist_ll_150_num, _hist_ll_150_den, bookScatter2D(3, 1, 3));
      divide(_hist_ll_240_num, _hist_ll_240_den, bookScatter2D(3, 1, 4));
    }


  private:

    /// Histograms
    Histo1D _hist_ee_100_num, _hist_ee_125_num, _hist_ee_150_num, _hist_ee_240_num;
    Histo1D _hist_ee_100_den, _hist_ee_125_den, _hist_ee_150_den, _hist_ee_240_den;
    Histo1D _hist_mm_100_num, _hist_mm_125_num, _hist_mm_150_num, _hist_mm_240_num;
    Histo1D _hist_mm_100_den, _hist_mm_125_den, _hist_mm_150_den, _hist_mm_240_den;
    Histo1D _hist_ll_100_num, _hist_ll_125_num, _hist_ll_150_num, _hist_ll_240_num;
    Histo1D _hist_ll_100_den, _hist_ll_125_den, _hist_ll_150_den, _hist_ll_240_den;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1122847);

}
#line 1 "CMS_2013_I1208923.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  // This analysis is a derived from the class Analysis:
  class CMS_2013_I1208923 : public Analysis {

  public:
    // Constructor
    CMS_2013_I1208923()
      : Analysis("CMS_2013_I1208923") {
      //setNeedsCrossSection(true);
    }

    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;
      declare(fs, "FS");

      // Initialize the projections
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "Jets");

      // Book histograms
      _h_sigma.addHistogram(0.0, 0.5, bookHisto1D(1, 1, 1));
      _h_sigma.addHistogram(0.5, 1.0, bookHisto1D(1, 1, 2));
      _h_sigma.addHistogram(1.0, 1.5, bookHisto1D(1, 1, 3));
      _h_sigma.addHistogram(1.5, 2.0, bookHisto1D(1, 1, 4));
      _h_sigma.addHistogram(2.0, 2.5, bookHisto1D(1, 1, 5));
      
      _h_invMass.addHistogram(0.0, 0.5, bookHisto1D(2, 1, 1));
      _h_invMass.addHistogram(0.5, 1.0, bookHisto1D(2, 1, 2));
      _h_invMass.addHistogram(1.0, 1.5, bookHisto1D(2, 1, 3));
      _h_invMass.addHistogram(1.5, 2.0, bookHisto1D(2, 1, 4));
      _h_invMass.addHistogram(2.0, 2.5, bookHisto1D(2, 1, 5));
    }

    // Analysis
    void analyze(const Event &event) {
      const double weight = event.weight();
      const FastJets &fJets = apply<FastJets>(event, "Jets");
      
      // Fill the jet pT spectra
      const Jets& jets = fJets.jetsByPt(Cuts::pt>100.*GeV && Cuts::absrap <2.5);
      foreach (const Jet &j, jets) {
        _h_sigma.fill(fabs(j.momentum().rapidity()), j.momentum().pT() / GeV, weight);
      }

      // Require two jets
      const Jets& dijets = fJets.jetsByPt(Cuts::pt>30.*GeV && Cuts::absrap < 2.5);
      if (dijets.size() > 1) {
        if (dijets[0].momentum().pT() / GeV > 60.) {
          // Fill the invariant mass histogram
          double ymax = max(dijets[0].momentum().absrapidity(), dijets[1].momentum().absrapidity());
          double invMass = FourMomentum(dijets[0].momentum() + dijets[1].momentum()).mass();
          _h_invMass.fill(fabs(ymax), invMass, weight);
        }
      } 

    }


    // Scale histograms by the production cross section
    void finalize() {
      _h_sigma.scale(  crossSection() / sumOfWeights() / 2.0, this);
      _h_invMass.scale(crossSection() / sumOfWeights() / 2.0, this);
    }

  private:
    BinnedHistogram<double> _h_sigma;
    BinnedHistogram<double> _h_invMass;
  };

  // This global object acts as a hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_2013_I1208923);
}
#line 1 "CMS_2013_I1209721.cc"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {

  


  /// CMS Z+jets delta(phi) and jet thrust measurement at 7 TeV
  class CMS_2013_I1209721 : public Analysis {
  public:

    CMS_2013_I1209721()
      : Analysis("CMS_2013_I1209721")
    {    }


    /// Book projections and histograms
    void init() {
      // Full final state
      const FinalState fs(-5.0,5.0);
      declare(fs, "FS");
      // Z finders for electrons and muons
      Cut cuts = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      const ZFinder zfe(fs, cuts, PID::ELECTRON, 71*GeV, 111*GeV);
      const ZFinder zfm(fs, cuts, PID::MUON,     71*GeV, 111*GeV);
      declare(zfe, "ZFE");
      declare(zfm, "ZFM");
      // Jets
      const FastJets jets(fs, FastJets::ANTIKT, 0.5);
      declare(jets, "JETS");

      // Book histograms from data
      for (size_t i = 0; i < 2; ++i) {
        _histDeltaPhiZJ1_1[i]  = bookHisto1D(1+i*9, 1, 1);
        _histDeltaPhiZJ1_2[i]  = bookHisto1D(2+i*9, 1, 1);
        _histDeltaPhiZJ1_3[i]  = bookHisto1D(4+i*9, 1, 1);
        _histDeltaPhiZJ2_3[i]  = bookHisto1D(5+i*9, 1, 1);
        _histDeltaPhiZJ3_3[i]  = bookHisto1D(3+i*9, 1, 1);
        _histDeltaPhiJ1J2_3[i] = bookHisto1D(6+i*9, 1, 1);
        _histDeltaPhiJ1J3_3[i] = bookHisto1D(7+i*9, 1, 1);
        _histDeltaPhiJ2J3_3[i] = bookHisto1D(8+i*9, 1, 1);
        _histTransvThrust[i]   = bookHisto1D(9+i*9, 1, 1);
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      // Apply the Z finders
      const ZFinder& zfe = apply<ZFinder>(event, "ZFE");
      const ZFinder& zfm = apply<ZFinder>(event, "ZFM");

      // Choose the Z candidate (there must be one)
      if (zfe.empty() && zfm.empty()) vetoEvent;
      const ParticleVector& z = !zfm.empty() ? zfm.bosons() : zfe.bosons();
      const ParticleVector& leptons = !zfm.empty() ? zfm.constituents() : zfe.constituents();

      // Determine whether we are in the boosted regime
      const bool is_boosted = (z[0].pT() > 150*GeV);

      // Build the jets
      const FastJets& jetfs = apply<FastJets>(event, "JETS");
      const Jets& jets = jetfs.jetsByPt(Cuts::pT > 50*GeV && Cuts::abseta < 2.5);

      // Clean the jets against the lepton candidates, as in the paper, with a deltaR cut of 0.4 against the clustered leptons
      vector<const Jet*> cleanedJets;
      for (size_t i = 0; i < jets.size(); ++i) {
        bool isolated = true;
        for (size_t j = 0; j < 2; ++j) {
          if (deltaR(leptons[j], jets[i]) < 0.4) {
            isolated = false;
            break;
          }
        }
        if (isolated) cleanedJets.push_back(&jets[i]);
      }

      // Require at least 1 jet
      const unsigned int Njets = cleanedJets.size();
      if (Njets < 1) vetoEvent;

      // Now compute the thrust
      // Collect Z and jets transverse momenta to calculate transverse thrust
      vector<Vector3> momenta;
      momenta.clear();
      Vector3 mom = z[0].p3();
      mom.setZ(0);
      momenta.push_back(mom);

      for (size_t i = 0; i < cleanedJets.size(); ++i) {
        Vector3 mj = cleanedJets[i]->momentum().p3();
        mj.setZ(0);
        momenta.push_back(mj);
      }

      if (momenta.size() <= 2){
        // We need to use a ghost so that Thrust.calc() doesn't return 1.
        momenta.push_back(Vector3(0.0000001,0.0000001,0.));
      }

      // Define a macro to appropriately fill both unboosted and boosted histo versions
      #define FILLx2(HNAME, VAL) do { double x = VAL; for (size_t i = 0; i < 2; ++i) { \
        if (i == 0 || is_boosted) HNAME[i]->fill(x, weight); } } while(0)

      Thrust thrust; thrust.calc(momenta);
      const double T = thrust.thrust();
      FILLx2(_histTransvThrust, log(max(1-T, 1e-6)));

      const double dphiZJ1 = deltaPhi(z[0], *cleanedJets[0]);
      FILLx2(_histDeltaPhiZJ1_1, dphiZJ1);
      if (Njets > 1) {
        FILLx2(_histDeltaPhiZJ1_2, dphiZJ1);
        if (Njets > 2) {
          FILLx2(_histDeltaPhiZJ1_3, dphiZJ1);
          FILLx2(_histDeltaPhiZJ2_3, deltaPhi(z[0], *cleanedJets[1]));
          FILLx2(_histDeltaPhiZJ3_3, deltaPhi(z[0], *cleanedJets[2]));
          FILLx2(_histDeltaPhiJ1J2_3, deltaPhi(*cleanedJets[0], *cleanedJets[1]));
          FILLx2(_histDeltaPhiJ1J3_3, deltaPhi(*cleanedJets[0], *cleanedJets[2]));
          FILLx2(_histDeltaPhiJ2J3_3, deltaPhi(*cleanedJets[1], *cleanedJets[2]));
        }
      }
    }


    /// Normalizations
    /// @note Most of these data normalizations neglect the overflow bins
    void finalize() {
      for (size_t i = 0; i < 2; ++i) {
        normalize(_histDeltaPhiZJ1_1[i], 1, false);
        normalize(_histDeltaPhiZJ1_2[i], 1, false);
        normalize(_histDeltaPhiZJ1_3[i], 1, false);
        normalize(_histDeltaPhiZJ2_3[i], 1, false);
        normalize(_histDeltaPhiZJ3_3[i], 1, false);
        normalize(_histDeltaPhiJ1J2_3[i], 1, false);
        normalize(_histDeltaPhiJ1J3_3[i], 1, false);
        normalize(_histDeltaPhiJ2J3_3[i], 1, false);
        normalize(_histTransvThrust[i]);
      }
    }


  private:

    // Arrays of unboosted/boosted histos
    Histo1DPtr _histDeltaPhiZJ1_1[2];
    Histo1DPtr _histDeltaPhiZJ1_2[2];
    Histo1DPtr _histDeltaPhiZJ1_3[2];
    Histo1DPtr _histDeltaPhiZJ2_3[2];
    Histo1DPtr _histDeltaPhiZJ3_3[2];
    Histo1DPtr _histDeltaPhiJ1J2_3[2];
    Histo1DPtr _histDeltaPhiJ1J3_3[2];
    Histo1DPtr _histDeltaPhiJ2J3_3[2];
    Histo1DPtr _histTransvThrust[2];

  };


  DECLARE_RIVET_PLUGIN(CMS_2013_I1209721);

}
#line 1 "CMS_2013_I1218372.cc"
// Samantha Dooling DESY
// February 2012
//
// -*- C++ -*-
// =============================
//
// Ratio of the energy deposited in the pseudorapidity range
// -6.6 < eta < -5.2 for events with a charged particle jet
//
// =============================
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  class CMS_2013_I1218372 : public Analysis {
  public:

  /// Constructor
  CMS_2013_I1218372()
    : Analysis("CMS_2013_I1218372")
    { }

    void init() {

      // gives the range of eta and min pT for the final state from which I get the jets
      FastJets jetpro (ChargedFinalState(-2.5, 2.5, 0.3*GeV), FastJets::ANTIKT, 0.5);
      declare(jetpro, "Jets");

      // skip Neutrinos and Muons
      VetoedFinalState fsv(FinalState(-7.0, -4.0, 0.*GeV));
      fsv.vetoNeutrinos();
      fsv.addVetoPairId(PID::MUON);
      declare(fsv, "fsv");

      // for the hadron level selection
      VetoedFinalState sfsv(FinalState(-MAXDOUBLE, MAXDOUBLE, 0.*GeV));
      sfsv.vetoNeutrinos();
      sfsv.addVetoPairId(PID::MUON);
      declare(sfsv, "sfsv");

      //counters
      passedSumOfWeights = 0.;
      inclEflow = 0.;

      // Temporary histograms to fill the energy flow for leading jet events.
      // Ratios are calculated in finalyze().
      int id = 0;
      if (fuzzyEquals(sqrtS()/GeV,  900, 1e-3)) id=1;
      if (fuzzyEquals(sqrtS()/GeV, 2760, 1e-3)) id=2;
      if (fuzzyEquals(sqrtS()/GeV, 7000, 1e-3)) id=3;
      _h_ratio  = bookScatter2D(id, 1, 1);
      _tmp_jet  = bookHisto1D  ("TMP/eflow_jet"  ,refData(id, 1, 1));  // Leading jet energy flow in pt
      _tmp_njet = bookHisto1D  ("TMP/number_jet" ,refData(id, 1, 1)); // Number of events in pt
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // Skip if the event is empty
      const FinalState& fsv = apply<FinalState>(event, "fsv");
      if (fsv.empty()) vetoEvent;

      // ====================== Minimum Bias selection

      const FinalState& sfsv = apply<FinalState>(event, "sfsv");
      Particles parts = sfsv.particles(cmpMomByRap);
      if (parts.empty()) vetoEvent;

      // find dymax
      double dymax = 0;
      int gap_pos  = -1;
      for (size_t i = 0; i < parts.size()-1; ++i) {
        double dy = parts[i+1].rapidity() - parts[i].rapidity();
        if (dy > dymax) {
          dymax = dy;
          gap_pos = i;
        }
      }

      // calculate mx2 and my2
      FourMomentum xmom;
      for (int i=0; i<=gap_pos; ++i) {
        xmom += parts[i].momentum();
      }
      double mx2 = xmom.mass2();
      if (mx2<0) vetoEvent;

      FourMomentum ymom;
      for (size_t i=gap_pos+1; i<parts.size(); ++i) {
        ymom += parts[i].momentum();
      }
      double my2 = ymom.mass2();
      if (my2<0) vetoEvent;

      // calculate xix and xiy and xidd
      double xix  = mx2 / sqr(sqrtS());
      double xiy  = my2 / sqr(sqrtS());
      double xidd = mx2*my2 / sqr(sqrtS()*0.938*GeV);

      // combine the selection: xi cuts
      bool passedHadronCuts = false;
      if (fuzzyEquals(sqrtS()/GeV,  900, 1e-3) && (xix > 0.1  || xiy > 0.4 || xidd > 0.5)) passedHadronCuts = true;
      if (fuzzyEquals(sqrtS()/GeV, 2760, 1e-3) && (xix > 0.07 || xiy > 0.2 || xidd > 0.5)) passedHadronCuts = true;
      if (fuzzyEquals(sqrtS()/GeV, 7000, 1e-3) && (xix > 0.04 || xiy > 0.1 || xidd > 0.5)) passedHadronCuts = true;
      if (!passedHadronCuts) vetoEvent;

      //  ============================== MINIMUM BIAS EVENTS

      // loop over particles to calculate the energy
      passedSumOfWeights += weight;

      foreach (const Particle& p, fsv.particles()) {
        if (-5.2 > p.eta() && p.eta() > -6.6) inclEflow += weight*p.E()/GeV;
      }

      //  ============================== JET EVENTS

      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets& jets = jetpro.jetsByPt(1.0*GeV);
      if (jets.size()<1) vetoEvent;

      if (fabs(jets[0].eta()) < 2.0) {
        _tmp_njet->fill(jets[0].pT()/GeV, weight);

        // energy flow
        foreach (const Particle& p, fsv.particles()) {
          if (p.eta() > -6.6 && p.eta() < -5.2) {  // ask for the CASTOR region
            _tmp_jet->fill(jets[0].pT()/GeV, weight * p.E()/GeV);
          }
        }
      }

    }// analysis

    void finalize() {
      scale(_tmp_jet, passedSumOfWeights/inclEflow);
      divide(_tmp_jet, _tmp_njet, _h_ratio);
    }

  private:
    // counters
    double passedSumOfWeights;
    double inclEflow;

    // histograms
    Scatter2DPtr _h_ratio;
    Histo1DPtr   _tmp_jet;
    Histo1DPtr   _tmp_njet;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1218372);

}
#line 1 "CMS_2013_I1223519.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalStates.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/Smearing.hh"
#include <bitset>

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2013_I1223519 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2013_I1223519);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState calofs(Cuts::abseta < 5.0);
      declare(calofs, "Clusters");

      MissingMomentum mm(calofs);
      declare(mm, "TruthMET");
      declare(SmearedMET(mm, MET_SMEAR_CMS_RUN2), "MET");

      FastJets fj(calofs, FastJets::ANTIKT, 0.5);
      declare(fj, "TruthJets");
      declare(SmearedJets(fj, JET_SMEAR_CMS_RUN2, [](const Jet& j) {
            if (j.abseta() > 2.4) return 0.;
            return j.bTagged() ? 0.65 : 0.01; }), "Jets"); ///< @note Charm mistag and exact b-tag eff not given

      FinalState ys(Cuts::abspid == PID::PHOTON && Cuts::abseta < 5.0);
      declare(ys, "TruthPhotons");
      declare(SmearedParticles(ys, PHOTON_EFF_CMS_RUN2 /*, PHOTON_SMEAR_CMS_RUN2 */), "Photons");

      FinalState es(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.5);
      declare(es, "TruthElectrons");
      declare(SmearedParticles(es, ELECTRON_EFF_CMS_RUN2, ELECTRON_SMEAR_CMS_RUN2), "Electrons");

      FinalState mus(Cuts::abspid == PID::MUON && Cuts::abseta < 2.4);
      declare(mus, "TruthMuons");
      declare(SmearedParticles(mus, MUON_EFF_CMS_RUN2, MUON_SMEAR_CMS_RUN2), "Muons");

      ChargedFinalState cfs(Cuts::abseta < 2.5);
      declare(cfs, "TruthTracks");
      declare(SmearedParticles(cfs, TRK_EFF_CMS_RUN2), "Tracks");


      // Book histograms
      _h_alphaT23 = bookHisto1D("alphaT23", 15, 0, 3);
      _h_alphaT4 = bookHisto1D("alphaT4", 15, 0, 3);
      /// @todo Add HT histograms

      // Book counters
      _h_srcounters.resize(8*7 + 3);
      for (size_t inj = 0; inj < 2; ++inj) {
        const size_t njmax = inj + 3;
        for (size_t nb = 0; nb < njmax; ++nb) {
          for (size_t iht = 0; iht < 8; ++iht) {
            const size_t i = 8 * ((inj == 0 ? 0 : 3) + nb) + iht;
            _h_srcounters[i] = bookCounter("srcount_j" + toString(njmax) + "_b" + toString(nb) + "_ht" + toString(iht+1));
          }
        }
      }
      // Special nj >= 4, nb >= 4 bins
      for (size_t iht = 0; iht < 3; ++iht) {
        _h_srcounters[8*7 + iht] = bookCounter("srcount_j4_b4_ht" + toString(iht+1));
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get baseline photons, electrons & muons
      Particles photons = apply<ParticleFinder>(event, "Photons").particles(Cuts::pT > 25*GeV);
      Particles elecs = apply<ParticleFinder>(event, "Electrons").particles(Cuts::pT > 10*GeV);
      Particles muons = apply<ParticleFinder>(event, "Muons").particles(Cuts::pT > 10*GeV);

      // Electron/muon isolation (guesswork/copied from other CMS analysis -- paper is unspecific)
      const Particles calofs = apply<ParticleFinder>(event, "Clusters").particles();
      ifilter_discard(photons, [&](const Particle& y) {
          double ptsum = -y.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,y) < 0.3) ptsum += p.pT();
          return ptsum / y.pT() > 0.1;
        });
      ifilter_discard(elecs, [&](const Particle& e) {
          double ptsum = -e.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,e) < 0.3) ptsum += p.pT();
          return ptsum / e.pT() > 0.1;
        });
      ifilter_discard(muons, [&](const Particle& m) {
          double ptsum = -m.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,m) < 0.3) ptsum += p.pT();
          return ptsum / m.pT() > 0.2;
        });

      // Veto the event if there are any remaining baseline photons or leptons
      if (!photons.empty()) vetoEvent;
      if (!elecs.empty()) vetoEvent;
      if (!muons.empty()) vetoEvent;


      // Get jets and apply jet-based event-selection cuts
      const JetAlg& jetproj = apply<JetAlg>(event, "Jets");
      const Jets alljets = jetproj.jetsByPt(Cuts::abseta < 3.0 && Cuts::Et > 37*GeV); //< most inclusive jets requirement
      if (filter_select(alljets, Cuts::Et > 73*GeV).size() < 2) vetoEvent; //< most inclusive lead jets requirement

      // Filter jets into different Et requirements & compute corresponding HTs
      /// @note It's not clear if different HTs are used to choose the HT bins
      const Jets jets37 = filter_select(alljets, Cuts::Et > 37*GeV);
      const Jets jets43 = filter_select(jets37, Cuts::Et > 43*GeV);
      const Jets jets50 = filter_select(jets43, Cuts::Et > 50*GeV);
      const double ht37 = sum(jets37, Et, 0.0);
      const double ht43 = sum(jets43, Et, 0.0);
      const double ht50 = sum(jets50, Et, 0.0);

      // Find the relevant HT bin and apply leading jet event-selection cuts
      static const vector<double> htcuts = { /* 275., 325., */ 375., 475., 575., 675., 775., 875.}; //< comment to avoid jets50 "fall-down"
      const int iht = inRange(ht37, 275*GeV, 325*GeV) ? 0 : inRange(ht43, 325*GeV, 375*GeV) ? 1 : (2+binIndex(ht50, htcuts, true));
      MSG_TRACE("HT = {" << ht37 << ", " << ht43 << ", " << ht50 << "} => IHT = " << iht);
      if (iht < 0) vetoEvent;
      if (iht == 1 && filter_select(jets43, Cuts::Et > 78*GeV).size() < 2) vetoEvent;
      if (iht >= 2 && filter_select(jets50, Cuts::Et > 100*GeV).size() < 2) vetoEvent;

      // Create references for uniform access to relevant set of jets & HT
      const double etcut = iht == 0 ? 37. : iht == 1 ? 43. : 50.;
      const double& ht = iht == 0 ? ht37 : iht == 1 ? ht43 : ht50;
      const Jets& jets = iht == 0 ? jets37 : iht == 1 ? jets43 : jets50;
      if (!jetproj.jets(Cuts::abseta > 3 && Cuts::Et > etcut*GeV).empty()) vetoEvent;
      const size_t nj = jets.size();
      const size_t nb = count_if(jets.begin(), jets.end(), [](const Jet& j) { return j.bTagged(Cuts::pT > 5*GeV); });

      // Compute HTmiss = pT of 4-vector sum of jet momenta
      const FourMomentum jsum = sum(jets, mom, FourMomentum());
      const double htmiss = jsum.pT();

      // Require HTmiss / ETmiss < 1.25
      const double etmiss = apply<SmearedMET>(event, "MET").met();
      if (htmiss/etmiss > 1.25) vetoEvent;

      // Compute DeltaHT = minimum difference of "dijet" ETs, i.e. max(|1+2-3|, |1+3-2|, |2+3-1|)
      double deltaht = -1;
      vector<double> jetets; transform(jets, jetets, Et);
      for (int i = 1; i < (1 << (jetets.size()-1)); ++i) { // count from 1 to 2**N-1, i.e. through all heterogeneous bitmasks with MSB(2**N)==0
        const bitset<10> bits(i); /// @warning There'd better not be more than 10 jets...
        const double htdiff = partition_diff(bits, jetets);
        // MSG_INFO(bits.to_string() << " => " << htdiff);
        if (deltaht < 0 || htdiff < deltaht) deltaht = htdiff;
      }
      MSG_DEBUG("dHT_bitmask = " << deltaht);

      // Cross-check calculation in 2- and 3-jet cases
      // if (jets.size() == 2) {
      //   MSG_INFO("dHT2 = " << fabs(jets[0].Et() - jets[1].Et()));
      // } else if (jets.size() == 3) {
      //   double deltaht_01_2 = fabs(jets[0].Et()+jets[1].Et()-jets[2].Et());
      //   double deltaht_02_1 = fabs(jets[0].Et()+jets[2].Et()-jets[1].Et());
      //   double deltaht_12_0 = fabs(jets[1].Et()+jets[2].Et()-jets[0].Et());
      //   MSG_INFO("dHT3 = " << min({deltaht_01_2, deltaht_02_1, deltaht_12_0}));
      // }

      // Compute alphaT from the above
      double alphaT = fabs(0.5*((ht-deltaht)/(sqrt((ht*ht)-(htmiss*htmiss)))));
      if (alphaT < 0.55) vetoEvent;

      /// @todo Need to include trigger efficiency sampling or weighting?

      // Fill histograms
      const double weight = event.weight();
      const size_t inj = nj < 4 ? 0 : 1;
      const size_t inb = nb < 4 ? nb : 4;
      if (iht >= 2)
        (inj == 0 ? _h_alphaT23 : _h_alphaT4)->fill(alphaT, weight);

      // Fill the appropriate counter -- after working out the irregular SR bin index! *sigh*
      size_t i = 8 * ((inj == 0 ? 0 : 3) + inb) + iht;
      if (inj == 1 && inb == 4) i = 8*7 + (iht < 3 ? iht : 2);
      MSG_INFO("inj = " << inj << ", inb = " << inb << ", i = " << i);
      _h_srcounters[i]->fill(weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = crossSection()/femtobarn*11.7/sumOfWeights();
      scale({_h_alphaT23,_h_alphaT4}, sf);
      for (size_t i = 0; i < 8*7+3; ++i)
        scale(_h_srcounters[i], sf);

    }

    //@}


    /// @name Utility functions for partitioning jet pTs into two groups and summing/diffing them
    //@{

    /// Sum the given values into two subsets according to the provided bitmask
    template <size_t N>
    pair<double, double> partition_sum(const bitset<N>& mask, const vector<double>& vals) const {
      pair<double, double> rtn(0., 0.);
      for (size_t i = 0; i < vals.size(); ++i) {
        (!mask[vals.size()-1-i] ? rtn.first : rtn.second) += vals[i];
      }
      return rtn;
    }

    /// Return the difference between summed subsets according to the provided bitmask
    template <size_t N>
    double partition_diff(const bitset<N>& mask, const vector<double>& vals) const {
      const pair<double, double> sums = partition_sum(mask, vals);
      const double diff = fabs(sums.first - sums.second);
      MSG_TRACE(mask.to_string() << ": " << sums.first << "/" << sums.second << " => " << diff);
      return diff;
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_alphaT23, _h_alphaT4;
    vector<CounterPtr> _h_srcounters;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1223519);


}
#line 1 "CMS_2013_I1224539_DIJET.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

namespace Rivet {


  class CMS_2013_I1224539_DIJET : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2013_I1224539_DIJET()
      : Analysis("CMS_2013_I1224539_DIJET"),
        _filter(fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3))),
        _trimmer(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03))),
        _pruner(fastjet::Pruner(fastjet::cambridge_algorithm, 0.1, 0.5))
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(-2.4, 2.4, 0*GeV);
      declare(fs, "FS");

      // Jet collections
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "JetsAK7");
      declare(FastJets(fs, FastJets::CAM, 0.8), "JetsCA8");
      declare(FastJets(fs, FastJets::CAM, 1.2), "JetsCA12");

      // Histograms
      for (size_t i = 0; i < N_PT_BINS_dj; ++i ) {
        _h_ungroomedAvgJetMass_dj[i] = bookHisto1D(i+1+0*N_PT_BINS_dj, 1, 1);
        _h_filteredAvgJetMass_dj[i] = bookHisto1D(i+1+1*N_PT_BINS_dj, 1, 1);
        _h_trimmedAvgJetMass_dj[i] = bookHisto1D(i+1+2*N_PT_BINS_dj, 1, 1);
        _h_prunedAvgJetMass_dj[i] = bookHisto1D(i+1+3*N_PT_BINS_dj, 1, 1);
      }
    }


    // Find the pT histogram bin index for value pt (in GeV), to hack a 2D histogram equivalent
    /// @todo Use a YODA axis/finder alg when available
    size_t findPtBin(double ptJ) {
      const double ptBins_dj[N_PT_BINS_dj+1] = { 220.0, 300.0, 450.0, 500.0, 600.0, 800.0, 1000.0, 1500.0};
      for (size_t ibin = 0; ibin < N_PT_BINS_dj; ++ibin) {
        if (inRange(ptJ, ptBins_dj[ibin], ptBins_dj[ibin+1])) return ibin;
      }
      return N_PT_BINS_dj;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Look at events with >= 2 jets
      const PseudoJets& psjetsAK7 = apply<FastJets>(event, "JetsAK7").pseudoJetsByPt( 50.0*GeV );
      if (psjetsAK7.size() < 2) vetoEvent;

      // Get the leading two jets and find their average pT
      const fastjet::PseudoJet& j0 = psjetsAK7[0];
      const fastjet::PseudoJet& j1 = psjetsAK7[1];
      double ptAvg = 0.5 * (j0.pt() + j1.pt());

      // Find the appropriate mean pT bin and escape if needed
      const size_t njetBin = findPtBin(ptAvg/GeV);
      if (njetBin >= N_PT_BINS_dj) vetoEvent;

      // Now run the substructure algs...
      fastjet::PseudoJet filtered0 = _filter(j0);
      fastjet::PseudoJet filtered1 = _filter(j1);
      fastjet::PseudoJet trimmed0 = _trimmer(j0);
      fastjet::PseudoJet trimmed1 = _trimmer(j1);
      fastjet::PseudoJet pruned0 = _pruner(j0);
      fastjet::PseudoJet pruned1 = _pruner(j1);

      // ... and fill the histograms
      _h_ungroomedAvgJetMass_dj[njetBin]->fill(0.5*(j0.m() + j1.m())/GeV, weight);
      _h_filteredAvgJetMass_dj[njetBin]->fill(0.5*(filtered0.m() + filtered1.m())/GeV, weight);
      _h_trimmedAvgJetMass_dj[njetBin]->fill(0.5*(trimmed0.m() + trimmed1.m())/GeV, weight);
      _h_prunedAvgJetMass_dj[njetBin]->fill(0.5*(pruned0.m() + pruned1.m())/GeV, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double normalizationVal = 1000;
      for (size_t i = 0; i < N_PT_BINS_dj; ++i) {
        normalize(_h_ungroomedAvgJetMass_dj[i], normalizationVal);
        normalize(_h_filteredAvgJetMass_dj[i], normalizationVal);
        normalize(_h_trimmedAvgJetMass_dj[i], normalizationVal);
        normalize(_h_prunedAvgJetMass_dj[i], normalizationVal);
      }
    }

    //@}


  private:

    /// @name FastJet grooming tools (configured in constructor init list)
    //@{
    const fastjet::Filter _filter;
    const fastjet::Filter _trimmer;
    const fastjet::Pruner _pruner;
    //@}


    /// @name Histograms
    //@{
    enum BINS_dj { PT_220_300_dj=0, PT_300_450_dj, PT_450_500_dj, PT_500_600_dj,
                   PT_600_800_dj, PT_800_1000_dj, PT_1000_1500_dj, N_PT_BINS_dj };
    Histo1DPtr _h_ungroomedJet0pt, _h_ungroomedJet1pt;
    Histo1DPtr _h_ungroomedAvgJetMass_dj[N_PT_BINS_dj];
    Histo1DPtr _h_filteredAvgJetMass_dj[N_PT_BINS_dj];
    Histo1DPtr _h_trimmedAvgJetMass_dj[N_PT_BINS_dj];
    Histo1DPtr _h_prunedAvgJetMass_dj[N_PT_BINS_dj];
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1224539_DIJET);


}
#line 1 "CMS_2013_I1224539_WJET.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

namespace Rivet {


  class CMS_2013_I1224539_WJET : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2013_I1224539_WJET()
      : Analysis("CMS_2013_I1224539_WJET"),
        _filter(fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3))),
        _trimmer(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03))),
        _pruner(fastjet::Pruner(fastjet::cambridge_algorithm, 0.1, 0.5))
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(-2.4, 2.4, 0*GeV);
      declare(fs, "FS");

      // Find W's with pT > 120, MET > 50
      WFinder wfinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 80*GeV, PID::ELECTRON, 50*GeV, 1000*GeV, 50.0*GeV,
                      0.2, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      declare(wfinder, "WFinder");

      // W+jet jet collections
      declare(FastJets(wfinder.remainingFinalState(), FastJets::ANTIKT, 0.7), "JetsAK7_wj");
      declare(FastJets(wfinder.remainingFinalState(), FastJets::CAM, 0.8), "JetsCA8_wj");
      declare(FastJets(wfinder.remainingFinalState(), FastJets::CAM, 1.2), "JetsCA12_wj");

      // Histograms
      /// @note These are 2D histos rendered into slices
      const int wjetsOffset = 51;
      for (size_t i = 0; i < N_PT_BINS_vj; ++i) {
        _h_ungroomedJetMass_AK7_wj[i] = bookHisto1D(wjetsOffset+i+1+0*N_PT_BINS_vj, 1, 1);
        _h_filteredJetMass_AK7_wj[i] = bookHisto1D(wjetsOffset+i+1+1*N_PT_BINS_vj, 1, 1);
        _h_trimmedJetMass_AK7_wj[i] = bookHisto1D(wjetsOffset+i+1+2*N_PT_BINS_vj, 1, 1);
        _h_prunedJetMass_AK7_wj[i] = bookHisto1D(wjetsOffset+i+1+3*N_PT_BINS_vj, 1, 1);
        _h_prunedJetMass_CA8_wj[i] = bookHisto1D(wjetsOffset+i+1+4*N_PT_BINS_vj, 1, 1);
        if (i > 0) _h_filteredJetMass_CA12_wj[i] = bookHisto1D(wjetsOffset+i+5*N_PT_BINS_vj, 1, 1);
      }
    }


    bool isBackToBack_wj(const WFinder& wf, const fastjet::PseudoJet& psjet) {
      const FourMomentum& w = wf.bosons()[0].momentum();
      const FourMomentum& l1 = wf.constituentLeptons()[0].momentum();
      const FourMomentum& l2 = wf.constituentNeutrinos()[0].momentum();
      /// @todo We should make FourMomentum know how to construct itself from a PseudoJet
      const FourMomentum jmom(psjet.e(), psjet.px(), psjet.py(), psjet.pz());
      return (deltaPhi(w, jmom) > 2.0 && deltaR(l1, jmom) > 1.0 && deltaPhi(l2, jmom) > 0.4);
    }


    // Find the pT histogram bin index for value pt (in GeV), to hack a 2D histogram equivalent
    /// @todo Use a YODA axis/finder alg when available
    size_t findPtBin(double ptJ) {
      const double ptBins_vj[N_PT_BINS_vj+1] = { 125.0, 150.0, 220.0, 300.0, 450.0 };
      for (size_t ibin = 0; ibin < N_PT_BINS_vj; ++ibin) {
        if (inRange(ptJ, ptBins_vj[ibin], ptBins_vj[ibin+1])) return ibin;
      }
      return N_PT_BINS_vj;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get the W
      const WFinder& wfinder = apply<WFinder>(event, "WFinder");
      if (wfinder.bosons().size() != 1) vetoEvent;
      const Particle& w = wfinder.bosons()[0];
      const Particle& l = wfinder.constituentLeptons()[0];

      // Require a fairly high-pT W and charged lepton
      if (l.pT() < 80*GeV || w.pT() < 120*GeV) vetoEvent;

      // Get the pseudojets.
      const PseudoJets& psjetsCA8_wj = apply<FastJets>(event, "JetsCA8_wj").pseudoJetsByPt( 50.0*GeV );
      const PseudoJets& psjetsCA12_wj = apply<FastJets>(event, "JetsCA12_wj").pseudoJetsByPt( 50.0*GeV );

      // AK7 jets
      const PseudoJets& psjetsAK7_wj = apply<FastJets>(event, "JetsAK7_wj").pseudoJetsByPt( 50.0*GeV );
      if (!psjetsAK7_wj.empty()) {
        // Get the leading jet and make sure it's back-to-back with the W
        const fastjet::PseudoJet& j0 = psjetsAK7_wj[0];
        if (isBackToBack_wj(wfinder, j0)) {
          const size_t njetBin = findPtBin(j0.pt()/GeV);
          if (njetBin < N_PT_BINS_vj) {
            fastjet::PseudoJet filtered0 = _filter(j0);
            fastjet::PseudoJet trimmed0 = _trimmer(j0);
            fastjet::PseudoJet pruned0 = _pruner(j0);
            _h_ungroomedJetMass_AK7_wj[njetBin]->fill(j0.m()/GeV, weight);
            _h_filteredJetMass_AK7_wj[njetBin]->fill(filtered0.m()/GeV, weight);
            _h_trimmedJetMass_AK7_wj[njetBin]->fill(trimmed0.m()/GeV, weight);
            _h_prunedJetMass_AK7_wj[njetBin]->fill(pruned0.m()/GeV, weight);
          }
        }
      }

      // CA8 jets
      if (!psjetsCA8_wj.empty()) {
        // Get the leading jet and make sure it's back-to-back with the W
        const fastjet::PseudoJet& j0 = psjetsCA8_wj[0];
        if (isBackToBack_wj(wfinder, j0)) {
          const size_t njetBin = findPtBin(j0.pt()/GeV);
          if (njetBin < N_PT_BINS_vj) {
            fastjet::PseudoJet pruned0 = _pruner(j0);
            _h_prunedJetMass_CA8_wj[njetBin]->fill(pruned0.m()/GeV, weight);
          }
        }
      }

      // CA12 jets
      if (!psjetsCA12_wj.empty()) {
        // Get the leading jet and make sure it's back-to-back with the W
        const fastjet::PseudoJet& j0 = psjetsCA12_wj[0];
        if (isBackToBack_wj(wfinder, j0)) {
          const size_t njetBin = findPtBin(j0.pt()/GeV);
          if (njetBin < N_PT_BINS_vj&&njetBin>0) {
            fastjet::PseudoJet filtered0 = _filter(j0);
            _h_filteredJetMass_CA12_wj[njetBin]->fill( filtered0.m() / GeV, weight);
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double normalizationVal = 1000;
      for (size_t i = 0; i < N_PT_BINS_vj; ++i) {
        normalize(_h_ungroomedJetMass_AK7_wj[i], normalizationVal);
        normalize(_h_filteredJetMass_AK7_wj[i], normalizationVal);
        normalize(_h_trimmedJetMass_AK7_wj[i], normalizationVal);
        normalize(_h_prunedJetMass_AK7_wj[i], normalizationVal);
        normalize(_h_prunedJetMass_CA8_wj[i], normalizationVal);
        if (i > 0) normalize( _h_filteredJetMass_CA12_wj[i], normalizationVal);
      }
    }

    //@}


  private:

    /// @name FastJet grooming tools (configured in constructor init list)
    //@{
    const fastjet::Filter _filter;
    const fastjet::Filter _trimmer;
    const fastjet::Pruner _pruner;
    //@}


    /// @name Histograms
    //@{
    enum BINS_vj { PT_125_150_vj=0, PT_150_220_vj, PT_220_300_vj, PT_300_450_vj, N_PT_BINS_vj };
    Histo1DPtr _h_ungroomedJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_trimmedJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_CA8_wj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_CA12_wj[N_PT_BINS_vj];
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1224539_WJET);

}
#line 1 "CMS_2013_I1224539_ZJET.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

namespace Rivet {

  


  class CMS_2013_I1224539_ZJET : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2013_I1224539_ZJET()
      : Analysis("CMS_2013_I1224539_ZJET"),
        _filter(fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3))),
        _trimmer(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03))),
        _pruner(fastjet::Pruner(fastjet::cambridge_algorithm, 0.1, 0.5))
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(Cuts::abseta < 2.4);
      declare(fs, "FS");

      // Find Zs with pT > 120 GeV
      ZFinder zfinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 30*GeV, PID::ELECTRON, 80*GeV, 100*GeV,
                      0.2, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);

      declare(zfinder, "ZFinder");

      // Z+jet jet collections
      declare(FastJets(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.7), "JetsAK7_zj");
      declare(FastJets(zfinder.remainingFinalState(), FastJets::CAM, 0.8), "JetsCA8_zj");
      declare(FastJets(zfinder.remainingFinalState(), FastJets::CAM, 1.2), "JetsCA12_zj");

      // Histograms
      /// @note These are 2D histos rendered into slices
      const int zjetsOffset = 28;
      for (size_t i = 0; i < N_PT_BINS_vj; ++i ) {
        _h_ungroomedJetMass_AK7_zj[i] = bookHisto1D(zjetsOffset+i+1+0*N_PT_BINS_vj, 1, 1);
        _h_filteredJetMass_AK7_zj[i] = bookHisto1D(zjetsOffset+i+1+1*N_PT_BINS_vj,1,1);
        _h_trimmedJetMass_AK7_zj[i] = bookHisto1D(zjetsOffset+i+1+2*N_PT_BINS_vj,1,1);
        _h_prunedJetMass_AK7_zj[i] = bookHisto1D(zjetsOffset+i+1+3*N_PT_BINS_vj,1,1);
        _h_prunedJetMass_CA8_zj[i] = bookHisto1D(zjetsOffset+i+1+4*N_PT_BINS_vj,1,1);
        if (i > 0) _h_filteredJetMass_CA12_zj[i] = bookHisto1D(zjetsOffset+i+5*N_PT_BINS_vj,1,1);
      }
    }


    bool isBackToBack_zj(const ZFinder& zf, const fastjet::PseudoJet& psjet) {
      const FourMomentum& z = zf.bosons()[0].momentum();
      const FourMomentum& l1 = zf.constituents()[0].momentum();
      const FourMomentum& l2 = zf.constituents()[1].momentum();
      /// @todo We should make FourMomentum know how to construct itself from a PseudoJet
      const FourMomentum jmom(psjet.e(), psjet.px(), psjet.py(), psjet.pz());
      return (deltaPhi(z, jmom) > 2.0 && deltaR(l1, jmom) > 1.0 && deltaR(l2, jmom) > 1.0);
    }


    // Find the pT histogram bin index for value pt (in GeV), to hack a 2D histogram equivalent
    /// @todo Use a YODA axis/finder alg when available
    size_t findPtBin(double ptJ) {
      const double ptBins_vj[N_PT_BINS_vj+1] = { 125.0, 150.0, 220.0, 300.0, 450.0 };
      for (size_t ibin = 0; ibin < N_PT_BINS_vj; ++ibin) {
        if (inRange(ptJ, ptBins_vj[ibin], ptBins_vj[ibin+1])) return ibin;
      }
      return N_PT_BINS_vj;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get the Z
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;
      const Particle& z = zfinder.bosons()[0];
      const Particle& l1 = zfinder.constituents()[0];
      const Particle& l2 = zfinder.constituents()[1];

      // Require a high-pT Z (and constituents)
      if (l1.pT() < 30*GeV || l2.pT() < 30*GeV || z.pT() < 120*GeV) vetoEvent;

      // AK7 jets
      const PseudoJets& psjetsAK7_zj = apply<FastJets>(event, "JetsAK7_zj").pseudoJetsByPt(50.0*GeV);
      if (!psjetsAK7_zj.empty()) {
        // Get the leading jet and make sure it's back-to-back with the Z
        const fastjet::PseudoJet& j0 = psjetsAK7_zj[0];
        if (isBackToBack_zj(zfinder, j0)) {
          const size_t njetBin = findPtBin(j0.pt()/GeV);
          if (njetBin < N_PT_BINS_vj) {
            fastjet::PseudoJet filtered0 = _filter(j0);
            fastjet::PseudoJet trimmed0 = _trimmer(j0);
            fastjet::PseudoJet pruned0 = _pruner(j0);
            _h_ungroomedJetMass_AK7_zj[njetBin]->fill(j0.m()/GeV, weight);
            _h_filteredJetMass_AK7_zj[njetBin]->fill(filtered0.m()/GeV, weight);
            _h_trimmedJetMass_AK7_zj[njetBin]->fill(trimmed0.m()/GeV, weight);
            _h_prunedJetMass_AK7_zj[njetBin]->fill(pruned0.m()/GeV, weight);
          }
        }
      }

      // CA8 jets
      const PseudoJets& psjetsCA8_zj = apply<FastJets>(event, "JetsCA8_zj").pseudoJetsByPt(50.0*GeV);
      if (!psjetsCA8_zj.empty()) {
        // Get the leading jet and make sure it's back-to-back with the Z
        const fastjet::PseudoJet& j0 = psjetsCA8_zj[0];
        if (isBackToBack_zj(zfinder, j0)) {
          const size_t njetBin = findPtBin(j0.pt()/GeV);
          if (njetBin < N_PT_BINS_vj) {
            fastjet::PseudoJet pruned0 = _pruner(j0);
            _h_prunedJetMass_CA8_zj[njetBin]->fill(pruned0.m()/GeV, weight);
          }
        }
      }

      // CA12 jets
      const PseudoJets& psjetsCA12_zj = apply<FastJets>(event, "JetsCA12_zj").pseudoJetsByPt(50.0*GeV);
      if (!psjetsCA12_zj.empty()) {
        // Get the leading jet and make sure it's back-to-back with the Z
        const fastjet::PseudoJet& j0 = psjetsCA12_zj[0];
        if (isBackToBack_zj(zfinder, j0)) {
          const size_t njetBin = findPtBin(j0.pt()/GeV);
          if (njetBin>0 && njetBin < N_PT_BINS_vj) {
            fastjet::PseudoJet filtered0 = _filter(j0);
            _h_filteredJetMass_CA12_zj[njetBin]->fill( filtered0.m() / GeV, weight);
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double normalizationVal = 1000;
      for (size_t i = 0; i < N_PT_BINS_vj; ++i ) {
        normalize( _h_ungroomedJetMass_AK7_zj[i], normalizationVal);
        normalize( _h_filteredJetMass_AK7_zj[i], normalizationVal);
        normalize( _h_trimmedJetMass_AK7_zj[i], normalizationVal);
        normalize( _h_prunedJetMass_AK7_zj[i], normalizationVal);
        normalize( _h_prunedJetMass_CA8_zj[i], normalizationVal);
        if (i > 0) normalize( _h_filteredJetMass_CA12_zj[i], normalizationVal);
      }
    }

    //@}


  private:

    /// @name FastJet grooming tools (configured in constructor init list)
    //@{
    const fastjet::Filter _filter;
    const fastjet::Filter _trimmer;
    const fastjet::Pruner _pruner;
    //@}


    /// @name Histograms
    //@{
    enum BINS_vj { PT_125_150_vj=0, PT_150_220_vj, PT_220_300_vj, PT_300_450_vj, N_PT_BINS_vj };
    Histo1DPtr _h_ungroomedJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_trimmedJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_CA8_zj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_CA12_zj[N_PT_BINS_vj];
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1224539_ZJET);

}
#line 1 "CMS_2013_I1256943.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// CMS cross-section and angular correlations in Z boson + b-hadrons events at 7 TeV
  class CMS_2013_I1256943 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1256943()
      : Analysis("CMS_2013_I1256943")
    {     }


    /// Add projections and book histograms
    void init() {
      _sumW = 0;
      _sumW50 = 0;
      _sumWpT = 0;

      FinalState fs(Cuts::abseta < 2.4 && Cuts::pT > 20*GeV);
      declare(fs, "FS");

      UnstableFinalState ufs(Cuts::abseta < 2 && Cuts::pT > 15*GeV);
      declare(ufs, "UFS");

      Cut zetacut = Cuts::abseta < 2.4;

      ZFinder zfindermu(fs, zetacut, PID::MUON, 81.0*GeV, 101.0*GeV, 0.1, ZFinder::NOCLUSTER, ZFinder::TRACK, 91.2*GeV);
      declare(zfindermu, "ZFinderMu");

      ZFinder zfinderel(fs, zetacut, PID::ELECTRON, 81.0*GeV, 101.0*GeV, 0.1, ZFinder::NOCLUSTER, ZFinder::TRACK, 91.2*GeV);
      declare(zfinderel, "ZFinderEl");


      // Histograms in non-boosted region of Z pT
      _h_dR_BB = bookHisto1D(1, 1, 1);
      _h_dphi_BB = bookHisto1D(2, 1, 1);
      _h_min_dR_ZB = bookHisto1D(3, 1, 1);
      _h_A_ZBB = bookHisto1D(4, 1, 1);

      // Histograms in boosted region of Z pT (pT > 50 GeV)
      _h_dR_BB_boost = bookHisto1D(5, 1, 1);
      _h_dphi_BB_boost = bookHisto1D(6, 1, 1);
      _h_min_dR_ZB_boost = bookHisto1D(7, 1, 1);
      _h_A_ZBB_boost = bookHisto1D(8, 1, 1);

      _h_min_ZpT = bookHisto1D(9,1,1);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      vector<FourMomentum> Bmom;

      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");
      const ZFinder& zfindermu = apply<ZFinder>(e, "ZFinderMu");
      const ZFinder& zfinderel = apply<ZFinder>(e, "ZFinderEl");

      // Look for a Z --> mu+ mu- event in the final state
      if (zfindermu.empty() && zfinderel.empty()) vetoEvent;

      const Particles& z = !zfindermu.empty() ? zfindermu.bosons() : zfinderel.bosons();
      const bool is_boosted = ( z[0].pT() > 50*GeV );

      // Loop over the unstable particles
      foreach (const Particle& p, ufs.particles()) {
        const PdgId pid = p.pid();

        // Look for particles with a bottom quark
        if (PID::hasBottom(pid)) {

          bool good_B = false;
          const GenParticle* pgen = p.genParticle();
          const GenVertex* vgen = pgen -> end_vertex();

          // Loop over the decay products of each unstable particle.
          // Look for a couple of B hadrons.
          for (GenVertex::particles_out_const_iterator it = vgen->particles_out_const_begin(); it !=  vgen->particles_out_const_end(); ++it) {
            // If the particle produced has a bottom quark do not count it and go to the next loop cycle.
            if (!( PID::hasBottom( (*it)->pdg_id() ) ) ) {
              good_B = true;
              continue;
            } else {
              good_B = false;
              break;
            }
          }
          if (good_B ) Bmom.push_back( p.momentum() );
        }
        else continue;
      }

      // If there are more than two B's in the final state veto the event
      if (Bmom.size() != 2 ) vetoEvent;

      // Calculate the observables
      double dphiBB = deltaPhi(Bmom[0], Bmom[1]);
      double dRBB = deltaR(Bmom[0], Bmom[1]);

      const FourMomentum& pZ = z[0].momentum();
      const bool closest_B = ( deltaR(pZ, Bmom[0]) < deltaR(pZ, Bmom[1]) );
      const double mindR_ZB = closest_B ? deltaR(pZ, Bmom[0]) : deltaR(pZ, Bmom[1]);
      const double maxdR_ZB = closest_B ? deltaR(pZ, Bmom[1]) : deltaR(pZ, Bmom[0]);
      const double AZBB = ( maxdR_ZB - mindR_ZB ) / ( maxdR_ZB + mindR_ZB );

      // Get event weight for histogramming
      const double weight = e.weight();

      // Fill the histograms in the non-boosted region
      _h_dphi_BB->fill(dphiBB, weight);
      _h_dR_BB->fill(dRBB, weight);
      _h_min_dR_ZB->fill(mindR_ZB, weight);
      _h_A_ZBB->fill(AZBB, weight);
      _sumW += weight;
      _sumWpT += weight;

      // Fill the histograms in the boosted region
      if (is_boosted) {
        _sumW50 += weight;
        _h_dphi_BB_boost->fill(dphiBB, weight);
        _h_dR_BB_boost->fill(dRBB, weight);
        _h_min_dR_ZB_boost->fill(mindR_ZB, weight);
        _h_A_ZBB_boost->fill(AZBB, weight);
      }

      // Fill Z pT (cumulative) histogram
      _h_min_ZpT->fill(0, weight);
      if (pZ.pT() > 40*GeV ) {
        _sumWpT += weight;
        _h_min_ZpT->fill(40, weight);
      }
      if (pZ.pT() > 80*GeV ) {
        _sumWpT += weight;
        _h_min_ZpT->fill(80, weight);
      }
      if (pZ.pT() > 120*GeV ) {
        _sumWpT += weight;
        _h_min_ZpT->fill(120, weight);
      }

      Bmom.clear();
    }


    /// Finalize
    void finalize() {

      // Normalize excluding overflow bins (d'oh)
      normalize(_h_dR_BB, 0.7*crossSection()*_sumW/sumOfWeights(), false);  // d01-x01-y01
      normalize(_h_dphi_BB, 0.53*crossSection()*_sumW/sumOfWeights(), false);   // d02-x01-y01
      normalize(_h_min_dR_ZB, 0.84*crossSection()*_sumW/sumOfWeights(), false); // d03-x01-y01
      normalize(_h_A_ZBB, 0.2*crossSection()*_sumW/sumOfWeights(), false);  // d04-x01-y01

      normalize(_h_dR_BB_boost, 0.84*crossSection()*_sumW50/sumOfWeights(), false); // d05-x01-y01
      normalize(_h_dphi_BB_boost, 0.63*crossSection()*_sumW50/sumOfWeights(), false);   // d06-x01-y01
      normalize(_h_min_dR_ZB_boost, 1*crossSection()*_sumW50/sumOfWeights(), false);    // d07-x01-y01
      normalize(_h_A_ZBB_boost, 0.25*crossSection()*_sumW50/sumOfWeights(), false); // d08-x01-y01

      normalize(_h_min_ZpT, 40*crossSection()*_sumWpT/sumOfWeights(), false);   // d09-x01-y01
    }


  private:

    /// @name Weight counters
    //@{
    double _sumW, _sumW50, _sumWpT;
    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _h_dphi_BB, _h_dR_BB, _h_min_dR_ZB, _h_A_ZBB;
    Histo1DPtr _h_dphi_BB_boost, _h_dR_BB_boost, _h_min_dR_ZB_boost, _h_A_ZBB_boost, _h_min_ZpT;
    //@}

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1256943);

}
#line 1 "CMS_2013_I1258128.cc"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {


  /// CMS Z rapidity measurement
  class CMS_2013_I1258128 : public Analysis {
  public:

    // Constructor
    CMS_2013_I1258128()
      : Analysis("CMS_2013_I1258128")
    {    }


    void init() {
      // Full final state
      const FinalState fs(Cuts::abseta < 5);
      declare(fs, "FS");

      // Z finders for electrons and muons
      Cut cuts = Cuts::abseta < 2.1 && Cuts::pT > 20*GeV;
      const ZFinder zfe(fs, cuts, PID::ELECTRON, 76*GeV, 106*GeV);
      const ZFinder zfm(fs, cuts, PID::MUON, 76*GeV, 106*GeV);
      declare(zfe, "ZFE");
      declare(zfm, "ZFM");

      // Try to get the leading photon
      LeadingParticlesFinalState photonfs(FinalState(-2.5, 2.5, 40.0*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // Jets
      const FastJets jets(fs, FastJets::ANTIKT, 0.5);
      declare(jets, "JETS");

      // Histograms
      _hist1YZ      = bookHisto1D(1, 1, 1);
      _hist1YJet    = bookHisto1D(2, 1, 1);
      _hist1YSum    = bookHisto1D(3, 1, 1);
	  _hist1YDif	= bookHisto1D(4, 1, 1);
	  _hist2YPhoton = bookHisto1D(5, 1, 1);
	  _hist2YJet	= bookHisto1D(6, 1, 1);
	  _hist2YSum	= bookHisto1D(7, 1, 1);
	  _hist2YDif	= bookHisto1D(8, 1, 1);
    }


    void makeZCut(const Event& event) {
      // Apply the Z finders and veto if no Z found
      const ZFinder& zfe = apply<ZFinder>(event, "ZFE");
      const ZFinder& zfm = apply<ZFinder>(event, "ZFM");
      if (zfe.empty() && zfm.empty()) vetoEvent;

      // Choose the Z candidate
      const ParticleVector& z = (!zfm.empty()) ? zfm.bosons() : zfe.bosons();
      const ParticleVector& clusteredConstituents = (!zfm.empty()) ? zfm.constituents() : zfe.constituents();

      // Insist that the Z is in a high-pT (boosted) regime
      if (z[0].pT() < 40*GeV) return;

      // Build the jets
      const FastJets& jetfs = apply<FastJets>(event, "JETS");
      Jets jets = jetfs.jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 2.4);
      if (jets.empty()) return;

      // Clean the jets against the lepton candidates with a DeltaR cut of 0.5
      vector<const Jet*> cleanedJets;
      foreach (const Jet& j, jets) {
        bool isolated = true;
        foreach (const Particle& p, clusteredConstituents) {
          if (deltaR(p, j) < 0.5) {
            isolated = false;
            break;
          }
        }
        if (isolated) cleanedJets.push_back(&j);
      }
      // Require exactly 1 isolated jet
      if (cleanedJets.size() != 1) return;

      // Fill histos
      const double weight = event.weight();
      const double yz = z[0].rapidity();
      const double yjet = cleanedJets[0]->momentum().rapidity();
      _hist1YZ->fill(fabs(yz), weight);
      _hist1YJet->fill(fabs(yjet), weight);
      _hist1YSum->fill(0.5*fabs(yz + yjet), weight);
      _hist1YDif->fill(0.5*fabs(yz - yjet), weight);
    }


    void makePhotonCut(const Event& event) {
        // Get the photon
        const FinalState& photonfs = apply<FinalState>(event, "LeadingPhoton");
        if (photonfs.particles().size() < 1) return;
        const Particle& photon = photonfs.particles().front();
        if (photon.pT() < 40*GeV) return;
        if (fabs(photon.eta()) > 1.4442 ) return;

      // Build the jets
      const FastJets& jetfs = apply<FastJets>(event, "JETS");
      Jets jets = jetfs.jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 2.4);
      if (jets.empty()) return;

      // Clean the jets against the photon candidate with a DeltaR cut of 0.5
      vector<const Jet*> cleanedJets;
      foreach (const Jet& j, jets)
        if (deltaR(photon, j) > 0.5)
          cleanedJets.push_back(&j);
      // Require exactly 1 jet
      if (cleanedJets.size() != 1) return;

      // Fill histos
      const double weight = event.weight();
      const double ypho = photon.rapidity();
      const double yjet = cleanedJets[0]->momentum().rapidity();
      _hist2YPhoton->fill(fabs(ypho), weight);
      _hist2YJet->fill(fabs(yjet), weight);
      _hist2YSum->fill(0.5*fabs(ypho + yjet), weight);
      _hist2YDif->fill(0.5*fabs(ypho - yjet), weight);
    }


    void analyze(const Event& event) {
      makeZCut(event);
      makePhotonCut(event);
    }


    void finalize() {
      normalizeByContents(_hist1YZ);
      normalizeByContents(_hist1YJet);
      normalizeByContents(_hist1YSum);
      normalizeByContents(_hist1YDif);
      normalizeByContents(_hist2YPhoton);
      normalizeByContents(_hist2YJet);
      normalizeByContents(_hist2YSum);
      normalizeByContents(_hist2YDif);
    }


    // The CMS normalization in this analysis is that the sum over bin contents
    // is equal to 1. This function normalizes to area = area*bin_width.  /
    // @note This is a strange definition... why?
    void normalizeByContents(Histo1DPtr h) {
      normalize(h, h->bin(0).xWidth());
    }


  private:

    Histo1DPtr _hist1YZ, _hist1YJet, _hist1YSum, _hist1YDif;
    Histo1DPtr _hist2YPhoton, _hist2YJet, _hist2YSum, _hist2YDif;

  };


  // Plugin system hook
  DECLARE_RIVET_PLUGIN(CMS_2013_I1258128);

}
#line 1 "CMS_2013_I1261026.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Jet and underlying event properties as a function of particle multiplicity
  class CMS_2013_I1261026 : public Analysis {
  public:

    CMS_2013_I1261026()
      : Analysis("CMS_2013_I1261026"), _jetStructNorm(5,0.), _multBinCent(5,0.),
	_jetCounter5GeV(5,0.), _jetCounter30GeV(5,0.), _passedEv(5,0.)
    {  }


    void init() {
      FastJets jetpro(ChargedFinalState(-2.4, 2.4, 0.25*GeV), FastJets::ANTIKT, 0.5);
      declare(jetpro, "Jets");

      const ChargedFinalState cfs(-2.4, 2.4, 0.25*GeV);
      declare(cfs, "CFS250");

      // For min bias trigger
      const ChargedFinalState cfsBSCplus(3.23, 4.65, 500*MeV);
      declare(cfsBSCplus, "cfsBSCplus");

      const ChargedFinalState cfsBSCminus(-4.65, -3.23, 500*MeV);
      declare(cfsBSCminus, "cfsBSCminus");

      // Histograms:
      _h_AllTrkMeanPt            = bookProfile1D(1, 1, 1);
      _h_SoftTrkMeanPt           = bookProfile1D(2, 1, 1);
      _h_IntrajetTrkMeanPt       = bookProfile1D(3, 1, 1);
      _h_IntrajetLeaderTrkMeanPt = bookProfile1D(4, 1, 1);
      _h_MeanJetPt               = bookProfile1D(5, 1, 1);
      _h_JetRate5GeV             = bookProfile1D(6, 1, 1);
      _h_JetRate30GeV            = bookProfile1D(7, 1, 1);

      for (int ihist = 0; ihist < 5; ++ihist) {
        _h_JetSpectrum[ihist] = bookHisto1D(ihist+8, 1, 1);
        _h_JetStruct[ihist]   = bookHisto1D(ihist+13, 1, 1);

        // Temp histograms for distribution parameters and SEM calculation
        _th_AllTrkSpectrum[ihist]  = Histo1D(200, 0.0, 20.0);
        _th_SoftTrkSpectrum[ihist] = Histo1D(100, 0.0, 15.0);
        _th_JetTrkSpectrum[ihist]  = Histo1D(100, 0.0, 20.0);
        _th_JetLTrkSpectrum[ihist] = Histo1D(100, 0.0, 20.0);
      }

      _multBinCent[0] = 20;
      _multBinCent[1] = 40;
      _multBinCent[2] = 65;
      _multBinCent[3] = 95;
      _multBinCent[4] = 125;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // MinBias trigger
      const ChargedFinalState& cfsBSCplus = apply<ChargedFinalState>(event, "cfsBSCplus");
      if (cfsBSCplus.empty()) vetoEvent;
      const ChargedFinalState& cfsBSCminus = apply<ChargedFinalState>(event, "cfsBSCminus");
      if (cfsBSCminus.empty()) vetoEvent;

      const ChargedFinalState& cfsp = apply<ChargedFinalState>(event, "CFS250");
      if (cfsp.empty()) vetoEvent;

      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets& jets = jetpro.jetsByPt(5.0*GeV);

      const int mult = cfsp.size();

      int multbin[6] = { 10, 30, 50, 80, 110, 140 };
      for (int ibin = 0; ibin < 5; ++ibin) {
        if (mult > multbin[ibin] && mult <= multbin[ibin + 1]) {
          _passedEv[ibin] += weight;
          eventDecomp(event, ibin, weight);

          for (size_t ijets = 0; ijets < jets.size(); ijets++) {
            if (jets[ijets].abseta() < 1.9) {
              _h_JetSpectrum[ibin]->fill(jets[ijets].pT()/GeV, weight);
              if (jets[ijets].pT() > 5*GeV) _jetCounter5GeV[ibin] += weight;
              if (jets[ijets].pT() > 30*GeV) _jetCounter30GeV[ibin] += weight;
            }
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t i = 0; i < 5; ++i) {
        // All trk mean pT vs Nch
        _h_AllTrkMeanPt->fill(_multBinCent[i], _th_AllTrkSpectrum[i].xMean(), getMeanError(_th_AllTrkSpectrum[i]));

        // Soft trk mean pT vs Nch
        _h_SoftTrkMeanPt->fill(_multBinCent[i], _th_SoftTrkSpectrum[i].xMean(), getMeanError(_th_SoftTrkSpectrum[i]));

        // Intrajet trk mean pT vs Nch
        _h_IntrajetTrkMeanPt->fill(_multBinCent[i], _th_JetTrkSpectrum[i].xMean(), getMeanError(_th_JetTrkSpectrum[i]));

        // Intrajet leader trk mean pT vs Nch
        _h_IntrajetLeaderTrkMeanPt->fill(_multBinCent[i], _th_JetLTrkSpectrum[i].xMean(), getMeanError(_th_JetLTrkSpectrum[i]));

        // Jet mean pT vs Nch
        const double sem = (_h_JetSpectrum[i]->xStdDev())/(sqrt(_h_JetSpectrum[i]->sumW())) / _h_JetSpectrum[i]->xMean();
        _h_MeanJetPt->fill(_multBinCent[i], _h_JetSpectrum[i]->xMean(), sem);

        // Jet rates
	double avJetRate5  = _jetCounter5GeV[i]  / _passedEv[i];
        double avJetRate30 = _jetCounter30GeV[i] / _passedEv[i];

        const double sem5 = (_jetCounter5GeV[i] != 0) ? 1 / sqrt(_jetCounter5GeV[i]) : 0;
        _h_JetRate5GeV->fill(_multBinCent[i], avJetRate5, sem5);

        const double sem30 = (_jetCounter30GeV[i] != 0) ?  1 / sqrt(_jetCounter30GeV[i]) : 0;
        _h_JetRate30GeV->fill(_multBinCent[i], avJetRate30, sem30);

        scale(_h_JetSpectrum[i], 4.0 / _jetCounter5GeV[i]);
        scale(_h_JetStruct[i], 0.08 / _jetStructNorm[i]);
      }

    }


    double getMeanError(const Histo1D& hist) {
      double sem = hist.xStdErr(); // Standard error of the mean
      return sem / hist.xMean(); // relative SEM
    }


    void eventDecomp(const Event& event, size_t ibin, double weight) {

      struct TrkInJet { double pt; double eta; double phi; double R; };
      TrkInJet jetConstituents[100][100]; //1-st index - the number of the jet, 2-nd index - track in the jet
      TrkInJet jetsEv[100];
      size_t j[100];
      size_t jCount = 0;

      for (size_t i = 0; i < 100; i++) {
        j[i] = 0;
        jetsEv[i].pt = 0;
        jetsEv[i].eta = 0;
        jetsEv[i].phi = 0;
        for (size_t k = 0; k < 100; k++) {
          jetConstituents[i][k].pt = 0;
          jetConstituents[i][k].phi = 0;
          jetConstituents[i][k].eta = 0;
          jetConstituents[i][k].R = 0;
        }
      }

      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets& jets = jetpro.jetsByPt(5.0*GeV);

      // Start event decomp

      for (size_t ijets = 0; ijets < jets.size(); ijets++) {
        jetsEv[ijets].pt = jets[ijets].pT();
        jetsEv[ijets].eta = jets[ijets].eta();
        jetsEv[ijets].phi = jets[ijets].phi();
        jCount++;
      }

      const ChargedFinalState& cfsp = apply<ChargedFinalState>(event, "CFS250");
      foreach (const Particle& p, cfsp.particles()) {
        _th_AllTrkSpectrum[ibin].fill(p.pT()/GeV, weight);
        int flag = 0;
        for (size_t i = 0; i < jCount; i++) {
          const double delta_phi = deltaPhi(jetsEv[i].phi, p.phi());
          const double delta_eta = jetsEv[i].eta - p.eta();
          const double R = sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
          if (R <= 0.5) {
            flag++;
            jetConstituents[i][j[i]].pt = p.pT();
            jetConstituents[i][j[i]].R = R;
            j[i]++;
          }
        }
        if (flag == 0) _th_SoftTrkSpectrum[ibin].fill(p.pT()/GeV, weight);
      }

      for (size_t i = 0; i < jCount; i++) {
        double ptInjetLeader = 0;
        if (!inRange(jetsEv[i].eta, -1.9, 1.9)) continue; // only fully reconstructed jets for internal jet studies
        for (size_t k = 0; k < j[i]; k++) {
          _th_JetTrkSpectrum[ibin].fill(jetConstituents[i][k].pt , weight);
          _h_JetStruct[ibin]->fill(jetConstituents[i][k].R, jetConstituents[i][k].pt/jetsEv[i].pt);
          _jetStructNorm[ibin] += jetConstituents[i][k].pt / jetsEv[i].pt;
          if (ptInjetLeader < jetConstituents[i][k].pt) ptInjetLeader = jetConstituents[i][k].pt;
        }
        if (ptInjetLeader != 0) _th_JetLTrkSpectrum[ibin].fill(ptInjetLeader, weight);
      }

    }


  private:

    // Counters etc.
    vector<double> _jetStructNorm;
    vector<double> _multBinCent;
    /// @todo Need to handle weights
    vector<double> _jetCounter5GeV, _jetCounter30GeV, _passedEv;

    Profile1DPtr _h_AllTrkMeanPt, _h_SoftTrkMeanPt;
    Profile1DPtr _h_IntrajetTrkMeanPt, _h_IntrajetLeaderTrkMeanPt;
    Profile1DPtr _h_MeanJetPt;
    Profile1DPtr _h_JetRate5GeV, _h_JetRate30GeV;

    Histo1DPtr _h_JetSpectrum[5];
    Histo1DPtr _h_JetStruct[5];

    // Temp histograms
    Histo1D _th_AllTrkSpectrum[5];
    Histo1D _th_SoftTrkSpectrum[5];
    Histo1D _th_JetTrkSpectrum[5];
    Histo1D _th_JetLTrkSpectrum[5];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1261026);

}
#line 1 "CMS_2013_I1265659.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CMS_2013_I1265659 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1265659()
      : Analysis("CMS_2013_I1265659")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {
      const FastJets jets(FinalState(-10, 10, 0.0*GeV), FastJets::ANTIKT, 0.5);
      declare(jets, "Jets");

      _h_hTotD = bookHisto1D(1, 1, 1);
      _h_hTotDF = bookHisto1D(1, 1, 2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(30.0*GeV);
      if (jets.size() < 3) vetoEvent;

      const FourMomentum jet1 = jets[0].momentum();
      const FourMomentum jet2 = jets[1].momentum();
      const FourMomentum jet3 = jets[2].momentum();

      // Cut on lead jet pT and lead/sublead jet centrality
      if (jet1.pT() < 100*GeV) vetoEvent;
      if (jet1.abseta() > 2.5 || jet2.abseta() > 2.5) vetoEvent;

      // Construct eta & phi distances between 2nd and 3rd jets
      double dEta23 = jet3.eta() - jet2.eta(); ///< Note not abs
      double dPhi23 = jet3.phi() - jet2.phi(); ///< Note not abs
      if (dPhi23 > M_PI)  dPhi23 -= 2*M_PI; ///< @todo Use mapTo... functions?
      if (dPhi23 < -M_PI) dPhi23 += 2*M_PI; ///< @todo Use mapTo... functions?

      // Cut on distance between 2nd and 3rd jets
      const double R23 = add_quad(dPhi23, dEta23);
      if (!inRange(R23, 0.5, 1.5)) vetoEvent;

      // Cut on dijet mass
      const FourMomentum diJet = jet1 + jet2;
      if (diJet.mass() < 220*GeV) vetoEvent;

      // Calc beta and fill histogram (choose central or fwd histo inline)
      double beta = fabs(atan2(dPhi23, sign(jet2.eta())*dEta23));
      ((jet2.abseta() < 0.8) ? _h_hTotD : _h_hTotDF)->fill(beta, event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double width = _h_hTotD->bin(0).xWidth();
      normalize(_h_hTotD, width);
      normalize(_h_hTotDF, width);
    }


  private:

    /// @name Histograms
    Histo1DPtr _h_hTotD;
    Histo1DPtr _h_hTotDF;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1265659);

}
#line 1 "CMS_2013_I1272853.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"

namespace Rivet {


  /// CMS W + 2 jet double parton scattering analysis
  class CMS_2013_I1272853 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1272853()
      : Analysis("CMS_2013_I1272853") { }


    /// Book histograms and initialise projections before the run
    void init() {

       const FinalState fs;
       declare(fs, "FS");

       /// @todo Use C++11 initialisation syntax
       vector<PdgIdPair> vidsW;
       vidsW += make_pair(PID::MUON, PID::NU_MUBAR), make_pair(PID::ANTIMUON, PID::NU_MU);
       InvMassFinalState invfsW(fs, vidsW, 20*GeV, 1e6*GeV);
       declare(invfsW, "INVFSW");

       VetoedFinalState vfs(fs);
       vfs.addVetoOnThisFinalState(invfsW);
       declare(vfs, "VFS");
       declare(FastJets(vfs, FastJets::ANTIKT, 0.5), "Jets");

       _h_deltaS_eq2jet_Norm = bookHisto1D(1,1,1);
       _h_rel_deltaPt_eq2jet_Norm = bookHisto1D(2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Find Ws
      const InvMassFinalState& invMassFinalStateW = apply<InvMassFinalState>(event, "INVFSW");
      if (invMassFinalStateW.empty()) vetoEvent;
      const Particles& WDecayProducts = invMassFinalStateW.particles();
      if (WDecayProducts.size() < 2) vetoEvent;

      // Cuts on W decay properties
      const int iNU_MU = (WDecayProducts[1].abspid() == PID::NU_MU) ? 1 : 0;
      const int iAN_MU = 1 - iNU_MU;
      const double pt1  = WDecayProducts[iAN_MU].pT();
      const double pt2  = WDecayProducts[iNU_MU].Et();
      const double eta1 = WDecayProducts[iAN_MU].abseta();
      const double phi1 = WDecayProducts[iAN_MU].phi();
      const double phi2 = WDecayProducts[iNU_MU].phi();
      const double mt   = sqrt(2 * pt1 * pt2 * (1 - cos(phi1-phi2)));
      if (mt < 50*GeV || pt1 < 35*GeV || eta1 > 2.1 || pt2 < 30*GeV) vetoEvent;

      // Get jets and make sure there are at least two of them in |y| < 2
      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      /// @todo Collapse this into jetpro.jetsByPt(ptGtr(20*GeV) & rapIn(2.0))
      vector<FourMomentum> jets;
      foreach (const Jet& jet, jetpro.jetsByPt(20*GeV))
        if (jet.absrap() < 2.0) jets.push_back(jet.momentum());
      if (jets.size() != 2) vetoEvent;

      const double mupx     = pt1 * cos(phi1);
      const double mupy     = pt1 * sin(phi1);
      const double met_x    = pt2 * cos(phi2);
      const double met_y    = pt2 * sin(phi2);
      const double dpt      = add_quad(jets[0].px() + jets[1].px(), jets[0].py() + jets[1].py());
      const double rel_dpt  = dpt / (jets[0].pT() + jets[1].pT());
      const double pT2      = sqr(mupx + met_x) + sqr(mupy + met_y);
      const double Px       = (mupx + met_x)*(jets[0].px() + jets[1].px());
      const double Py       = (mupy + met_y)*(jets[0].py() + jets[1].py());
      const double p1p2_mag = dpt * sqrt(pT2);
      const double dS       = acos((Px+Py) / p1p2_mag);

      const double weight = event.weight();
      _h_rel_deltaPt_eq2jet_Norm->fill(rel_dpt, weight);
      _h_deltaS_eq2jet_Norm->fill(dS, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double rel_dpt_bw = 1.0002 / 30.0;
      const double dphi_bw = 3.14160 / 30.0;
      normalize(_h_rel_deltaPt_eq2jet_Norm, rel_dpt_bw);
      normalize(_h_deltaS_eq2jet_Norm, dphi_bw);
    }


  private:

    Histo1DPtr _h_rel_deltaPt_eq2jet_Norm;
    Histo1DPtr _h_deltaS_eq2jet_Norm;

  };



  DECLARE_RIVET_PLUGIN(CMS_2013_I1272853);

}
#line 1 "CMS_2013_I1273574.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CMS 4-jet production at 7 TeV
  class CMS_2013_I1273574 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1273574()
      : Analysis("CMS_2013_I1273574")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState cnfs(-4.7, 4.7);
      declare(FastJets(cnfs, FastJets::ANTIKT, 0.5), "Jets");

      _h_jetetas[0]     = bookHisto1D(1,1,1);
      _h_jetpts[0]      = bookHisto1D(2,1,1);
      _h_DeltaS         = bookHisto1D(3,1,1);
      _h_DeltaPhiSoft   = bookHisto1D(4,1,1);
      _h_DeltaPtRelSoft = bookHisto1D(5,1,1);
      _h_jetetas[2]     = bookHisto1D(6,1,1);
      _h_jetpts[2]      = bookHisto1D(7,1,1);
      _h_jetetas[3]     = bookHisto1D(8,1,1);
      _h_jetpts[3]      = bookHisto1D(9,1,1);
      _h_jetetas[1]     = bookHisto1D(10,1,1);
      _h_jetpts[1]      = bookHisto1D(11,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      /// @todo Use jetsByPt(ptGtr(20*GeV) & absetaIn(4.7)), then no need for the lower loop;
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(20*GeV);
      if (jets.size() < 4) vetoEvent;

      // Ensure that there are exactly 4 jets > 20 GeV, with two above 50 GeV
      Jets hardjets, alljets;
      foreach (const Jet& j, jets) {
        if (j.abseta() > 4.7) continue;
        if (j.pT() > 50*GeV) hardjets.push_back(j);
        if (j.pT() > 20*GeV) alljets.push_back(j);
      }
      if (hardjets.size() < 2 || alljets.size() != 4) vetoEvent;
      const double weight = event.weight();

      // Histogram pT and eta of all 4 jets
      for (size_t i = 0; i < 4; ++i) {
        _h_jetpts[i]->fill(alljets[i].pT()/GeV, weight);
        _h_jetetas[i]->fill(alljets[i].eta(), weight);
      }

      // Create vector sums of the hard and soft pairs of jets
      const FourMomentum p12 = alljets[0].momentum() + alljets[1].momentum();
      const FourMomentum p34 = alljets[2].momentum() + alljets[3].momentum();

      // Fill the delta(phi) between the soft jets
      const double dphisoft = deltaPhi(alljets[2], alljets[3]);
      _h_DeltaPhiSoft->fill(dphisoft, weight);

      // Fill the pT balance between the soft jets
      const double ptbalanceSoft = p34.pT() / (alljets[2].pT() + alljets[3].pT());
      _h_DeltaPtRelSoft->fill(ptbalanceSoft, weight);

      // Fill the azimuthal angle difference between the two jet pairs
      const double p12p34_trans = p12.px()*p34.px() + p12.py()*p34.py();
      const double DeltaS = acos( p12p34_trans / p12.pT() / p34.pT() );
      _h_DeltaS->fill(DeltaS, weight);
    }


    /// Normalise histograms (mostly to cross-section)
    void finalize() {
      const double invlumi = crossSection()/picobarn/sumOfWeights();
      for (size_t i = 0; i < 4; ++i) {
        scale(_h_jetpts[i], invlumi);
        scale(_h_jetetas[i], invlumi);
      }
      normalize(_h_DeltaPtRelSoft);
      normalize(_h_DeltaPhiSoft);
      normalize(_h_DeltaS);
    }


  private:

    Histo1DPtr _h_jetpts[4], _h_jetetas[4];
    Histo1DPtr _h_DeltaS, _h_DeltaPhiSoft, _h_DeltaPtRelSoft;

  };


  DECLARE_RIVET_PLUGIN(CMS_2013_I1273574);

}
#line 1 "CMS_2014_I1298810.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Ratios of jet pT spectra, related to ratios of differential jet cross sections
  class CMS_2014_I1298810 : public Analysis {
  public:

    /// Constructor
    CMS_2014_I1298810()
      : Analysis("CMS_2014_I1298810")
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      // Projections
      FastJets jetsak5(FinalState(), FastJets::ANTIKT, 0.5);
      declare(jetsak5, "JetsAK5");
      FastJets jetsak7(FinalState(), FastJets::ANTIKT, 0.7);
      declare(jetsak7, "JetsAK7");

      // Histograms
      _h_pt_05_ak5    = bookHisto1D(1, 1, 1);
      _h_pt_05_10_ak5 = bookHisto1D(2, 1, 1);
      _h_pt_10_15_ak5 = bookHisto1D(3, 1, 1);
      _h_pt_15_20_ak5 = bookHisto1D(4, 1, 1);
      _h_pt_20_25_ak5 = bookHisto1D(5, 1, 1);
      _h_pt_25_30_ak5 = bookHisto1D(6, 1, 1);

      _h_pt_05_ak7    = bookHisto1D(7, 1, 1);
      _h_pt_05_10_ak7 = bookHisto1D(8, 1, 1);
      _h_pt_10_15_ak7 = bookHisto1D(9, 1, 1);
      _h_pt_15_20_ak7 = bookHisto1D(10, 1, 1);
      _h_pt_20_25_ak7 = bookHisto1D(11, 1, 1);
      _h_pt_25_30_ak7 = bookHisto1D(12, 1, 1);

      _h_pt_05_ratio    = bookScatter2D(13, 1, 1);
      _h_pt_05_10_ratio = bookScatter2D(14, 1, 1);
      _h_pt_10_15_ratio = bookScatter2D(15, 1, 1);
      _h_pt_15_20_ratio = bookScatter2D(16, 1, 1);
      _h_pt_20_25_ratio = bookScatter2D(17, 1, 1);
      _h_pt_25_30_ratio = bookScatter2D(18, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Jets& jetsak5 = apply<FastJets>(event, "JetsAK5").jetsByPt(56*GeV);
      const Jets& jetsak7 = apply<FastJets>(event, "JetsAK7").jetsByPt(56*GeV);
      if (jetsak5.size() < 1 && jetsak7.size() < 1) vetoEvent;

      const double weight = event.weight();

      // Filling R = 0.5 jets
      foreach(const Jet& jet, jetsak5) {
        if (jet.absrapidity() < 0.5) {
          _h_pt_05_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 1.0) {
          _h_pt_05_10_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 1.5) {
          _h_pt_10_15_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 2.0) {
          _h_pt_15_20_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 2.5) {
          _h_pt_20_25_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 3.0) {
          _h_pt_25_30_ak5->fill(jet.pT()/GeV, weight);
        }
      }


      // Filling R = 0.7 jets
      foreach(const Jet& jet, jetsak7) {
        if (jet.absrapidity() < 0.5) {
          _h_pt_05_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 1.0) {
          _h_pt_05_10_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 1.5) {
          _h_pt_10_15_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 2.0) {
          _h_pt_15_20_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 2.5) {
          _h_pt_20_25_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 3.0) {
          _h_pt_25_30_ak7->fill(jet.pT() * GeV, weight);
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pt_05_ak5,    crossSection()/sumOfWeights());
      scale(_h_pt_05_10_ak5, crossSection()/sumOfWeights());
      scale(_h_pt_10_15_ak5, crossSection()/sumOfWeights());
      scale(_h_pt_15_20_ak5, crossSection()/sumOfWeights());
      scale(_h_pt_20_25_ak5, crossSection()/sumOfWeights());
      scale(_h_pt_25_30_ak5, crossSection()/sumOfWeights());

      scale(_h_pt_05_ak7,    crossSection()/sumOfWeights());
      scale(_h_pt_05_10_ak7, crossSection()/sumOfWeights());
      scale(_h_pt_10_15_ak7, crossSection()/sumOfWeights());
      scale(_h_pt_15_20_ak7, crossSection()/sumOfWeights());
      scale(_h_pt_20_25_ak7, crossSection()/sumOfWeights());
      scale(_h_pt_25_30_ak7, crossSection()/sumOfWeights());

      divide(_h_pt_05_ak5,    _h_pt_05_ak7,    _h_pt_05_ratio);
      divide(_h_pt_05_10_ak5, _h_pt_05_10_ak7, _h_pt_05_10_ratio);
      divide(_h_pt_10_15_ak5, _h_pt_10_15_ak7, _h_pt_10_15_ratio);
      divide(_h_pt_15_20_ak5, _h_pt_15_20_ak7, _h_pt_15_20_ratio);
      divide(_h_pt_20_25_ak5, _h_pt_20_25_ak7, _h_pt_20_25_ratio);
      divide(_h_pt_25_30_ak5, _h_pt_25_30_ak7, _h_pt_25_30_ratio);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_pt_05_ak5, _h_pt_05_10_ak5, _h_pt_10_15_ak5, _h_pt_15_20_ak5, _h_pt_20_25_ak5, _h_pt_25_30_ak5;
    Histo1DPtr _h_pt_05_ak7, _h_pt_05_10_ak7, _h_pt_10_15_ak7, _h_pt_15_20_ak7, _h_pt_20_25_ak7, _h_pt_25_30_ak7;
    Scatter2DPtr _h_pt_05_ratio, _h_pt_05_10_ratio, _h_pt_10_15_ratio, _h_pt_15_20_ratio, _h_pt_20_25_ratio, _h_pt_25_30_ratio;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2014_I1298810);

}
#line 1 "CMS_2014_I1303894.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// @brief Differential cross-section of W bosons + jets in pp collisions at sqrt(s)=7 TeV
  /// @author Darin Baumgartel (darinb@cern.ch)
  ///
  /// Based on Rivet analysis originally created by Anil Singh (anil@cern.ch), Lovedeep Saini (lovedeep@cern.ch)
  class CMS_2014_I1303894 : public Analysis {
  public:

    /// Constructor
    CMS_2014_I1303894()
      : Analysis("CMS_2014_I1303894")
    {   }


    // Book histograms and initialise projections before the run
    void init() {
      // Projections
      const FinalState fs;
      declare(fs, "FS");

      MissingMomentum missing(fs);
      declare(missing, "MET");

      IdentifiedFinalState bareMuons(fs);
      bareMuons.acceptIdPair(PID::MUON);
      DressedLeptons muonClusters(fs, bareMuons, 0.1, Cuts::open(), false, false);
      declare(muonClusters, "muonClusters");

      IdentifiedFinalState neutrinos;
      neutrinos.acceptIdPair(PID::NU_MU);
      declare(neutrinos, "neutrinos");

      VetoedFinalState jetFS(fs);
      jetFS.addVetoOnThisFinalState(muonClusters);
      jetFS.addVetoOnThisFinalState(neutrinos);
      jetFS.vetoNeutrinos();
      FastJets JetProjection(jetFS, FastJets::ANTIKT, 0.5);
      JetProjection.useInvisibles(false);
      declare(JetProjection, "Jets");

      // Histograms
      _histDPhiMuJet1 = bookHisto1D(1,1,1);
      _histDPhiMuJet2 = bookHisto1D(2,1,1);
      _histDPhiMuJet3 = bookHisto1D(3,1,1);
      _histDPhiMuJet4 = bookHisto1D(4,1,1);

      _histEtaJet1 = bookHisto1D(5,1,1);
      _histEtaJet2 = bookHisto1D(6,1,1);
      _histEtaJet3 = bookHisto1D(7,1,1);
      _histEtaJet4 = bookHisto1D(8,1,1);

      _histHT1JetInc = bookHisto1D(9,1,1);
      _histHT2JetInc = bookHisto1D(10,1,1);
      _histHT3JetInc = bookHisto1D(11,1,1);
      _histHT4JetInc = bookHisto1D(12,1,1);

      _histJet30MultExc  = bookHisto1D(13,1,1);
      _histJet30MultInc  = bookHisto1D(14,1,1);

      _histPtJet1 = bookHisto1D(15,1,1);
      _histPtJet2 = bookHisto1D(16,1,1);
      _histPtJet3 = bookHisto1D(17,1,1);
      _histPtJet4 = bookHisto1D(18,1,1);

      // Counters
      _n_1jet = 0.0;
      _n_2jet = 0.0;
      _n_3jet = 0.0;
      _n_4jet = 0.0;
      _n_inclusivebinsummation = 0.0;
    }


    void analyze(const Event& event) {
      // Get the dressed muon
      const DressedLeptons& muonClusters = apply<DressedLeptons>(event, "muonClusters");
      int nmu = muonClusters.dressedLeptons().size();
      if (nmu < 1) vetoEvent;
      DressedLepton dressedmuon = muonClusters.dressedLeptons()[0];
      if (dressedmuon.momentum().abseta() > 2.1) vetoEvent;
      if (dressedmuon.momentum().pT() < 25.0*GeV) vetoEvent;

      // Get the muon neutrino
      const Particles& neutrinos = apply<FinalState>(event, "neutrinos").particlesByPt();
      if (neutrinos.empty()) vetoEvent;

      // Check that the muon and neutrino are not decay products of tau
      if (dressedmuon.constituentLepton().hasAncestor( PID::TAU)) vetoEvent;
      if (dressedmuon.constituentLepton().hasAncestor(-PID::TAU)) vetoEvent;
      if (dressedmuon.constituentLepton().hasAncestor( PID::NU_TAU)) vetoEvent;
      if (dressedmuon.constituentLepton().hasAncestor(-PID::NU_TAU)) vetoEvent;
      if (neutrinos[0].hasAncestor( PID::TAU)) vetoEvent;
      if (neutrinos[0].hasAncestor(-PID::TAU)) vetoEvent;
      if (neutrinos[0].hasAncestor( PID::NU_TAU)) vetoEvent;
      if (neutrinos[0].hasAncestor(-PID::NU_TAU)) vetoEvent;

      // Recording of event weight and numbers
      const double weight = event.weight();

      // Get the missing momentum
      const MissingMomentum& met = apply<MissingMomentum>(event, "MET");
      const double ptmet = met.visibleMomentum().pT();
      const double phimet = (-met.visibleMomentum()).phi();

      // Calculate MET and MT(mu,MET), and remove events with MT < 50 GeV
      const double ptmuon = dressedmuon.pT();
      const double phimuon = dressedmuon.phi();
      const double mt_mumet = sqrt(2*ptmuon*ptmet*(1.0 - cos(phimet-phimuon)));

      // Remove events in MT < 50 region
      if (mt_mumet < 50*GeV) vetoEvent;

      // Loop over jets and fill pt/eta/phi quantities in vectors
      const Jets& jets_filtered = apply<FastJets>(event, "Jets").jetsByPt(0.0*GeV);
      vector<float> finaljet_pT_list, finaljet_eta_list, finaljet_phi_list;
      double htjets = 0.0;
      for (size_t ii = 0; ii < jets_filtered.size(); ++ii) {
        // Jet pT/eta/phi
        double jet_pt = jets_filtered[ii].pT();
        double jet_eta = jets_filtered[ii].eta();
        double jet_phi = jets_filtered[ii].phi();

        // Kinemetic cuts for jet acceptance
        if (fabs(jet_eta) > 2.4) continue;
        if (jet_pt < 30.0*GeV) continue;
        if (deltaR(dressedmuon, jets_filtered[ii]) < 0.5) continue;

        // Add jet to jet list and increases the HT variable
        finaljet_pT_list.push_back(jet_pt);
        finaljet_eta_list.push_back(jet_eta);
        finaljet_phi_list.push_back(jet_phi);
        htjets += fabs(jet_pt);
      }


      // Filling of histograms:
      // Fill as many jets as there are into the exclusive jet multiplicity
      if (!finaljet_pT_list.empty())
        _histJet30MultExc->fill(finaljet_pT_list.size(), weight);

      for (size_t ij = 0; ij < finaljet_pT_list.size(); ++ij) {
        _histJet30MultInc->fill(ij+1, weight);
        _n_inclusivebinsummation += weight;
      }

      if (finaljet_pT_list.size() >= 1) {
        _histPtJet1->fill(finaljet_pT_list[0],weight);
        _histEtaJet1->fill(fabs(finaljet_eta_list[0]), weight);
        _histDPhiMuJet1->fill(deltaPhi(finaljet_phi_list[0], phimuon),weight);
        _histHT1JetInc->fill(htjets, weight);
        _n_1jet +=weight;
      }

      if (finaljet_pT_list.size() >= 2) {
        _histPtJet2->fill(finaljet_pT_list[1], weight);
        _histEtaJet2->fill(fabs(finaljet_eta_list[1]), weight);
        _histDPhiMuJet2->fill(deltaPhi(finaljet_phi_list[1], phimuon), weight);
        _histHT2JetInc->fill(htjets, weight);
        _n_2jet += weight;
      }

      if (finaljet_pT_list.size() >= 3) {
        _histPtJet3->fill(finaljet_pT_list[2], weight);
        _histEtaJet3->fill(fabs(finaljet_eta_list[2]), weight);
        _histDPhiMuJet3->fill(deltaPhi(finaljet_phi_list[2], phimuon), weight);
        _histHT3JetInc->fill(htjets, weight);
        _n_3jet += weight;
      }

      if (finaljet_pT_list.size() >=4 ) {
        _histPtJet4->fill(finaljet_pT_list[3], weight);
        _histEtaJet4->fill(fabs(finaljet_eta_list[3]), weight);
        _histDPhiMuJet4->fill(deltaPhi(finaljet_phi_list[3], phimuon), weight);
        _histHT4JetInc-> fill(htjets, weight);
        _n_4jet += weight;
      }

    }


    // Finalize the histograms.
    void finalize() {

      const double inclusive_cross_section = crossSection();
      const double norm_1jet_histo = inclusive_cross_section*_n_1jet/sumOfWeights();
      const double norm_2jet_histo = inclusive_cross_section*_n_2jet/sumOfWeights();
      const double norm_3jet_histo = inclusive_cross_section*_n_3jet/sumOfWeights();
      const double norm_4jet_histo = inclusive_cross_section*_n_4jet/sumOfWeights();
      const double norm_incmultiplicity = inclusive_cross_section*_n_inclusivebinsummation/sumOfWeights();

      normalize(_histJet30MultExc, norm_1jet_histo);
      normalize(_histJet30MultInc, norm_incmultiplicity);

      normalize(_histPtJet1, norm_1jet_histo);
      normalize(_histHT1JetInc, norm_1jet_histo);
      normalize(_histEtaJet1, norm_1jet_histo);
      normalize(_histDPhiMuJet1, norm_1jet_histo);

      normalize(_histPtJet2, norm_2jet_histo);
      normalize(_histHT2JetInc, norm_2jet_histo);
      normalize(_histEtaJet2, norm_2jet_histo);
      normalize(_histDPhiMuJet2, norm_2jet_histo);

      normalize(_histPtJet3, norm_3jet_histo);
      normalize(_histHT3JetInc, norm_3jet_histo);
      normalize(_histEtaJet3, norm_3jet_histo);
      normalize(_histDPhiMuJet3, norm_3jet_histo);

      normalize(_histPtJet4, norm_4jet_histo);
      normalize(_histHT4JetInc, norm_4jet_histo);
      normalize(_histEtaJet4, norm_4jet_histo);
      normalize(_histDPhiMuJet4, norm_4jet_histo);
    }


  private:

    Histo1DPtr  _histJet30MultExc, _histJet30MultInc;
    Histo1DPtr  _histPtJet1, _histPtJet2, _histPtJet3, _histPtJet4;
    Histo1DPtr  _histEtaJet1, _histEtaJet2, _histEtaJet3, _histEtaJet4;
    Histo1DPtr  _histDPhiMuJet1, _histDPhiMuJet2, _histDPhiMuJet3, _histDPhiMuJet4;
    Histo1DPtr  _histHT1JetInc, _histHT2JetInc, _histHT3JetInc, _histHT4JetInc;

    double _n_1jet, _n_2jet, _n_3jet, _n_4jet, _n_inclusivebinsummation;

  };


  DECLARE_RIVET_PLUGIN(CMS_2014_I1303894);

}
#line 1 "CMS_2014_I1305624.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  namespace {

    /// Number of event shape variables
    /// @todo Move into the EventShape class
    const int NEVTVAR = 5;

    /// Number of leading jet pT thresholds
    /// @todo Move into the analysis class
    const int NJETPTMN = 5;
    /// Leading jet pT thresholds
    /// @todo Move into the analysis class
    const double LEADINGPTTHRESHOLD[NJETPTMN] = { 110.0, 170.0, 250.0, 320.0, 390.0 };


    // Helpers for event shape calculations in hidden namespace; implementation at bottom of file
    /// @todo Why a class? Improve/remove this junk
    class EventShape {
    public:

      /// Constructor from vectors of four-vectors as input objects in the event to calculate the event shapes
      EventShape(const vector<double>& px_vector, const vector<double>& py_vector, const vector<double>& pz_vector,
                 const vector<double>& e_vector, double eta_central, int irap, int nmn)
        : _object_px(px_vector), _object_py(py_vector), _object_pz(pz_vector),
          _object_e(e_vector), _eta_c(eta_central), _irap(irap), _nmnjet(nmn)
      {   }

      /// @brief Returns the values of the five event shapes
      ///
      /// Event shape indices:
      /// 0. central transverse thrust
      /// 1. central total jet broadening
      /// 2. central total jet mass
      /// 3. central total transverse jet mass
      /// 4. central three-jet resolution threshold
      vector<double> getEventShapes() {
        _calculate(); ///< @todo There should be some test for success/failure!!
        return _event_shapes;
      }

      /// Returns the global thrust axis Nx, Ny, Nz=0
      vector<double> getThrustAxis() {
        _calculate(); ///< @todo There should be some test for success/failure!!
        return _thrust_axis;
      }

      /// Returns the central thrust axis Nx, Ny, Nz=0
      vector<double> getThrustAxisC() {
        _calculate(); ///< @todo There should be some test for success/failure!!
        return _thrust_axis_c;
      }

      // /// @brief Choice of the central region
      // void setEtaC(double eta_central) { _eta_c = eta_central; }

      // // Whether to use the rapidity y (rap==1)  or the pseudorapidity eta (rap==0)
      // void setRapType(int irap) { _irap = irap; }


    private:

      /// Calculate everything
      int _calculate();

      /// Returns the difference in phi between two vectors
      double _delta_phi(double, double);

      /// The Lorentz scalar product
      double _lorentz_sp(const vector<double>&, const vector<double>&);

      // Calculates the three-jet resolutions
      double _three_jet_res(const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, int);

      // Calculates the thrust axis and the tau values
      vector<double> _thrust(const vector<double>&, const vector<double>&);


      vector<double> _object_px, _object_py, _object_pz, _object_p;
      vector<double> _object_pt, _object_e, _object_phi, _object_eta;
      vector<double> _event_shapes;
      vector<double> _thrust_axis, _thrust_axis_c;

      double _eta_c;
      int _irap;
      size_t _nmnjet;

    };

  }




  class CMS_2014_I1305624 : public Analysis {
  public:

    /// Constructor
    CMS_2014_I1305624()
      : Analysis("CMS_2014_I1305624")
    {    }


    /// @name Analysis methods

    /// Book histograms and initialise projections before the run
    void init() {
      const FastJets jets(FinalState(Cuts::abseta < 2.6), FastJets::ANTIKT, 0.5);
      declare(jets, "Jets");

      for (int ij=0; ij < NJETPTMN; ij++) {
        _h_thrustc[ij] = bookHisto1D(1, 1, ij+1);
        _h_broadt[ij] = bookHisto1D(1, 2, ij+1);
        _h_tot3dmass[ij] = bookHisto1D(1, 3, ij+1);
        _h_tottrnsmass[ij] = bookHisto1D(1, 4, ij+1);
        _h_y23c[ij] = bookHisto1D(1, 5, ij+1);
        //
        _alow1[ij] = _h_thrustc[ij]->xMin();
        _alow2[ij] = _h_broadt[ij]->xMin();
        _alow3[ij] = _h_tot3dmass[ij]->xMin();
        _alow4[ij] = _h_tottrnsmass[ij]->xMin();
        _alow5[ij] = _h_y23c[ij]->xMin();
        //
        _ahgh1[ij] = _h_thrustc[ij]->xMax();
        _ahgh2[ij] = _h_broadt[ij]->xMax();
        _ahgh3[ij] = _h_tot3dmass[ij]->xMax();
        _ahgh4[ij] = _h_tottrnsmass[ij]->xMax();
        _ahgh5[ij] = _h_y23c[ij]->xMax();
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(30.0*GeV);
      if (jets.size() < 2) vetoEvent;
      if (jets[0].abseta() > 2.4 || jets[1].abseta() > 2.4) vetoEvent;

      const double leadingpt = jets[0].pT();
      if (leadingpt < 110*GeV) vetoEvent;

      vector<double> jtpx, jtpy, jtpz, jten;
      foreach (const Jet& j, jets) {
        if (j.abseta() < 2.4) {
          jtpx.push_back(j.px());
          jtpy.push_back(j.py());
          jtpz.push_back(j.pz());
          jten.push_back(j.E());
        }
      }

      EventShape eventshape(jtpx, jtpy, jtpz, jten, 2.4, 0, 2);
      const vector<double> eventvar = eventshape.getEventShapes();
      if (eventvar[NEVTVAR] < 0) vetoEvent; // Jets are not only one hemisphere

      const double weight = event.weight();
      for (int ij = NJETPTMN-1; ij >= 0; --ij) {
        if (leadingpt/GeV > LEADINGPTTHRESHOLD[ij]) {
          if (inRange(eventvar[0], _alow1[ij], _ahgh1[ij])) _h_thrustc[ij]->fill(eventvar[0], weight);
          if (inRange(eventvar[2], _alow3[ij], _ahgh3[ij])) _h_tot3dmass[ij]->fill(eventvar[2], weight);
          if (inRange(eventvar[3], _alow4[ij], _ahgh4[ij])) _h_tottrnsmass[ij]->fill(eventvar[3], weight);
          if (eventvar[NEVTVAR] >= 3) {
            if (inRange(eventvar[1], _alow2[ij], _ahgh2[ij])) _h_broadt[ij]->fill(eventvar[1], weight);
            if (inRange(eventvar[4], _alow5[ij], _ahgh5[ij])) _h_y23c[ij]->fill(eventvar[4], weight);
          }
          break;
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (int ij = 0; ij < NJETPTMN; ij++) {
        normalize(_h_thrustc[ij]);
        normalize(_h_broadt[ij]);
        normalize(_h_tot3dmass[ij]);
        normalize(_h_tottrnsmass[ij]);
        normalize(_h_y23c[ij]);
      }
    }


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_thrustc[NJETPTMN];
    Histo1DPtr _h_broadt[NJETPTMN];
    Histo1DPtr _h_tot3dmass[NJETPTMN];
    Histo1DPtr _h_tottrnsmass[NJETPTMN];
    Histo1DPtr _h_y23c[NJETPTMN];
    //@}

    // Data members
    double _alow1[NJETPTMN], _alow2[NJETPTMN], _alow3[NJETPTMN], _alow4[NJETPTMN], _alow5[NJETPTMN];
    double _ahgh1[NJETPTMN], _ahgh2[NJETPTMN], _ahgh3[NJETPTMN], _ahgh4[NJETPTMN], _ahgh5[NJETPTMN];

  };


  DECLARE_RIVET_PLUGIN(CMS_2014_I1305624);



  /////////////////////


  namespace {

    // EventShape helper class method implementations:

    int EventShape::_calculate() {
      if (!_event_shapes.empty() && !_thrust_axis.empty() && !_thrust_axis_c.empty())
        return 1; //< return success if this appears to already have been run

      const size_t length = (size_t) _object_px.size();

      if (((size_t) _object_py.size() != length) ||
          ((size_t) _object_pz.size() != length) ||
          ((size_t) _object_e.size() != length)) {
        /// @todo Change to exception or assert
        // cout << "ERROR!!!! Input vectors differ in size! Change that please!" << endl;
        // cout<<"py_size: "<<_object_py.size()<<" ,pz_size: "<<_object_pz.size()
        //     <<" ,px_size: "<<_object_px.size()<<" ,E_size: "<<_object_e.size()<<endl;
        return 0;
      }

      if (!_object_p.empty()) {
        _object_p.clear();
        _object_pt.clear();
        _object_eta.clear();
        _object_phi.clear();
        _event_shapes.clear();
        _thrust_axis.clear();
        _thrust_axis_c.clear();
      }

      for (size_t j = 0; j < length; j++) {
        _object_p.push_back(0.);
        _object_pt.push_back(0.);
        _object_eta.push_back(0.);
        _object_phi.push_back(0.);
      }

      for (int j = 0; j < NEVTVAR; j++) {
        _event_shapes.push_back(-50.);
      }

      _event_shapes.push_back(double(_object_px.size())); //< WTF?

      for (size_t j = 0; j < 3; j++) {
        _thrust_axis.push_back(0.);
        _thrust_axis_c.push_back(0.);
      }

      double theta = 0;

      for (size_t k = 0; k < length; k++) {
        _object_p[k] = sqrt(pow(_object_px[k],2) + pow(_object_py[k],2) + pow(_object_pz[k],2));
        _object_pt[k] = sqrt(pow(_object_px[k],2) + pow(_object_py[k],2));
        if (_object_p[k] > _object_e[k] + 1e-4) {
          /// @todo Change to exception or assert
          // cout << "ERROR!!! object " << k <<" has P = " << _object_p[k]
          //      << " which is bigger than E = " << _object_e[k] <<" "
          //      << _object_px[k] <<" "<< _object_py[k] <<" "
          //      << _object_pz[k] <<" of total length "<< length
          //      << endl;
          return 0;
        }

        //to prevent a division by zero
        if (_irap == 0) {
          if (fabs(_object_pz[k]) > 1e-5) {
            theta = atan(_object_pt[k]/(_object_pz[k]));
          } else {
            theta = M_PI/2;
          }
          if (theta < 0.) theta = theta + M_PI;
          _object_eta[k] = -log(tan(0.5*theta));
        }
        if (_irap == 1) {
          if (_object_pz[k] == _object_e[k]) {
            /// @todo Change to exception or assert
            // cout << "ERROR!!! object "<<k<<" has Pz "<< _object_pz[k] <<" which is equal to E = "<< _object_e[k] <<endl;
            return 0;
          }
          _object_eta[k]=0.5*log((_object_e[k]+_object_pz[k])/(_object_e[k]-_object_pz[k]));
        }
        if (_irap != 0 && _irap != 1) {
          /// @todo Change to exception or assert
          // cout << "ERROR!!!, The choice to use the rapidity y or the pseudorapidity eta is not set correctly! Change that please!" << endl;
          return 0;
        }
        _object_phi[k] = atan2(_object_py[k], _object_px[k]);
      }

      vector<double> object_px_in, object_py_in, object_pz_in, object_pt_in, object_e_in, object_et_in, object_eta_in;
      vector<double> object_px_out, object_py_out, object_pz_out, object_e_out, object_pt_out, object_eta_out;
      if (!object_px_in.empty()) { //< FFS, this is impossible: it's only just been created!
        object_px_in.clear();
        object_py_in.clear();
        object_pz_in.clear();
        object_pt_in.clear();
        object_e_in.clear();
        object_et_in.clear();
        object_eta_in.clear();
        object_px_out.clear();
        object_py_out.clear();
        object_pz_out.clear();
        object_pt_out.clear();
        object_e_out.clear();
        object_eta_out.clear();
      }

      size_t nin = 0;

      for (size_t j = 0; j < length; j++) {
        if (fabs(_object_eta[j]) < _eta_c) {
          object_px_in.push_back(_object_px[j]);
          object_py_in.push_back(_object_py[j]);
          object_pz_in.push_back(_object_pz[j]);
          object_e_in.push_back(_object_e[j]);
          object_pt_in.push_back(sqrt(pow(_object_px[j],2)+pow(_object_py[j],2)));
          object_et_in.push_back(sqrt((pow(_object_e[j],2)*pow(_object_pt[j],2))/(pow(_object_pt[j],2)+pow(_object_pz[j],2))));
          object_eta_in.push_back(_object_eta[j]);
          nin += 1;
      } else {
          object_px_out.push_back(_object_px[j]);
          object_py_out.push_back(_object_py[j]);
          object_pz_out.push_back(_object_pz[j]);
          object_e_out.push_back(_object_e[j]);
          object_pt_out.push_back(sqrt(pow(_object_px[j],2)+pow(_object_py[j],2)));
          object_eta_out.push_back(_object_eta[j]);
        }
      }

      if (object_px_in.size() != nin) {
        /// @todo Change to exception or assert
        cout<<"ERROR!!! wrong dimension of 'in' momenta"<<endl;
        //return 0; ///< @todo Why not do this?
      }
      const size_t nout = length - nin;

      if (nin < _nmnjet) {
        for (int i = 0; i < NEVTVAR; i++) {
          _event_shapes[i] = -50.0;
        }
      }

      _event_shapes[NEVTVAR] = nin;

      if (nin >= _nmnjet) {
        double p_sum_c = 0; //GMA
        double pt_sum_c = 0;
        double eta_cw=0;
        double px_sum_in = 0;
        double py_sum_in = 0;
        for (size_t j = 0; j < nin; j++) {
          pt_sum_c += object_pt_in[j];
          p_sum_c += sqrt(pow(object_pt_in[j],2.) + pow(object_pz_in[j], 2.0)); //GMA
          eta_cw += object_pt_in[j]*object_eta_in[j];
          px_sum_in += object_px_in[j];
          py_sum_in += object_py_in[j];
        }
        eta_cw /= pt_sum_c;

        double expTerm = 0;
        for (size_t j = 0; j < nout; j++) {
          expTerm += object_pt_out[j] * exp(-fabs(object_eta_out[j]-eta_cw));
        }
        expTerm /= pt_sum_c;

        //the central global transverse thrust centrthr is calculated
        double centrthr = 0;
        vector<double> thrust_central = _thrust(object_px_in, object_py_in);

        for (size_t l=0; l<3; l++) _thrust_axis_c[l] = thrust_central[l];
        //the variable which gets resummed is not thrust
        //but tau = 1 - thrust - see calculation
        centrthr = thrust_central[3];
        _event_shapes[0] = centrthr;

        double alpha_c = atan2(_thrust_axis_c[1], _thrust_axis_c[0]);
        //central jet masses
        //define two jet masses in region U and D
        double cenjm_up = 0;
        double cenjm_down= 0;
        double dot_product = 0;

        vector<double> up_sum;
        vector<double> down_sum;
        for (size_t j=0; j<4;j++) {
          up_sum.push_back(0.);
          down_sum.push_back(0.);
        }
        for (size_t i=0;i<nin;i++) {
          dot_product = object_px_in[i] * _thrust_axis_c[0] + object_py_in[i] * _thrust_axis_c[1];
          if (dot_product >= 0) {
            up_sum[0]+=object_px_in[i];
            up_sum[1]+=object_py_in[i];
            up_sum[2]+=object_pz_in[i];
            up_sum[3]+=object_e_in[i];
          } else {
            down_sum[0]+=object_px_in[i];
            down_sum[1]+=object_py_in[i];
            down_sum[2]+=object_pz_in[i];
            down_sum[3]+=object_e_in[i];
          }
        }
        cenjm_up = _lorentz_sp(up_sum, up_sum) / pow(p_sum_c, 2.); //GMA pow(pt_sum_c,2);
        cenjm_down = _lorentz_sp(down_sum, down_sum) / pow(p_sum_c, 2.); //GMA pow(pt_sum_c,2);

        //central total jet mass centotjm
        double centotjm=0;
        centotjm = cenjm_up + cenjm_down;

        _event_shapes[2]=centotjm;

        double centrjm_up=0, centrjm_down=0;
        vector<double> upsum;
        vector<double> downsum;
        for (size_t j = 0; j < 3; j++) {
          upsum.push_back(0.);
          downsum.push_back(0.);
        }
        for (size_t i = 0; i < nin; i++) {
          dot_product = object_px_in[i]*_thrust_axis_c[0]+object_py_in[i]*_thrust_axis_c[1];
          if (dot_product >= 0) {
            upsum[0] += object_px_in[i];
            upsum[1] += object_py_in[i];
            upsum[2] += object_et_in[i];
          } else {
            downsum[0] += object_px_in[i];
            downsum[1] += object_py_in[i];
            downsum[2] += object_et_in[i];
          }
        }
        centrjm_up = _lorentz_sp(upsum, upsum) / pow(pt_sum_c, 2);
        centrjm_down = _lorentz_sp(downsum, downsum) / pow(pt_sum_c, 2);
        double centottrjm = centrjm_up + centrjm_down;

        _event_shapes[3] = centottrjm;

        //central three-jet resolution threshold
        double ceny3=0;
        if (nin < 3) {
          ceny3 = -1.0;
        } else {
          ceny3 = _three_jet_res(object_px_in, object_py_in, object_pz_in, object_e_in, _irap);
        }

        _event_shapes[4] = ceny3;

        //the central jet broadenings in the up and down region
        double cenbroad_up=0;
        double cenbroad_down=0;

        double eta_up=0;
        size_t num_up=0;
        double eta_down =0;
        size_t num_down =0;
        double phi_temp =0;
        double phi_up_aver =0;
        double phi_down_aver =0;
        double pt_sum_up =0;
        double pt_sum_down =0;
        double dot_product_b =0;
        vector<double> phi_up;
        vector<double> phi_down;
        double py_rot =0;
        double px_rot =0;

        for (size_t j = 0; j < 4; j++) {
          up_sum.push_back(0.);
          down_sum.push_back(0.);
        }

        for (size_t i=0;i<nin;i++) {
          dot_product_b =sqrt(object_px_in[i]*_thrust_axis_c[0] + object_py_in[i]*_thrust_axis_c[1]);
          if (dot_product_b>=0){
            pt_sum_up += object_pt_in[i];
            //rotate the coordinate system so that
            //the central thrust axis is e_x
            px_rot = cos(alpha_c)*object_px_in[i]+sin(alpha_c)*object_py_in[i];
            py_rot = - sin(alpha_c)*object_px_in[i]+cos(alpha_c)*object_py_in[i];
            //calculate the eta and phi in the rotated system
            eta_up += object_pt_in[i]*object_eta_in[i];
            phi_temp = atan2(py_rot,px_rot);

            if(phi_temp > M_PI/2){
              phi_temp = phi_temp - M_PI/2;
            }
            if (phi_temp < -M_PI/2){
              phi_temp = phi_temp + M_PI/2;
            }
            phi_up.push_back(phi_temp);
            phi_up_aver += object_pt_in[i]*phi_temp;
            num_up += 1;
          } else {
            eta_down += object_pt_in[i]*object_eta_in[i];
            pt_sum_down += object_pt_in[i];
            px_rot = cos(alpha_c)*object_px_in[i]+sin(alpha_c)*object_py_in[i];
            py_rot = - sin(alpha_c)*object_px_in[i]+cos(alpha_c)*object_py_in[i];
            phi_temp = atan2(py_rot,px_rot);
            if (phi_temp > M_PI/2) {
              //if phi is bigger than pi/2 in the new system calculate
              //the difference to the thrust axis
              phi_temp = M_PI -phi_temp;
            }
            if (phi_temp<-M_PI/2) {
              //if phi is smaller than
              phi_temp = -M_PI-phi_temp;
            }
            phi_down.push_back(phi_temp);
            //calculate the pt-weighted phi
            phi_down_aver += object_pt_in[i]*phi_temp;
            num_down += 1;
          }
        }
        if (num_up!=0){
          eta_up = eta_up/pt_sum_up;
          phi_up_aver = phi_up_aver/pt_sum_up;
        }
        if (num_down!=0) {
          eta_down = eta_down/pt_sum_down;
          phi_down_aver = phi_down_aver/pt_sum_down;
        }

        size_t index_up=0, index_down=0;
        for (size_t i = 0; i < nin; i++) {
          dot_product_b = object_px_in[i]*_thrust_axis_c[0] + object_py_in[i]*_thrust_axis_c[1];
          if (dot_product_b >= 0) {
            //calculate the broadenings of the regions with the rotated system
            //and the pt-weighted average of phi in the rotated system
            cenbroad_up += object_pt_in[i]*sqrt(pow(object_eta_in[i]-eta_up, 2) +
                                                pow(_delta_phi(phi_up[index_up], phi_up_aver), 2));
            index_up += 1;
          } else {
            cenbroad_down += object_pt_in[i]*sqrt(pow(object_eta_in[i]-eta_down, 2)+
                                                  pow(_delta_phi(phi_down[index_down], phi_down_aver), 2));
            index_down += 1;
          }
        }

        if (index_up == 0 || index_down ==0) _event_shapes[NEVTVAR] *= -1;

        cenbroad_up=cenbroad_up/(2*pt_sum_c);
        cenbroad_down=cenbroad_down/(2*pt_sum_c);

        //central total jet broadening
        double centotbroad = 0;
        centotbroad = cenbroad_up + cenbroad_down;

        _event_shapes[1] = centotbroad;

        for (int ij = 0; ij < 5; ij++) {
          if (_event_shapes[ij] < 1.e-20) _event_shapes[ij] = 1.e-20;
          _event_shapes[ij] = log(_event_shapes[ij]);
        }
      }

      return 1;
    }


    double EventShape::_three_jet_res(const vector<double>& in_object_px, const vector<double>& in_object_py, const vector<double>& in_object_pz, const vector<double>& in_object_e, int irap) {

      size_t y3_length = (size_t)in_object_px.size();
      if (((size_t) in_object_py.size()!=y3_length) ||
          ((size_t) in_object_pz.size()!=y3_length) ||
          (in_object_e.size()!=y3_length)) {
        // cout << "ERROR!!!! Input vectors differ in size! Change that please!" << endl;
        // cout<<"py_size: "<<in_object_py.size()<<" ,pz_size: "<<in_object_pz.size()
        //     <<" ,px_size: "<<in_object_px.size()<<" , E_size: "<<in_object_e.size() <<endl;
        return 0.0;
      }

      vector<double> in_object_p, in_object_pt, in_object_eta, in_object_phi;
      if (!in_object_p.empty()) {
        in_object_p.clear();
        in_object_pt.clear();
        in_object_eta.clear();
        in_object_phi.clear();
      }
      for (size_t j = 0; j < y3_length; j++) {
        in_object_p.push_back(0.);
        in_object_pt.push_back(0.);
        in_object_eta.push_back(0.);
        in_object_phi.push_back(0.);
      }
      double theta_y3_1st = 0;
      for (size_t k =0; k<y3_length; k++) {
        in_object_p[k] = sqrt(pow(in_object_px[k],2) + pow(in_object_py[k],2) + pow(in_object_pz[k],2));
        in_object_pt[k] = sqrt(pow(in_object_px[k],2) + pow(in_object_py[k],2));

        //calculates the pseudorapidity to prevent a division by zero
        if (irap == 0) {
          if (fabs(in_object_pz[k]) > 1E-5) {
            theta_y3_1st = atan(in_object_pt[k]/(in_object_pz[k]));
          } else {
            theta_y3_1st = M_PI/2;
          }
          if (theta_y3_1st<0.) theta_y3_1st = theta_y3_1st + M_PI;
          in_object_eta[k] = - log(tan(0.5*theta_y3_1st));
        }
        //calculates the real rapidity
        if (irap == 1) {
          in_object_eta[k]=0.5*log((in_object_e[k]+in_object_pz[k])/(in_object_e[k]-in_object_pz[k]));
        }
        in_object_phi[k] = atan2(in_object_py[k], in_object_px[k]);
      }

      //the three-jet resolution
      //threshold y3
      double y3 = 0;

      //vector which will be filled with the
      //minimum of the distances
      double max_dmin_temp=0;

      double max_dmin = 0;

      //distance input object k, beam
      double distance_jB = 0;
      double distance_jB_min = 0;
      //distance of input object k to l
      double distance_jk = 0;
      double distance_jk_min = 0;
      //as we search the minimum of the distances
      //give them values which are for sure higher
      //than those we evaluate first in the for-loups

      size_t index_jB = 0;
      size_t index_j_jk = 0;
      size_t index_k_jk = 0;

      //to decide later if the minmum is a jB or jk
      int decide_jB = -1;

      vector<double> input_pt, input_px, input_py, input_pz;
      vector<double> input_p, input_e, input_phi, input_eta;

      if (!input_pt.empty()) {
        input_pt.clear();
        input_px.clear();
        input_px.clear();
        input_pz.clear();
        input_p.clear();
        input_e.clear();
        input_phi.clear();
        input_eta.clear();
      }

      for (size_t j = 0; j < y3_length; j++){
        input_pt.push_back(in_object_pt[j]);
        input_px.push_back(in_object_px[j]);
        input_py.push_back(in_object_py[j]);
        input_pz.push_back(in_object_pz[j]);
        input_p.push_back(in_object_p[j]);
        input_e.push_back(in_object_e[j]);
        input_phi.push_back(in_object_phi[j]);
        input_eta.push_back(in_object_eta[j]);
      }
      if (y3_length<3) {
        return -1;
      } else {
        size_t rest = y3_length;
        for (size_t i = 0; i<y3_length; i++) {
          //make the minima at the initialization step
          //of each looping bigger than the first values
          distance_jB_min = 0.36*pow(input_pt[0],2) + 10;
          //DELTA PHIs wanted not the pure difference
          distance_jk_min = min(pow(input_pt[1], 2), pow(input_pt[0], 2)) *
            (pow(input_eta[1]-input_eta[0], 2) +
             pow(_delta_phi(input_phi[1], input_phi[0]), 2)) + 10;
          //do the procedure only until we have only 2 objects left anymore
          if (rest > 2) {
            for (size_t j=0; j<rest;j++) {
              //calculate the distance between object j and the beam
              distance_jB = 0.36*pow(input_pt[j], 2);
              if(distance_jB < distance_jB_min){
                distance_jB_min = distance_jB;
                index_jB = j;
              }
              if (j > 0) {
                for(size_t k=0; k<j;k++){
                  //calculate the distance in delta eta and delta phi between object i and object j
                  distance_jk = min(pow(input_pt[j], 2),pow(input_pt[k], 2))*
                    (pow(input_eta[j]-input_eta[k], 2)+
                     pow(_delta_phi(input_phi[j],input_phi[k]), 2));
                  if (distance_jk<distance_jk_min) {
                    distance_jk_min = distance_jk;
                    index_j_jk = j;
                    index_k_jk =k;
                  }
                }
              }
            }
            //decide if the minimum is from a jB or jk combination
            if (distance_jk_min<distance_jB_min) {
              max_dmin_temp = max(distance_jk_min,max_dmin_temp);
              decide_jB = 0;
            } else {
              max_dmin_temp = max(distance_jB_min,max_dmin_temp);
              decide_jB=1;
            }
            //if we have only three jets left calculate
            //the maxima of the dmin's
            //if the minimum is a jB eliminate the input object
            if (decide_jB == 1) {
              //if index_jB is the last one nothing is to do
              if (index_jB != rest-1) {
                for (size_t i=index_jB; i<rest-1;i++) {
                  input_pt[i]=input_pt[i+1];
                  input_phi[i]=input_phi[i+1];
                  input_eta[i]=input_eta[i+1];
                  input_px[i]=input_px[i+1];
                  input_py[i]=input_py[i+1];
                  input_pz[i]=input_pz[i+1];
                  input_e[i]=input_e[i+1];
                }
              }
            }
            //if the minimum is a jk combine both input objects
            if(decide_jB==0) {
              input_px[index_k_jk] = input_px[index_k_jk]+input_px[index_j_jk];
              input_py[index_k_jk] = input_py[index_k_jk]+input_py[index_j_jk];
              input_pz[index_k_jk] = input_pz[index_k_jk]+input_pz[index_j_jk];
              input_e[index_k_jk] = input_e[index_k_jk]+input_e[index_j_jk];
              input_p[index_k_jk] = sqrt(pow(input_px[index_k_jk], 2)+
                                         pow(input_py[index_k_jk], 2)+
                                         pow(input_pz[index_k_jk], 2));
              //calculate the pt, eta and phi of the new combined momenta k_jk
              input_pt[index_k_jk] = sqrt(pow(input_px[index_k_jk], 2)+
                                          pow(input_py[index_k_jk], 2));
              //in the case of pseudorapidity
              if (irap == 0) {
                double theta_new =0;
                if (fabs(input_pz[index_k_jk]) > 1E-5){
                  theta_new = atan(input_pt[index_k_jk]/(input_pz[index_k_jk]));
                } else {
                  theta_new = M_PI/2;
                }
                if (theta_new < 0) {
                  theta_new = theta_new + M_PI;
                }
                input_eta[index_k_jk] = - log(tan(0.5*theta_new));
              }
              //in the real rapidity y is wanted
              if (irap == 1) {
                input_eta[index_k_jk] = 0.5 * log((input_e[index_k_jk]+
                                                   input_pz[index_k_jk]) /
                                                  (input_e[index_k_jk] -
                                                   input_pz[index_k_jk]));
              }
              input_phi[index_k_jk] = atan2(input_py[index_k_jk], input_px[index_k_jk]);
              if (index_j_jk != rest-1) {
                for (size_t i = index_j_jk; i<rest-1;i++) {
                  input_pt[i] = input_pt[i+1];
                  input_phi[i] = input_phi[i+1];
                  input_eta[i] = input_eta[i+1];
                  input_px[i] = input_px[i+1];
                  input_py[i] = input_py[i+1];
                  input_pz[i] = input_pz[i+1];
                  input_e[i] = input_e[i+1];
                }
              }
            }
          }
          if (rest == 3) max_dmin = max_dmin_temp;
          rest = rest-1;
        }
      }

      double et2 = 0;
      et2 = input_pt[0] + input_pt[1];
      y3 = max_dmin/pow(et2,2);

      return y3;
    }


    vector<double> EventShape::_thrust(const vector<double>& input_px, const vector<double>& input_py) {

      double thrustmax_calc = 0;
      double temp_calc = 0;
      size_t length_thrust_calc = 0;
      vector<double> thrust_values, thrust_axis_calc;
      vector<double> p_thrust_max_calc, p_dec_1_calc,  p_dec_2_calc, p_pt_beam_calc;

      if (!thrust_values.empty()){
        thrust_values.clear();
        thrust_axis_calc.clear();
        p_thrust_max_calc.clear();
        p_dec_1_calc.clear();
        p_dec_2_calc.clear();
        p_pt_beam_calc.clear();
      }

      for (size_t j = 0; j < 3; j++){
        p_pt_beam_calc.push_back(0.);
        p_dec_1_calc.push_back(0.);
        p_dec_2_calc.push_back(0.);
        p_thrust_max_calc.push_back(0.);
        thrust_axis_calc.push_back(0.);
      }

      for (size_t j = 0; j < 4; j++) {
        thrust_values.push_back(0.);
      }

      length_thrust_calc = input_px.size();
      if (input_py.size() != length_thrust_calc) {
        /// @todo Change to exception or assert
        cout<<"ERROR in thrust calculation!!! Size of input vectors differs. Change that please!"<<endl;
        return thrust_values;
      }

      double pt_sum_calc =0;
      for(size_t k=0;k<length_thrust_calc;k++){
        pt_sum_calc+=sqrt(pow(input_px[k],2)+pow(input_py[k],2));
        for(size_t j = 0; j < 3; j++){
          p_thrust_max_calc[j]=0;
        }
        //get a vector perpendicular to the beam axis and
        //perpendicular to the momentum of particle k
        //per default beam axis b = (0,0,1)
        p_pt_beam_calc[0] = input_py[k]*1;
        p_pt_beam_calc[1] = - input_px[k]*1;
        p_pt_beam_calc[2] = 0.; // GMA p_pt_beam_calc[3] = 0.;
        for(size_t i=0;i<length_thrust_calc;i++){
          if(i!=k){
            if((input_px[i]*p_pt_beam_calc[0]+input_py[i]*p_pt_beam_calc[1])>=0){
              p_thrust_max_calc[0]= p_thrust_max_calc[0] + input_px[i];
              p_thrust_max_calc[1]= p_thrust_max_calc[1] + input_py[i];
            }else{
              p_thrust_max_calc[0]= p_thrust_max_calc[0] - input_px[i];
              p_thrust_max_calc[1]= p_thrust_max_calc[1] - input_py[i];
            }
          }
        }
        p_dec_1_calc[0] = p_thrust_max_calc[0] + input_px[k];
        p_dec_1_calc[1] = p_thrust_max_calc[1] + input_py[k];
        p_dec_1_calc[2] = 0;
        p_dec_2_calc[0] = p_thrust_max_calc[0] - input_px[k];
        p_dec_2_calc[1] = p_thrust_max_calc[1] - input_py[k];
        p_dec_2_calc[2] = 0;
        temp_calc = pow(p_dec_1_calc[0], 2) + pow(p_dec_1_calc[1], 2);

        if (temp_calc>thrustmax_calc) {
          thrustmax_calc =temp_calc;
          for (size_t i=0; i<3; i++) {
            thrust_axis_calc[i] = p_dec_1_calc[i]/sqrt(thrustmax_calc);
          }
        }
        temp_calc = pow(p_dec_2_calc[0], 2)+pow(p_dec_2_calc[1], 2);
        if (temp_calc > thrustmax_calc) {
          thrustmax_calc =temp_calc;
          for (size_t i=0; i<3; i++) {
            thrust_axis_calc[i] = p_dec_2_calc[i]/sqrt(thrustmax_calc);
          }
        }
      }
      for (size_t j = 0; j < 3; j++) thrust_values[j] = thrust_axis_calc[j];
      const double thrust_calc = sqrt(thrustmax_calc)/pt_sum_calc;

      // the variable which gets returned is not the thrust but tau=1-thrust
      thrust_values[3] = 1 - thrust_calc;
      if (thrust_values[3] < 1e-20) thrust_values[3] = 1e-20;

      return thrust_values;
    }


    double EventShape::_delta_phi(double phi1, double phi2) {
      double dphi = fabs(phi2 - phi1);
      if (dphi > M_PI) dphi = 2*M_PI - dphi;
      return dphi;
    }


    // Returns the scalar product between two 4 momenta
    double EventShape::_lorentz_sp(const vector<double>& a, const vector<double>& b) {
      size_t dim = (size_t) a.size();
      if (a.size()!=b.size()) {
        cout<<"ERROR!!! Dimension of input vectors are different! Change that please!"<<endl;
        return 0;
      } else {
        double l_dot_product=a[dim-1]*b[dim-1];
        for(size_t i=0; i<dim-1;i++){
          l_dot_product-=a[i]*b[i];
        }
        return l_dot_product;
      }
    }


  }


}
#line 1 "CMS_2015_I1310737.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  class CMS_2015_I1310737 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2015_I1310737);


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs; ///< @todo No cuts?
      VisibleFinalState visfs(fs);

      ZFinder zeeFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::ELECTRON, 71.0*GeV, 111.0*GeV);
      declare(zeeFinder, "ZeeFinder");

      ZFinder zmumuFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::MUON, 71.0*GeV, 111.0*GeV);
      declare(zmumuFinder, "ZmumuFinder");

      VetoedFinalState jetConstits(visfs);
      jetConstits.addVetoOnThisFinalState(zeeFinder);
      jetConstits.addVetoOnThisFinalState(zmumuFinder);

      FastJets akt05Jets(jetConstits, FastJets::ANTIKT, 0.5);
      declare(akt05Jets, "AntiKt05Jets");


      _h_excmult_jets_tot = bookHisto1D(1, 1, 1);
      _h_incmult_jets_tot = bookHisto1D(2, 1, 1);
      _h_leading_jet_pt_tot = bookHisto1D(3, 1, 1);
      _h_second_jet_pt_tot = bookHisto1D(4, 1, 1);
      _h_third_jet_pt_tot = bookHisto1D(5, 1, 1);
      _h_fourth_jet_pt_tot = bookHisto1D(6, 1, 1);
      _h_leading_jet_eta_tot = bookHisto1D(7, 1, 1);
      _h_second_jet_eta_tot = bookHisto1D(8, 1, 1);
      _h_third_jet_eta_tot = bookHisto1D(9, 1, 1);
      _h_fourth_jet_eta_tot = bookHisto1D(10, 1, 1);
      _h_ht1_tot = bookHisto1D(11, 1, 1);
      _h_ht2_tot = bookHisto1D(12, 1, 1);
      _h_ht3_tot = bookHisto1D(13, 1, 1);
      _h_ht4_tot = bookHisto1D(14, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {;

      const ZFinder& zeeFS = apply<ZFinder>(event, "ZeeFinder");
      const ZFinder& zmumuFS = apply<ZFinder>(event, "ZmumuFinder");

      const Particles& zees = zeeFS.bosons();
      const Particles& zmumus = zmumuFS.bosons();

      // We did not find exactly one Z. No good.
      if (zees.size() + zmumus.size() != 1) {
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }

      // Find the (dressed!) leptons
      const Particles& dressedLeptons = zees.size() ? zeeFS.constituents() : zmumuFS.constituents();

      // Cluster jets
      // NB. Veto has already been applied on leptons and photons used for dressing
      const FastJets& fj = apply<FastJets>(event, "AntiKt05Jets");
      const Jets& jets = fj.jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 30*GeV);

      // Perform lepton-jet overlap and HT calculation
      double ht = 0;
      Jets goodjets;
      foreach (const Jet& j, jets) {
        // Decide if this jet is "good", i.e. isolated from the leptons
        /// @todo Nice use-case for any() and a C++11 lambda
        bool overlap = false;
        foreach (const Particle& l, dressedLeptons) {
          if (Rivet::deltaR(j, l) < 0.5) {
            overlap = true;
            break;
          }
        }

        // Fill HT and good-jets collection
        if (overlap) continue;
        goodjets.push_back(j);
        ht += j.pT();
      }

      // We don't care about events with no isolated jets
      if (goodjets.empty()) {
        MSG_DEBUG("No jets in event");
        vetoEvent;
      }


      /////////////////


      // Weight to be used for histo filling
      const double w = 0.5 * event.weight();

      // Fill jet number integral histograms
      _h_excmult_jets_tot->fill(goodjets.size(), w);
      /// @todo Could be better computed by toIntegral transform on exclusive histo
      for (size_t iJet = 1; iJet <= goodjets.size(); iJet++ )
        _h_incmult_jets_tot->fill(iJet, w);

      // Fill leading jet histograms
      const Jet& j1 = goodjets[0];
      _h_leading_jet_pt_tot->fill(j1.pT()/GeV, w);
      _h_leading_jet_eta_tot->fill(j1.abseta(), w);
      _h_ht1_tot->fill(ht/GeV, w);

      // Fill 2nd jet histograms
      if (goodjets.size() < 2) return;
      const Jet& j2 = goodjets[1];
      _h_second_jet_pt_tot->fill(j2.pT()/GeV, w);
      _h_second_jet_eta_tot->fill(j2.abseta(), w);
      _h_ht2_tot->fill(ht/GeV, w);

      // Fill 3rd jet histograms
      if (goodjets.size() < 3) return;
      const Jet& j3 = goodjets[2];
      _h_third_jet_pt_tot->fill(j3.pT()/GeV, w);
      _h_third_jet_eta_tot->fill(j3.abseta(), w);
      _h_ht3_tot->fill(ht/GeV, w);

      // Fill 4th jet histograms
      if (goodjets.size() < 4) return;
      const Jet& j4 = goodjets[3];
      _h_fourth_jet_pt_tot->fill(j4.pT()/GeV, w);
      _h_fourth_jet_eta_tot->fill(j4.abseta(), w);
      _h_ht4_tot->fill(ht/GeV, w);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = (sumOfWeights() != 0) ? crossSection()/sumOfWeights() : 1.0;

      MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
      MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
      MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());
      MSG_INFO("Norm factor   = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(6) << norm);

      scale(_h_excmult_jets_tot, norm );
      scale(_h_incmult_jets_tot, norm );
      scale(_h_leading_jet_pt_tot, norm );
      scale(_h_second_jet_pt_tot, norm );
      scale(_h_third_jet_pt_tot, norm );
      scale(_h_fourth_jet_pt_tot, norm );
      scale(_h_leading_jet_eta_tot, norm );
      scale(_h_second_jet_eta_tot, norm );
      scale(_h_third_jet_eta_tot, norm );
      scale(_h_fourth_jet_eta_tot, norm );
      scale(_h_ht1_tot, norm );
      scale(_h_ht2_tot, norm );
      scale(_h_ht3_tot, norm );
      scale(_h_ht4_tot, norm );
    }


  private:

    /// @name Histograms

    Histo1DPtr _h_excmult_jets_tot,  _h_incmult_jets_tot;
    Histo1DPtr _h_leading_jet_pt_tot, _h_second_jet_pt_tot, _h_third_jet_pt_tot, _h_fourth_jet_pt_tot;
    Histo1DPtr _h_leading_jet_eta_tot, _h_second_jet_eta_tot, _h_third_jet_eta_tot, _h_fourth_jet_eta_tot;
    Histo1DPtr _h_ht1_tot, _h_ht2_tot, _h_ht3_tot, _h_ht4_tot;

  };


  DECLARE_RIVET_PLUGIN(CMS_2015_I1310737);


}
#line 1 "CMS_2015_I1327224.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CMS_2015_I1327224 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2015_I1327224);


    void init() {
      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");
      _h_chi_dijet.addHistogram(4200., 8000., bookHisto1D(1, 1, 1));
      _h_chi_dijet.addHistogram(3600., 4200., bookHisto1D(2, 1, 1));
      _h_chi_dijet.addHistogram(3000., 3600., bookHisto1D(3, 1, 1));
      _h_chi_dijet.addHistogram(2400., 3000., bookHisto1D(4, 1, 1));
      _h_chi_dijet.addHistogram(1900., 2400., bookHisto1D(5, 1, 1));
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      FourMomentum j0(jets[0].momentum());
      FourMomentum j1(jets[1].momentum());
      double y0 = j0.rapidity();
      double y1 = j1.rapidity();
      if (fabs(y0 + y1) / 2. > 1.11) vetoEvent;

      double mjj = FourMomentum(j0 + j1).mass();
      if (mjj/GeV <1900) vetoEvent;

      double chi = exp(fabs(y0 - y1));
      if (chi >= 16.) vetoEvent;

      // Fill the histogram
      _h_chi_dijet.fill(mjj/GeV, chi, weight);
    }

    void finalize() {
      foreach (Histo1DPtr hist, _h_chi_dijet.getHistograms()) {
        normalize(hist);
      }
    }

  private:
    BinnedHistogram<double> _h_chi_dijet;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1327224);
}
#line 1 "CMS_2015_I1346843.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// Differential cross-section of FSR photons in Z decays
  class CMS_2015_I1346843 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2015_I1346843);

    /// Book histograms and initialise projections before the run
    void init() {

      Cut c_photons = Cuts::pT >= 5.0*GeV && (Cuts::etaIn(-2.5, 1.4) || Cuts::etaIn(1.6, 2.5));
      IdentifiedFinalState photons(c_photons);
      photons.acceptId(PID::PHOTON);
      declare(photons, "PHOTFS");

      Cut c_muons   = Cuts::pT > 9*GeV && Cuts::abseta < 2.4;
      IdentifiedFinalState muons(c_muons);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "MUFS");


      _hist_pho_et           = bookHisto1D(1, 1, 1);  // photon transverse energy
      _hist_pho_et_wide      = bookHisto1D(1, 2, 1);  // photon transverse energy (0.5 < dr < 3.0)
      _hist_pho_et_close     = bookHisto1D(1, 3, 1);  // photon transverse energy (0.05 < dr < 0.5)
      _hist_pho_et_lqt       = bookHisto1D(1, 4, 1);  // photon transverse energy (q_T < 10)
      _hist_pho_et_hqt       = bookHisto1D(1, 5, 1);  // photon transverse energy (q_T > 50)
      _hist_pho_dr           = bookHisto1D(2, 1, 1);  // delta_R
      _hist_pho_dr_lqt       = bookHisto1D(2, 2, 1);  // delta_R (q_T < 10)
      _hist_pho_dr_hqt       = bookHisto1D(2, 3, 1);  // delta_R  (q_T > 50)
    }


    // Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles muons = apply<IdentifiedFinalState>(event, "MUFS").particlesByPt();

      if (muons.size() < 2) vetoEvent;
      if (muons[0].pT()/GeV < 31) vetoEvent;
      if (muons[0].charge()*muons[1].charge() > 0) vetoEvent;
      const double mZ = (muons[0].momentum() + muons[1].momentum()).mass();
      if (!inRange(mZ, 30*GeV, 87*GeV)) vetoEvent;

      const Particles photons = apply<IdentifiedFinalState>(event, "PHOTFS").particlesByPt();
      // We want the photon with the highest pT that does not come from a decay
      foreach(const Particle& p, photons) {
        if (p.fromDecay() || !p.isStable()) continue;

        const double dR = std::min(deltaR(p, muons[0]), deltaR(p, muons[1]) );
        if (!inRange(dR, 0.05, 3.0)) continue;

        // Calculate the three-body (mu,mu,gamma) transverse momentum
        const double qT = (muons[0].mom() + muons[1].mom() + p.mom()).pT();

        // Fill the analysis histograms
        _hist_pho_et->fill(p.pT()/GeV, event.weight());
        _hist_pho_dr->fill(dR, event.weight());

        (dR <= 0.5 ? _hist_pho_et_close : _hist_pho_et_wide)->fill(p.pT()/GeV, event.weight());

        if (qT / GeV < 10.) {
          _hist_pho_et_lqt->fill(p.pT()/GeV, event.weight());
          _hist_pho_dr_lqt->fill(dR, event.weight());
        }

        if (qT / GeV > 50.) {
          _hist_pho_et_hqt->fill(p.pT()/GeV, event.weight());
          _hist_pho_dr_hqt->fill(dR, event.weight());
        }

        break; // Exit the loop since we found the highest pT lepton already
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_pho_et,       crossSection() / sumOfWeights());
      scale(_hist_pho_et_wide,  crossSection() / sumOfWeights());
      scale(_hist_pho_et_close, crossSection() / sumOfWeights());
      scale(_hist_pho_et_lqt,   crossSection() / sumOfWeights());
      scale(_hist_pho_et_hqt,   crossSection() / sumOfWeights());
      scale(_hist_pho_dr,       crossSection() / sumOfWeights());
      scale(_hist_pho_dr_lqt,   crossSection() / sumOfWeights());
      scale(_hist_pho_dr_hqt,   crossSection() / sumOfWeights());
    }


  private:

    Histo1DPtr _hist_pho_et;
    Histo1DPtr _hist_pho_et_wide, _hist_pho_et_close;
    Histo1DPtr _hist_pho_et_lqt,  _hist_pho_et_hqt;
    Histo1DPtr _hist_pho_dr;
    Histo1DPtr _hist_pho_dr_lqt, _hist_pho_dr_hqt;

  };


  DECLARE_RIVET_PLUGIN(CMS_2015_I1346843);

}
#line 1 "CMS_2015_I1356998.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class CMS_2015_I1356998 : public Analysis {
  public:

    CMS_2015_I1356998()
      : Analysis("CMS_2015_I1356998"), edge(4.7)
    {    }


    void init() {

      declare(FinalState(),"FS");

      _h_noCASTORtag = bookHisto1D(1, 1, 1);
      _h_CASTORtag   = bookHisto1D(2, 1, 1);
      _h_centralGap  = bookHisto1D(3, 1, 1);
      _h_sigmaVis    = bookHisto1D(4, 1, 1);
      _h_maxFwdGap   = bookHisto1D(5, 1, 1);

    }


    void analyze(const Event& event) {

      const double weight = event.weight();
      const FinalState& fs = apply<FinalState>(event, "FS");

      // A vector containing a lot of eta values
      vector<double> detparticles;
      detparticles.push_back(-edge);
      foreach (const Particle& p, fs.particles(Cuts::pT > 0.2*GeV && Cuts::abseta<edge, cmpMomByEta) ) {
        detparticles.push_back(p.momentum().eta());
      }
      detparticles.push_back(edge);

      // Find maximum gap size
      vector <double>::iterator iter;
      vector<double> detgaps;
      for (iter = detparticles.begin()+1; iter != detparticles.end(); ++iter) {
        const double detgap = *iter - *(iter-1);
        detgaps.push_back(detgap);
      }
      double detgapbwd = detgaps.front();
      double detgapfwd = detgaps.back();
      double detfmax = max(detgapbwd, detgapfwd);

      // Fill rapidity gap histo
      if (detfmax != 2*edge ) {
        _h_maxFwdGap->fill(detfmax, weight);
      }
      // Everything that follows has to do with the cross-section measurements

      if (fs.size() < 2) vetoEvent;

      // Gap center calculations
      const ParticleVector particlesByRapidity = fs.particles(cmpMomByRap); //ByRapidity();

      vector<double> gaps;
      vector<double> midpoints;
      for (size_t ip = 1; ip < particlesByRapidity.size(); ++ip) {
        const Particle& p1 = particlesByRapidity[ip-1];
        const Particle& p2 = particlesByRapidity[ip];
        const double gap = p2.momentum().rapidity()  - p1.momentum().rapidity();
        const double mid = (p2.momentum().rapidity() + p1.momentum().rapidity()) / 2.;
        gaps.push_back(gap);
        midpoints.push_back(mid);
      }

      int imid = std::distance(gaps.begin(), max_element(gaps.begin(), gaps.end()));
      double gapcenter = midpoints[imid];

      // Calculations for cross-sections
      FourMomentum MxFourVector(0.,0.,0.,0.);
      FourMomentum MyFourVector(0.,0.,0.,0.);

      foreach(const Particle& p, fs.particles(cmpMomByEta)) {
        if (p.momentum().rapidity() > gapcenter) {
          MxFourVector += p.momentum();
        }
        else {
          MyFourVector += p.momentum();
        }
      }

      double Mx = MxFourVector.mass();
      double My = MyFourVector.mass();

      const double xix = (Mx*Mx)/(sqrtS()/GeV * sqrtS()/GeV);

      if (log10(My) < 0.5) {
        _h_noCASTORtag->fill(log10(xix), weight);
        if (log10(xix) > -5.5 && log10(xix) < -2.5) _h_sigmaVis->fill(0.5, weight);
      }
      else if (log10(My) < 1.1) {
        _h_CASTORtag->fill(log10(xix), weight);
        if (log10(xix) > -5.5 && log10(xix) < -2.5) _h_sigmaVis->fill(1.5, weight);
      }

      // Central gap x-section
      double xigen = (Mx*Mx) * (My*My) / (sqrtS()/GeV * sqrtS()/GeV * 0.93827 * 0.93827); // Proton masses...
      double dy0 = -log(xigen);

      if (dy0 > 3.) {
        if (log10(My) > 1.1 && log10(Mx) > 1.1) {
          _h_centralGap->fill(dy0, weight);
          _h_sigmaVis->fill(2.5, weight);
        }
      }

    }

    void finalize() {

      double xs = crossSection()/millibarn/sumOfWeights();
      scale(_h_noCASTORtag, xs);
      scale(_h_CASTORtag  , xs);
      scale(_h_centralGap , xs);
      scale(_h_sigmaVis   , xs);
      scale(_h_maxFwdGap  , xs);

    }

  private:

    Histo1DPtr _h_noCASTORtag;
    Histo1DPtr _h_CASTORtag;
    Histo1DPtr _h_centralGap;
    Histo1DPtr _h_sigmaVis;
    Histo1DPtr _h_maxFwdGap;
    double edge;

  };


  DECLARE_RIVET_PLUGIN(CMS_2015_I1356998);

}
#line 1 "CMS_2015_I1370682.cc"
#include "Rivet/Analysis.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  namespace { //< only visible in this compilation unit

    /// @brief Pseudo top finder
    ///
    /// Find top quark in the particle level.
    /// The definition is based on the agreement at the LHC working group.
    class PseudoTop : public FinalState {
    public:
      /// @name Standard constructors and destructors.
      //@{

      /// The default constructor. May specify the minimum and maximum
      /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
      PseudoTop(double lepR = 0.1, double lepMinPt = 20, double lepMaxEta = 2.4,
                double jetR = 0.4, double jetMinPt = 30, double jetMaxEta = 4.7)
        : FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV),
          _lepR(lepR), _lepMinPt(lepMinPt), _lepMaxEta(lepMaxEta),
          _jetR(jetR), _jetMinPt(jetMinPt), _jetMaxEta(jetMaxEta)
      {
        setName("PseudoTop");
      }

      enum TTbarMode {CH_NONE=-1, CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON};
      enum DecayMode {CH_HADRON = 0, CH_MUON, CH_ELECTRON};

      TTbarMode mode() const {
        if (!_isValid) return CH_NONE;
        if (_mode1 == CH_HADRON && _mode2 == CH_HADRON) return CH_FULLHADRON;
        else if ( _mode1 != CH_HADRON && _mode2 != CH_HADRON) return CH_FULLLEPTON;
        else return CH_SEMILEPTON;
      }
      DecayMode mode1() const {return _mode1;}
      DecayMode mode2() const {return _mode2;}

      /// Clone on the heap.
      virtual unique_ptr<Projection> clone() const {
        return unique_ptr<Projection>(new PseudoTop(*this));
      }

      //@}

    public:
      Particle t1() const {return _t1;}
      Particle t2() const {return _t2;}
      Particle b1() const {return _b1;}
      Particle b2() const {return _b2;}
      ParticleVector wDecays1() const {return _wDecays1;}
      ParticleVector wDecays2() const {return _wDecays2;}
      Jets jets() const {return _jets;}
      Jets bjets() const {return _bjets;}
      Jets ljets() const {return _ljets;}

    protected:
      // Apply the projection to the event
      void project(const Event& e); // override; ///< @todo Re-enable when C++11 allowed
      void cleanup(std::map<double, std::pair<size_t, size_t> >& v, const bool doCrossCleanup=false) const;

    private:
      const double _lepR, _lepMinPt, _lepMaxEta;
      const double _jetR, _jetMinPt, _jetMaxEta;

      //constexpr ///< @todo Re-enable when C++11 allowed
      static double _tMass; // = 172.5*GeV; ///< @todo Re-enable when C++11 allowed
      //constexpr ///< @todo Re-enable when C++11 allowed
      static double _wMass; // = 80.4*GeV; ///< @todo Re-enable when C++11 allowed

    private:
      bool _isValid;
      DecayMode _mode1, _mode2;

      Particle _t1, _t2;
      Particle _b1, _b2;
      ParticleVector _wDecays1, _wDecays2;
      Jets _jets, _bjets, _ljets;

    };

    // More implementation below the analysis code

  }



  /// Pseudo-top analysis from CMS
  class CMS_2015_I1370682 : public Analysis {
  public:

    CMS_2015_I1370682()
      : Analysis("CMS_2015_I1370682"),
        _applyCorrection(true),
        _doShapeOnly(false)
    {    }


    void init() {
      declare(PseudoTop(0.1, 20, 2.4, 0.5, 30, 2.4), "ttbar");

      // Lepton + Jet channel
      _hSL_topPt         = bookHisto1D("d15-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hSL_topPtTtbarSys = bookHisto1D("d16-x01-y01"); // 1/sigma dsigma/dpt*(top)
      _hSL_topY          = bookHisto1D("d17-x01-y01"); // 1/sigma dsigma/dy(top)
      _hSL_ttbarDelPhi   = bookHisto1D("d18-x01-y01"); // 1/sigma dsigma/ddeltaphi(t,tbar)
      _hSL_topPtLead     = bookHisto1D("d19-x01-y01"); // 1/sigma dsigma/dpt(t1)
      _hSL_topPtSubLead  = bookHisto1D("d20-x01-y01"); // 1/sigma dsigma/dpt(t2)
      _hSL_ttbarPt       = bookHisto1D("d21-x01-y01"); // 1/sigma dsigma/dpt(ttbar)
      _hSL_ttbarY        = bookHisto1D("d22-x01-y01"); // 1/sigma dsigma/dy(ttbar)
      _hSL_ttbarMass     = bookHisto1D("d23-x01-y01"); // 1/sigma dsigma/dm(ttbar)

      // Dilepton channel
      _hDL_topPt         = bookHisto1D("d24-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hDL_topPtTtbarSys = bookHisto1D("d25-x01-y01"); // 1/sigma dsigma/dpt*(top)
      _hDL_topY          = bookHisto1D("d26-x01-y01"); // 1/sigma dsigma/dy(top)
      _hDL_ttbarDelPhi   = bookHisto1D("d27-x01-y01"); // 1/sigma dsigma/ddeltaphi(t,tbar)
      _hDL_topPtLead     = bookHisto1D("d28-x01-y01"); // 1/sigma dsigma/dpt(t1)
      _hDL_topPtSubLead  = bookHisto1D("d29-x01-y01"); // 1/sigma dsigma/dpt(t2)
      _hDL_ttbarPt       = bookHisto1D("d30-x01-y01"); // 1/sigma dsigma/dpt(ttbar)
      _hDL_ttbarY        = bookHisto1D("d31-x01-y01"); // 1/sigma dsigma/dy(ttbar)
      _hDL_ttbarMass     = bookHisto1D("d32-x01-y01"); // 1/sigma dsigma/dm(ttbar)

    }


    void analyze(const Event& event) {

      // Get the ttbar candidate
      const PseudoTop& ttbar = apply<PseudoTop>(event, "ttbar");
      if ( ttbar.mode() == PseudoTop::CH_NONE ) vetoEvent;

      const FourMomentum& t1P4 = ttbar.t1().momentum();
      const FourMomentum& t2P4 = ttbar.t2().momentum();
      const double pt1 = std::max(t1P4.pT(), t2P4.pT());
      const double pt2 = std::min(t1P4.pT(), t2P4.pT());
      const double dPhi = deltaPhi(t1P4, t2P4);
      const FourMomentum ttP4 = t1P4 + t2P4;
      const FourMomentum t1P4AtCM = LorentzTransform::mkFrameTransformFromBeta(ttP4.betaVec()).transform(t1P4);

      const double weight = event.weight();

      if ( ttbar.mode() == PseudoTop::CH_SEMILEPTON ) {
        const Particle lCand1 = ttbar.wDecays1()[0]; // w1 dau0 is the lepton in the PseudoTop
        if (lCand1.pT() < 33*GeV || lCand1.abseta() > 2.1) vetoEvent;
        _hSL_topPt->fill(t1P4.pT(), weight);
        _hSL_topPt->fill(t2P4.pT(), weight);
        _hSL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
        _hSL_topY->fill(t1P4.rapidity(), weight);
        _hSL_topY->fill(t2P4.rapidity(), weight);
        _hSL_ttbarDelPhi->fill(dPhi, weight);
        _hSL_topPtLead->fill(pt1, weight);
        _hSL_topPtSubLead->fill(pt2, weight);
        _hSL_ttbarPt->fill(ttP4.pT(), weight);
        _hSL_ttbarY->fill(ttP4.rapidity(), weight);
        _hSL_ttbarMass->fill(ttP4.mass(), weight);
      }
      else if ( ttbar.mode() == PseudoTop::CH_FULLLEPTON ) {
        const Particle lCand1 = ttbar.wDecays1()[0]; // dau0 are the lepton in the PseudoTop
        const Particle lCand2 = ttbar.wDecays2()[0]; // dau0 are the lepton in the PseudoTop
        if (lCand1.pT() < 20*GeV || lCand1.abseta() > 2.4) vetoEvent;
        if (lCand2.pT() < 20*GeV || lCand2.abseta() > 2.4) vetoEvent;
        _hDL_topPt->fill(t1P4.pT(), weight);
        _hDL_topPt->fill(t2P4.pT(), weight);
        _hDL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
        _hDL_topY->fill(t1P4.rapidity(), weight);
        _hDL_topY->fill(t2P4.rapidity(), weight);
        _hDL_ttbarDelPhi->fill(dPhi, weight);
        _hDL_topPtLead->fill(pt1, weight);
        _hDL_topPtSubLead->fill(pt2, weight);
        _hDL_ttbarPt->fill(ttP4.pT(), weight);
        _hDL_ttbarY->fill(ttP4.rapidity(), weight);
        _hDL_ttbarMass->fill(ttP4.mass(), weight);
      }

    }


    void finalize() {
      if ( _applyCorrection ) {
        // Correction functions for TOP-12-028 paper, (parton bin height)/(pseudotop bin height)
        const double ch15[] = { 5.473609, 4.941048, 4.173346, 3.391191, 2.785644, 2.371346, 2.194161, 2.197167, };
        const double ch16[] = { 5.470905, 4.948201, 4.081982, 3.225532, 2.617519, 2.239217, 2.127878, 2.185918, };
        const double ch17[] = { 10.003667, 4.546519, 3.828115, 3.601018, 3.522194, 3.524694, 3.600951, 3.808553, 4.531891, 9.995370, };
        const double ch18[] = { 4.406683, 4.054041, 3.885393, 4.213646, };
        const double ch19[] = { 6.182537, 5.257703, 4.422280, 3.568402, 2.889408, 2.415878, 2.189974, 2.173210, };
        const double ch20[] = { 5.199874, 4.693318, 3.902882, 3.143785, 2.607877, 2.280189, 2.204124, 2.260829, };
        const double ch21[] = { 6.053523, 3.777506, 3.562251, 3.601356, 3.569347, 3.410472, };
        const double ch22[] = { 11.932351, 4.803773, 3.782709, 3.390775, 3.226806, 3.218982, 3.382678, 3.773653, 4.788191, 11.905338, };
        const double ch23[] = { 7.145255, 5.637595, 4.049882, 3.025917, 2.326430, 1.773824, 1.235329, };

        const double ch24[] = { 2.268193, 2.372063, 2.323975, 2.034655, 1.736793, };
        const double ch25[] = { 2.231852, 2.383086, 2.341894, 2.031318, 1.729672, 1.486993, };
        const double ch26[] = { 3.993526, 2.308249, 2.075136, 2.038297, 2.036302, 2.078270, 2.295817, 4.017713, };
        const double ch27[] = { 2.205978, 2.175010, 2.215376, 2.473144, };
        const double ch28[] = { 2.321077, 2.371895, 2.338871, 2.057821, 1.755382, };
        const double ch29[] = { 2.222707, 2.372591, 2.301688, 1.991162, 1.695343, };
        const double ch30[] = { 2.599677, 2.026855, 2.138620, 2.229553, };
        const double ch31[] = { 5.791779, 2.636219, 2.103642, 1.967198, 1.962168, 2.096514, 2.641189, 5.780828, };
        const double ch32[] = { 2.006685, 2.545525, 2.477745, 2.335747, 2.194226, 2.076500, };

        applyCorrection(_hSL_topPt, ch15);
        applyCorrection(_hSL_topPtTtbarSys, ch16);
        applyCorrection(_hSL_topY, ch17);
        applyCorrection(_hSL_ttbarDelPhi, ch18);
        applyCorrection(_hSL_topPtLead, ch19);
        applyCorrection(_hSL_topPtSubLead, ch20);
        applyCorrection(_hSL_ttbarPt, ch21);
        applyCorrection(_hSL_ttbarY, ch22);
        applyCorrection(_hSL_ttbarMass, ch23);

        applyCorrection(_hDL_topPt, ch24);
        applyCorrection(_hDL_topPtTtbarSys, ch25);
        applyCorrection(_hDL_topY, ch26);
        applyCorrection(_hDL_ttbarDelPhi, ch27);
        applyCorrection(_hDL_topPtLead, ch28);
        applyCorrection(_hDL_topPtSubLead, ch29);
        applyCorrection(_hDL_ttbarPt, ch30);
        applyCorrection(_hDL_ttbarY, ch31);
        applyCorrection(_hDL_ttbarMass, ch32);
      }

      if ( _doShapeOnly ) {
        normalize(_hSL_topPt        );
        normalize(_hSL_topPtTtbarSys);
        normalize(_hSL_topY         );
        normalize(_hSL_ttbarDelPhi  );
        normalize(_hSL_topPtLead    );
        normalize(_hSL_topPtSubLead );
        normalize(_hSL_ttbarPt      );
        normalize(_hSL_ttbarY       );
        normalize(_hSL_ttbarMass    );

        normalize(_hDL_topPt        );
        normalize(_hDL_topPtTtbarSys);
        normalize(_hDL_topY         );
        normalize(_hDL_ttbarDelPhi  );
        normalize(_hDL_topPtLead    );
        normalize(_hDL_topPtSubLead );
        normalize(_hDL_ttbarPt      );
        normalize(_hDL_ttbarY       );
        normalize(_hDL_ttbarMass    );
      }
      else {
        const double s = 1./sumOfWeights();
        scale(_hSL_topPt        , s);
        scale(_hSL_topPtTtbarSys, s);
        scale(_hSL_topY         , s);
        scale(_hSL_ttbarDelPhi  , s);
        scale(_hSL_topPtLead    , s);
        scale(_hSL_topPtSubLead , s);
        scale(_hSL_ttbarPt      , s);
        scale(_hSL_ttbarY       , s);
        scale(_hSL_ttbarMass    , s);
        scale(_hDL_topPt        , s);
        scale(_hDL_topPtTtbarSys, s);
        scale(_hDL_topY         , s);
        scale(_hDL_ttbarDelPhi  , s);
        scale(_hDL_topPtLead    , s);
        scale(_hDL_topPtSubLead , s);
        scale(_hDL_ttbarPt      , s);
        scale(_hDL_ttbarY       , s);
        scale(_hDL_ttbarMass    , s);
      }

    }


    void applyCorrection(Histo1DPtr h, const double* cf) {
      vector<YODA::HistoBin1D>& bins = h->bins();
      for (size_t i=0, n=bins.size(); i<n; ++i ) {
        const double s = cf[i];
        YODA::HistoBin1D& bin = bins[i];
        bin.scaleW(s);
      }
    }


  private:

    const bool _applyCorrection, _doShapeOnly;
    Histo1DPtr _hSL_topPt, _hSL_topPtTtbarSys, _hSL_topY, _hSL_ttbarDelPhi, _hSL_topPtLead,
      _hSL_topPtSubLead, _hSL_ttbarPt, _hSL_ttbarY, _hSL_ttbarMass;
    Histo1DPtr _hDL_topPt, _hDL_topPtTtbarSys, _hDL_topY, _hDL_ttbarDelPhi, _hDL_topPtLead,
      _hDL_topPtSubLead, _hDL_ttbarPt, _hDL_ttbarY, _hDL_ttbarMass;

  };



  DECLARE_RIVET_PLUGIN(CMS_2015_I1370682);


  ///////////////

  // More PseudoTop implementation
  namespace {


    double PseudoTop::_tMass = 172.5*GeV;
    double PseudoTop::_wMass = 80.4*GeV;


    void PseudoTop::cleanup(map<double, pair<size_t, size_t> >& v, const bool doCrossCleanup) const {
      vector<map<double, pair<size_t, size_t> >::iterator> toErase;
      set<size_t> usedLeg1, usedLeg2;
      if ( !doCrossCleanup ) {
        /// @todo Reinstate when C++11 allowed: for (auto key = v.begin(); key != v.end(); ++key) {
        for (map<double, pair<size_t, size_t> >::iterator key = v.begin(); key != v.end(); ++key) {
          const size_t leg1 = key->second.first;
          const size_t leg2 = key->second.second;
          if (usedLeg1.find(leg1) == usedLeg1.end() and
              usedLeg2.find(leg2) == usedLeg2.end()) {
            usedLeg1.insert(leg1);
            usedLeg2.insert(leg2);
          } else {
            toErase.push_back(key);
          }
        }
      }
      else {
        /// @todo Reinstate when C++11 allowed: for (auto key = v.begin(); key != v.end(); ++key) {
        for (map<double, pair<size_t, size_t> >::iterator key = v.begin(); key != v.end(); ++key) {
          const size_t leg1 = key->second.first;
          const size_t leg2 = key->second.second;
          if (usedLeg1.find(leg1) == usedLeg1.end() and
              usedLeg1.find(leg2) == usedLeg1.end()) {
            usedLeg1.insert(leg1);
            usedLeg1.insert(leg2);
          } else {
            toErase.push_back(key);
          }
        }
      }
      /// @todo Reinstate when C++11 allowed:  for (auto& key : toErase) v.erase(key);
      for (size_t i = 0; i < toErase.size(); ++i) v.erase(toErase[i]);
    }


    void PseudoTop::project(const Event& e) {
      // Leptons : do the lepton clustering anti-kt R=0.1 using stable photons and leptons not from hadron decay
      // Neutrinos : neutrinos not from hadron decay
      // MET : vector sum of all invisible particles in x-y plane
      // Jets : anti-kt R=0.4 using all particles excluding neutrinos and particles used in lepton clustering
      //        add ghost B hadrons during the jet clustering to identify B jets.

      // W->lv : dressed lepton and neutrino pairs
      // W->jj : light flavored dijet
      // W candidate : select lv or jj pairs which minimise |mW1-80.4|+|mW2-80.4|
      //               lepton-neutrino pair will be selected with higher priority

      // t->Wb : W candidate + b jet
      // t candidate : select Wb pairs which minimise |mtop1-172.5|+|mtop2-172.5|

      _isValid = false;
      _theParticles.clear();
      _wDecays1.clear();
      _wDecays2.clear();
      _jets.clear();
      _bjets.clear();
      _ljets.clear();
      _mode1 = _mode2 = CH_HADRON;

      // Collect final state particles
      Particles pForLep, pForJet;
      Particles neutrinos; // Prompt neutrinos
      /// @todo Avoid this unsafe jump into HepMC -- all this can be done properly via VisibleFS and HeavyHadrons projections
      for (const GenParticle* p : Rivet::particles(e.genEvent())) {
        const int status = p->status();
        const int pdgId = p->pdg_id();
        if (status == 1) {
          Particle rp = *p;
          if (!PID::isHadron(pdgId) && !rp.fromHadron()) {
            // Collect particles not from hadron decay
            if (rp.isNeutrino()) {
              // Prompt neutrinos are kept in separate collection
              neutrinos.push_back(rp);
            } else if (pdgId == 22 || rp.isLepton()) {
              // Leptons and photons for the dressing
              pForLep.push_back(rp);
            }
          } else if (!rp.isNeutrino()) {
            // Use all particles from hadron decay
            pForJet.push_back(rp);
          }
        } else if (PID::isHadron(pdgId) && PID::hasBottom(pdgId)) {
          // NOTE: Consider B hadrons with pT > 5GeV - not in CMS proposal
          //if ( p->momentum().perp() < 5 ) continue;

          // Do unstable particles, to be used in the ghost B clustering
          // Use last B hadrons only
          bool isLast = true;
          for (GenParticle* pp : Rivet::particles(p->end_vertex(), HepMC::children)) {
            if (PID::hasBottom(pp->pdg_id())) {
              isLast = false;
              break;
            }
          }
          if (!isLast) continue;

          // Rescale momentum by 10^-20
          Particle ghost(pdgId, FourMomentum(p->momentum())*1e-20/p->momentum().rho());
          pForJet.push_back(ghost);
        }
      }

      // Start object building from trivial thing - prompt neutrinos
      sortByPt(neutrinos);

      // Proceed to lepton dressing
      const PseudoJets lep_pjs = mkPseudoJets(pForLep);
      const fastjet::JetDefinition lep_jdef(fastjet::antikt_algorithm, _lepR);
      const Jets leps_all = mkJets(fastjet::ClusterSequence(lep_pjs, lep_jdef).inclusive_jets());
      const Jets leps_sel = sortByPt(filterBy(leps_all, Cuts::pT > _lepMinPt));
      // FastJets fjLep(FastJets::ANTIKT, _lepR);
      // fjLep.calc(pForLep);

      Jets leptons;
      vector<int> leptonsId;
      set<int> dressedIdxs;
      for (const Jet& lep : leps_sel) {
        if (lep.abseta() > _lepMaxEta) continue;
        double leadingPt = -1;
        int leptonId = 0;
        for (const Particle& p : lep.particles()) {
          /// @warning Barcodes aren't future-proof in HepMC
          dressedIdxs.insert(p.genParticle()->barcode());
          if (p.isLepton() && p.pT() > leadingPt) {
            leadingPt = p.pT();
            leptonId = p.pid();
          }
        }
        if (leptonId == 0) continue;
        leptons.push_back(lep);
        leptonsId.push_back(leptonId);
      }

      // Re-use particles not used in lepton dressing
      for (const Particle& rp : pForLep) {
        /// @warning Barcodes aren't future-proof in HepMC
        const int barcode = rp.genParticle()->barcode();
        // Skip if the particle is used in dressing
        if (dressedIdxs.find(barcode) != dressedIdxs.end()) continue;
        // Put back to be used in jet clustering
        pForJet.push_back(rp);
      }

      // Then do the jet clustering
      const PseudoJets jet_pjs = mkPseudoJets(pForJet);
      const fastjet::JetDefinition jet_jdef(fastjet::antikt_algorithm, _jetR);
      const Jets jets_all = mkJets(fastjet::ClusterSequence(jet_pjs, jet_jdef).inclusive_jets());
      const Jets jets_sel = sortByPt(filterBy(jets_all, Cuts::pT > _jetMinPt));
      // FastJets fjJet(FastJets::ANTIKT, _jetR);
      //fjJet.useInvisibles(); // NOTE: CMS proposal to remove neutrinos (AB: wouldn't work anyway, since they were excluded from clustering inputs)
      // fjJet.calc(pForJet);
      for (const Jet& jet : jets_sel) {
        if (jet.abseta() > _jetMaxEta) continue;
        _jets.push_back(jet);
        bool isBJet = false;
        for (const Particle& rp : jet.particles()) {
          if (PID::hasBottom(rp.pdgId())) {
            isBJet = true;
            break;
          }
        }
        if ( isBJet ) _bjets.push_back(jet);
        else _ljets.push_back(jet);
      }

      // Every building blocks are ready. Continue to pseudo-W and pseudo-top combination

      if (_bjets.size() < 2) return; // Ignore single top for now
      map<double, pair<size_t, size_t> > wLepCandIdxs;
      map<double, pair<size_t, size_t> > wHadCandIdxs;

      // Collect leptonic-decaying W's
      for (size_t iLep = 0, nLep = leptons.size(); iLep < nLep; ++iLep) {
        const Jet& lep = leptons.at(iLep);
        for (size_t iNu = 0, nNu = neutrinos.size(); iNu < nNu; ++iNu) {
          const Particle& nu = neutrinos.at(iNu);
          const double m = (lep.momentum()+nu.momentum()).mass();
          const double dm = std::abs(m-_wMass);
          wLepCandIdxs[dm] = make_pair(iLep, iNu);
        }
      }

      // Continue to hadronic decaying W's
      for (size_t i = 0, nLjet = _ljets.size(); i < nLjet; ++i) {
        const Jet& ljet1 = _ljets[i];
        for (size_t j = i+1; j < nLjet; ++j) {
          const Jet& ljet2 = _ljets[j];
          const double m = (ljet1.momentum()+ljet2.momentum()).mass();
          const double dm = std::abs(m-_wMass);
          wHadCandIdxs[dm] = make_pair(i, j);
        }
      }

      // Cleanup W candidate, choose pairs with minimum dm if they share decay products
      cleanup(wLepCandIdxs);
      cleanup(wHadCandIdxs, true);
      const size_t nWLepCand = wLepCandIdxs.size();
      const size_t nWHadCand = wHadCandIdxs.size();

      if (nWLepCand + nWHadCand < 2) return; // We skip single top

      int w1Q = 1, w2Q = -1;
      int w1dau1Id = 1, w2dau1Id = -1;
      FourMomentum w1dau1LVec, w1dau2LVec;
      FourMomentum w2dau1LVec, w2dau2LVec;
      if (nWLepCand == 0) { // Full hadronic case
        const pair<size_t, size_t>& idPair1 = wHadCandIdxs.begin()->second;
        const pair<size_t, size_t>& idPair2 = (++wHadCandIdxs.begin())->second;  ///< @todo Reinstate std::next
        const Jet& w1dau1 = _ljets[idPair1.first];
        const Jet& w1dau2 = _ljets[idPair1.second];
        const Jet& w2dau1 = _ljets[idPair2.first];
        const Jet& w2dau2 = _ljets[idPair2.second];
        w1dau1LVec = w1dau1.momentum();
        w1dau2LVec = w1dau2.momentum();
        w2dau1LVec = w2dau1.momentum();
        w2dau2LVec = w2dau2.momentum();
      } else if (nWLepCand == 1) { // Semi-leptonic case
        const pair<size_t, size_t>& idPair1 = wLepCandIdxs.begin()->second;
        const pair<size_t, size_t>& idPair2 = wHadCandIdxs.begin()->second;
        const Jet& w1dau1 = leptons[idPair1.first];
        const Particle& w1dau2 = neutrinos[idPair1.second];
        const Jet& w2dau1 = _ljets[idPair2.first];
        const Jet& w2dau2 = _ljets[idPair2.second];
        w1dau1LVec = w1dau1.momentum();
        w1dau2LVec = w1dau2.momentum();
        w2dau1LVec = w2dau1.momentum();
        w2dau2LVec = w2dau2.momentum();
        w1dau1Id = leptonsId[idPair1.first];
        w1Q = w1dau1Id > 0 ? -1 : 1;
        w2Q = -w1Q;
        switch (w1dau1Id) {
        case 13: case -13: _mode1 = CH_MUON; break;
        case 11: case -11: _mode1 = CH_ELECTRON; break;
        }
      } else { // Full leptonic case
        const pair<size_t, size_t>& idPair1 = wLepCandIdxs.begin()->second;
        const pair<size_t, size_t>& idPair2 = (++wLepCandIdxs.begin())->second;  ///< @todo Reinstate std::next
        const Jet& w1dau1 = leptons[idPair1.first];
        const Particle& w1dau2 = neutrinos[idPair1.second];
        const Jet& w2dau1 = leptons[idPair2.first];
        const Particle& w2dau2 = neutrinos[idPair2.second];
        w1dau1LVec = w1dau1.momentum();
        w1dau2LVec = w1dau2.momentum();
        w2dau1LVec = w2dau1.momentum();
        w2dau2LVec = w2dau2.momentum();
        w1dau1Id = leptonsId[idPair1.first];
        w2dau1Id = leptonsId[idPair2.first];
        w1Q = w1dau1Id > 0 ? -1 : 1;
        w2Q = w2dau1Id > 0 ? -1 : 1;
        switch (w1dau1Id) {
        case 13: case -13: _mode1 = CH_MUON; break;
        case 11: case -11: _mode1 = CH_ELECTRON; break;
        }
        switch (w2dau1Id) {
        case 13: case -13: _mode2 = CH_MUON; break;
        case 11: case -11: _mode2 = CH_ELECTRON; break;
        }
      }
      const FourMomentum w1LVec = w1dau1LVec+w1dau2LVec;
      const FourMomentum w2LVec = w2dau1LVec+w2dau2LVec;

      // Combine b jets
      double sumDm = 1e9;
      FourMomentum b1LVec, b2LVec;
      for (size_t i = 0, n = _bjets.size(); i < n; ++i) {
        const Jet& bjet1 = _bjets[i];
        const double mtop1 = (w1LVec+bjet1.momentum()).mass();
        const double dmtop1 = std::abs(mtop1-_tMass);
        for (size_t j=0; j<n; ++j) {
          if (i == j) continue;
          const Jet& bjet2 = _bjets[j];
          const double mtop2 = (w2LVec+bjet2.momentum()).mass();
          const double dmtop2 = std::abs(mtop2-_tMass);

          if (sumDm <= dmtop1+dmtop2) continue;

          sumDm = dmtop1+dmtop2;
          b1LVec = bjet1.momentum();
          b2LVec = bjet2.momentum();
        }
      }
      if (sumDm >= 1e9) return; // Failed to make top, but this should not happen.

      const FourMomentum t1LVec = w1LVec + b1LVec;
      const FourMomentum t2LVec = w2LVec + b2LVec;

      // Put all of them into candidate collection
      _t1 = Particle(w1Q*6, t1LVec);
      _b1 = Particle(w1Q*5, b1LVec);
      _wDecays1.push_back(Particle(w1dau1Id, w1dau1LVec));
      _wDecays1.push_back(Particle(-w1dau1Id+w1Q, w1dau2LVec));

      _t2 = Particle(w2Q*6, t2LVec);
      _b2 = Particle(w2Q*5, b2LVec);
      _wDecays2.push_back(Particle(w2dau1Id, w2dau1LVec));
      _wDecays2.push_back(Particle(-w2dau1Id+w2Q, w2dau2LVec));

      _isValid = true;
    }

  }


}
#line 1 "CMS_2015_I1384119.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class CMS_2015_I1384119 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2015_I1384119);


    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fsa(Cuts::abseta < 20);
      declare(fsa, "FSA");
      const ChargedFinalState cfs(Cuts::abseta < 2);
      declare(cfs, "CFS");

      _hist_dNch_dEta_inel = bookHisto1D(1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Apply inelastic selection (veto pp -> pp elastic events)
      const FinalState& fsa = apply<FinalState>(event, "FSA");
      if (fsa.size() <= 2) vetoEvent;

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      foreach (const Particle& p, cfs.particles()) {
        const int id = p.abspid();
        // continue if particle is a proton, a kaon or a pion
        if (id == 211 || id == 321 || id == 2212) ///< @todo Use PID:: ID constants
          _hist_dNch_dEta_inel->fill(p.eta(), event.weight());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_dNch_dEta_inel,  1/sumOfWeights());
    }


  private:

    /// Histograms
    Histo1DPtr _hist_dNch_dEta_inel;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1384119);

}
#line 1 "CMS_2015_I1385107.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// CMS UE charged particles vs. leading jet at 2.76 TeV
  class CMS_2015_I1385107 : public Analysis {
  public:
    /// Constructor
    CMS_2015_I1385107() : Analysis("CMS_2015_I1385107"),
                          ETACUT(2.0),
                          AREATOT(2*ETACUT * 2*M_PI),
                          AREA3(AREATOT / 3.),
                          AREA6(AREATOT / 6.)
    {   }


    /// Book histograms and initialise projections before the run
    void init() {

      const ChargedFinalState cfs(Cuts::abseta < 2 && Cuts::pT > 500*MeV);
      declare(cfs, "CFS");

      const ChargedFinalState cfsforjet(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV);
      const FastJets jetpro(cfsforjet, FastJets::SISCONE, 0.5);
      declare(jetpro, "Jets");

      _h_Nch_TransAVE_vs_pT = bookProfile1D(1, 1, 1); // Nch vs. pT_max      (TransAVE)
      _h_Sum_TransAVE_vs_pT = bookProfile1D(2, 1, 1); // sum(pT) vs. pT_max  (TransAVE)
      _h_Nch_TransMAX_vs_pT = bookProfile1D(3, 1, 1); // Nch vs. pT_max      (TransMAX)
      _h_Sum_TransMAX_vs_pT = bookProfile1D(4, 1, 1); // sum(pT) vs. pT_max  (TransMAX)
      _h_Nch_TransMIN_vs_pT = bookProfile1D(5, 1, 1); // Nch vs. pT_max      (TransMIN)
      _h_Sum_TransMIN_vs_pT = bookProfile1D(6, 1, 1); // sum(pT) vs. pT_max  (TransMIN)
      _h_Nch_TransDIF_vs_pT = bookProfile1D(7, 1, 1); // Nch vs. pT_max      (TransDIF)
      _h_Sum_TransDIF_vs_pT = bookProfile1D(8, 1, 1); // sum(pT) vs. pT_max  (TransDIF)
    }


    /// Local definition of a signed dphi, for use in differentating L and R trans regions
    double signedDeltaPhi(double jetphi, double partphi) {
      double delta = partphi - jetphi;
      while (delta <= -PI) delta += 2 * PI;
      while (delta > PI) delta -= 2 * PI;
      return delta;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Find the lead jet, applying a restriction that the jets must be within |eta| < 2.
      FourMomentum p_lead;
      foreach (const Jet& j, apply<FastJets>(event, "Jets").jetsByPt(1*GeV)) {
        if (j.abseta() < 2.0) {
          p_lead = j.momentum();
          break;
        }
      }
      if (p_lead.isZero()) vetoEvent;
      const double phi_lead = p_lead.phi();
      const double pT_lead  = p_lead.pT();

      // Loop on charged particles and separate Left and Right transverse regions
      Particles particles = apply<ChargedFinalState>(event, "CFS").particlesByPt();
      int nch_TransLeft = 0, nch_TransRight = 0;
      double ptSum_TransLeft = 0., ptSum_TransRight = 0.;
      foreach (const Particle& p, particles) {
        const double dphi = signedDeltaPhi(phi_lead, p.momentum().phi());
        if (!inRange(fabs(dphi), PI/3, 2*PI/3.)) continue; //< only fill trans regions
        if (dphi < 0) {  // Transverse Right region
          nch_TransRight += 1;
          ptSum_TransRight += p.pT() / GeV;
        } else if (dphi > 0) {  // Transverse Left region
          nch_TransLeft += 1;
          ptSum_TransLeft += p.pT() / GeV;
        }
      }

      // Translate to min and max (+sum and diff) Transverse regions
      const int nch_TransMIN = std::min(nch_TransLeft, nch_TransRight);
      const int nch_TransMAX = std::max(nch_TransLeft, nch_TransRight);
      const int nch_TransSUM = nch_TransMAX + nch_TransMIN;
      const int nch_TransDIF = nch_TransMAX - nch_TransMIN;
      //
      const double ptSum_TransMIN = std::min(ptSum_TransLeft, ptSum_TransRight);
      const double ptSum_TransMAX = std::max(ptSum_TransLeft, ptSum_TransRight);
      const double ptSum_TransSUM = ptSum_TransMAX + ptSum_TransMIN;
      const double ptSum_TransDIF = ptSum_TransMAX - ptSum_TransMIN;

      // Fill profiles
      const double weight = event.weight();
      _h_Nch_TransMIN_vs_pT->fill(pT_lead/GeV, 1/AREA6 * nch_TransMIN, weight);
      _h_Sum_TransMIN_vs_pT->fill(pT_lead/GeV, 1/AREA6 * ptSum_TransMIN, weight);
      //
      _h_Nch_TransMAX_vs_pT->fill(pT_lead/GeV, 1/AREA6 * nch_TransMAX, weight);
      _h_Sum_TransMAX_vs_pT->fill(pT_lead/GeV, 1/AREA6 * ptSum_TransMAX, weight);
      //
      _h_Nch_TransAVE_vs_pT->fill(pT_lead/GeV, 1/AREA3 * nch_TransSUM, weight);
      _h_Sum_TransAVE_vs_pT->fill(pT_lead/GeV, 1/AREA3 * ptSum_TransSUM, weight);
      //
      _h_Nch_TransDIF_vs_pT->fill(pT_lead/GeV, 1/AREA6 * nch_TransDIF, weight);
      _h_Sum_TransDIF_vs_pT->fill(pT_lead/GeV, 1/AREA6 * ptSum_TransDIF, weight);
    }


  private:

    // Data members like post-cuts event weight counters go here
    const double ETACUT, AREATOT, AREA3, AREA6;

    /// Histograms
    Profile1DPtr _h_Nch_TransAVE_vs_pT, _h_Sum_TransAVE_vs_pT;
    Profile1DPtr _h_Nch_TransDIF_vs_pT, _h_Sum_TransDIF_vs_pT;
    Profile1DPtr _h_Nch_TransMIN_vs_pT, _h_Sum_TransMIN_vs_pT;
    Profile1DPtr _h_Nch_TransMAX_vs_pT, _h_Sum_TransMAX_vs_pT;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1385107);

}
#line 1 "CMS_2015_I1397174.cc"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Fully leptonic partonic ttbar analysis
  class CMS_2015_I1397174 : public Analysis {
  public:

    /// Minimal constructor
    CMS_2015_I1397174()
      : Analysis("CMS_2015_I1397174") { }


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {

      // Parton level top quarks
      addProjection(PartonicTops(PartonicTops::E_MU, false), "PartonTops");

      // Find jets not related to the top/W decays
      VetoedFinalState vfs;
      vfs.addDecayProductsVeto(PID::WPLUSBOSON);
      vfs.addDecayProductsVeto(PID::WMINUSBOSON);
      FastJets fj(vfs, FastJets::ANTIKT, 0.5, JetAlg::ALL_MUONS, JetAlg::ALL_INVISIBLES);
      addProjection(fj, "Jets");

      // Book histograms
      _hVis_nJet30_abs       = bookHisto1D( 1, 1, 1);
      _hVis_nJet30           = bookHisto1D( 2, 1, 1);
      _hVis_nJet60_abs       = bookHisto1D( 3, 1, 1);
      _hVis_nJet60           = bookHisto1D( 4, 1, 1);
      _hVis_nJet100_abs      = bookHisto1D( 5, 1, 1);
      _hVis_nJet100          = bookHisto1D( 6, 1, 1);

      _hVis_addJet1Pt_abs    = bookHisto1D( 7, 1, 1);
      _hVis_addJet1Pt        = bookHisto1D( 8, 1, 1);
      _hVis_addJet1Eta_abs   = bookHisto1D( 9, 1, 1);
      _hVis_addJet1Eta       = bookHisto1D(10, 1, 1);
      _hVis_addJet2Pt_abs    = bookHisto1D(11, 1, 1);
      _hVis_addJet2Pt        = bookHisto1D(12, 1, 1);
      _hVis_addJet2Eta_abs   = bookHisto1D(13, 1, 1);
      _hVis_addJet2Eta       = bookHisto1D(14, 1, 1);
      _hVis_addJJMass_abs    = bookHisto1D(15, 1, 1);
      _hVis_addJJMass        = bookHisto1D(16, 1, 1);
      _hVis_addJJDR_abs      = bookHisto1D(17, 1, 1);
      _hVis_addJJDR          = bookHisto1D(18, 1, 1);
      _hVis_addJJHT_abs      = bookHisto1D(19, 1, 1);
      _hVis_addJJHT          = bookHisto1D(20, 1, 1);

      _hFull_addJet1Pt_abs   = bookHisto1D(21, 1, 1);
      _hFull_addJet1Pt       = bookHisto1D(22, 1, 1);
      _hFull_addJet1Eta_abs  = bookHisto1D(23, 1, 1);
      _hFull_addJet1Eta      = bookHisto1D(24, 1, 1);
      _hFull_addJet2Pt_abs   = bookHisto1D(25, 1, 1);
      _hFull_addJet2Pt       = bookHisto1D(26, 1, 1);
      _hFull_addJet2Eta_abs  = bookHisto1D(27, 1, 1);
      _hFull_addJet2Eta      = bookHisto1D(28, 1, 1);
      _hFull_addJJMass_abs   = bookHisto1D(29, 1, 1);
      _hFull_addJJMass       = bookHisto1D(30, 1, 1);
      _hFull_addJJDR_abs     = bookHisto1D(31, 1, 1);
      _hFull_addJJDR         = bookHisto1D(32, 1, 1);
      _hFull_addJJHT_abs     = bookHisto1D(33, 1, 1);
      _hFull_addJJHT         = bookHisto1D(34, 1, 1);

      _hVis_addBJet1Pt_abs   = bookHisto1D(35, 1, 1);
      _hVis_addBJet1Pt       = bookHisto1D(36, 1, 1);
      _hVis_addBJet1Eta_abs  = bookHisto1D(37, 1, 1);
      _hVis_addBJet1Eta      = bookHisto1D(38, 1, 1);
      _hVis_addBJet2Pt_abs   = bookHisto1D(39, 1, 1);
      _hVis_addBJet2Pt       = bookHisto1D(40, 1, 1);
      _hVis_addBJet2Eta_abs  = bookHisto1D(41, 1, 1);
      _hVis_addBJet2Eta      = bookHisto1D(42, 1, 1);
      _hVis_addBBMass_abs    = bookHisto1D(43, 1, 1);
      _hVis_addBBMass        = bookHisto1D(44, 1, 1);
      _hVis_addBBDR_abs      = bookHisto1D(45, 1, 1);
      _hVis_addBBDR          = bookHisto1D(46, 1, 1);

      _hFull_addBJet1Pt_abs  = bookHisto1D(47, 1, 1);
      _hFull_addBJet1Pt      = bookHisto1D(48, 1, 1);
      _hFull_addBJet1Eta_abs = bookHisto1D(49, 1, 1);
      _hFull_addBJet1Eta     = bookHisto1D(50, 1, 1);
      _hFull_addBJet2Pt_abs  = bookHisto1D(51, 1, 1);
      _hFull_addBJet2Pt      = bookHisto1D(52, 1, 1);
      _hFull_addBJet2Eta_abs = bookHisto1D(53, 1, 1);
      _hFull_addBJet2Eta     = bookHisto1D(54, 1, 1);
      _hFull_addBBMass_abs   = bookHisto1D(55, 1, 1);
      _hFull_addBBMass       = bookHisto1D(56, 1, 1);
      _hFull_addBBDR_abs     = bookHisto1D(57, 1, 1);
      _hFull_addBBDR         = bookHisto1D(58, 1, 1);

      _h_gap_addJet1Pt       = bookProfile1D(59, 1, 1);
      _h_gap_addJet1Pt_eta0  = bookProfile1D(60, 1, 1);
      _h_gap_addJet1Pt_eta1  = bookProfile1D(61, 1, 1);
      _h_gap_addJet1Pt_eta2  = bookProfile1D(62, 1, 1);
      _h_gap_addJet2Pt       = bookProfile1D(63, 1, 1);
      _h_gap_addJet2Pt_eta0  = bookProfile1D(64, 1, 1);
      _h_gap_addJet2Pt_eta1  = bookProfile1D(65, 1, 1);
      _h_gap_addJet2Pt_eta2  = bookProfile1D(66, 1, 1);
      _h_gap_addJetHT        = bookProfile1D(67, 1, 1);
      _h_gap_addJetHT_eta0   = bookProfile1D(68, 1, 1);
      _h_gap_addJetHT_eta1   = bookProfile1D(69, 1, 1);
      _h_gap_addJetHT_eta2   = bookProfile1D(70, 1, 1);
    }


    void analyze(const Event& event) {

      // The objects used in the PAPER 12-041 are defined as follows (see p.16 for details):
      //
      //   * Leptons    : from the W boson decays after FSR
      //   * Jets       : anti-kT R=0.5 to all stable particles
      //                               exclude W->enu, munu, taunu
      //   * B jet      : B-Ghost matched
      //   * B from top : B hadron from top->b decay
      //
      // Visible phase space definition:
      //
      //   * Leptons         : pT > 20, |eta| < 2.4
      //   * B jets from top : pT > 30, |eta| < 2.4
      //     Additional jets : pT > 20, |eta| < 2.4
      //   *
      // Full phase space definition:
      //
      //   * Correction to dilepton BR from W boson BR
      //   * No cut on top decay products
      //   * Additional jets : pT > 20, |eta| < 2.4

      // Do the analysis only for the ttbar full leptonic channel, removing tau decays
      const Particles partontops = apply<ParticleFinder>(event, "PartonTops").particlesByPt();
      if (partontops.size() != 2) vetoEvent;
      const Particle& t1 = partontops[0];
      const Particle& t2 = partontops[1];

      // Apply acceptance cuts on top-decay leptons (existence should be guaranteed)
      const auto isPromptChLepton = [](const Particle& p){return isChargedLepton(p) && !fromDecay(p);};
      const Particle lep1 = t1.allDescendants(lastParticleWith(isPromptChLepton)).front();
      const Particle lep2 = t2.allDescendants(lastParticleWith(isPromptChLepton)).front();
      if (lep1.pT() < 1e-9*GeV || lep2.pT() < 1e-9*GeV) vetoEvent; // sanity check?

      const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.4);
      int nJet30 = 0, nJet60 = 0, nJet100 = 0;
      Jets topBJets, addJets, addBJets, addJets_eta0, addJets_eta1, addJets_eta2;
      for (const Jet& jet : jets) {
        if (jet.pT() >  30*GeV) nJet30 += 1;
        if (jet.pT() >  60*GeV) nJet60 += 1;
        if (jet.pT() > 100*GeV) nJet100 += 1;

        const bool isBtagged = jet.bTagged();
        const bool isBFromTop = any(jet.bTags(), hasParticleAncestorWith(Cuts::abspid == PID::TQUARK));

        if (isBFromTop) {
          if (jet.pT() > 30*GeV) topBJets.push_back(jet);
        } else {
          addJets.push_back(jet);
          if (isBtagged) addBJets.push_back(jet);
          if      (jet.abseta() < 0.8 ) addJets_eta0.push_back(jet);
          else if (jet.abseta() < 1.5 ) addJets_eta1.push_back(jet);
          else if (jet.abseta() < 2.4 ) addJets_eta2.push_back(jet);
        }
      }


      const bool isVisiblePS = topBJets.size() >= 2
        && lep1.pT() > 20*GeV && lep1.abseta() < 2.4 && lep2.pT() > 20*GeV && lep2.abseta() < 2.4;
      MSG_DEBUG(isVisiblePS << ": #b(top) = " << topBJets.size()
                << "; l1 = " << lep1.pT() << ", " << lep1.abseta()
                << "; l2 = " << lep2.pT() << ", " << lep2.abseta());

      const double weight = event.weight();


      if (isVisiblePS) {
        fillWithOF(_hVis_nJet30_abs,  nJet30, weight);
        fillWithOF(_hVis_nJet30,      nJet30, weight);
        fillWithOF(_hVis_nJet60_abs,  nJet60, weight);
        fillWithOF(_hVis_nJet60,      nJet60, weight);
        fillWithOF(_hVis_nJet100_abs, nJet100, weight);
        fillWithOF(_hVis_nJet100,     nJet100, weight);

        fillGapFractions(addJets, _h_gap_addJet1Pt, _h_gap_addJet2Pt, _h_gap_addJetHT, weight);
        fillGapFractions(addJets_eta0, _h_gap_addJet1Pt_eta0, _h_gap_addJet2Pt_eta0, _h_gap_addJetHT_eta0, weight);
        fillGapFractions(addJets_eta1, _h_gap_addJet1Pt_eta1, _h_gap_addJet2Pt_eta1, _h_gap_addJetHT_eta1, weight);
        fillGapFractions(addJets_eta2, _h_gap_addJet1Pt_eta2, _h_gap_addJet2Pt_eta2, _h_gap_addJetHT_eta2, weight);
      }

      // Plots with two additional jets
      if (addJets.size() >= 1) {
        const double ht = sum(addJets, pT, 0.0);
        _hFull_addJJHT_abs->fill(ht/GeV, weight);
        _hFull_addJJHT    ->fill(ht/GeV, weight);
        if (isVisiblePS) {
          _hVis_addJJHT_abs->fill(ht/GeV, weight);
          _hVis_addJJHT    ->fill(ht/GeV, weight);
        }

        const Jet& j1 = addJets[0];
        _hFull_addJet1Pt_abs ->fill(j1.pT()/GeV, weight);
        _hFull_addJet1Pt     ->fill(j1.pT()/GeV, weight);
        _hFull_addJet1Eta_abs->fill(j1.abseta(), weight);
        _hFull_addJet1Eta    ->fill(j1.abseta(), weight);
        if (isVisiblePS) {
          _hVis_addJet1Pt_abs ->fill(j1.pT()/GeV, weight);
          _hVis_addJet1Pt     ->fill(j1.pT()/GeV, weight);
          _hVis_addJet1Eta_abs->fill(j1.abseta(), weight);
          _hVis_addJet1Eta    ->fill(j1.abseta(), weight);
        }

        if (addJets.size() >= 2) {
          const Jet& j2 = addJets[1];

          _hFull_addJet2Pt_abs ->fill(j2.pT()/GeV, weight);
          _hFull_addJet2Pt     ->fill(j2.pT()/GeV, weight);
          _hFull_addJet2Eta_abs->fill(j2.abseta(), weight);
          _hFull_addJet2Eta    ->fill(j2.abseta(), weight);
          if (isVisiblePS) {
            _hVis_addJet2Pt_abs ->fill(j2.pT()/GeV, weight);
            _hVis_addJet2Pt     ->fill(j2.pT()/GeV, weight);
            _hVis_addJet2Eta_abs->fill(j2.abseta(), weight);
            _hVis_addJet2Eta    ->fill(j2.abseta(), weight);
          }

          const double jjmass = (j1.mom() + j2.mom()).mass();
          const double jjdR = deltaR(j1, j2);
          _hFull_addJJMass_abs->fill(jjmass/GeV, weight);
          _hFull_addJJMass    ->fill(jjmass/GeV, weight);
          _hFull_addJJDR_abs  ->fill(jjdR, weight);
          _hFull_addJJDR      ->fill(jjdR, weight);
          if (isVisiblePS) {
            _hVis_addJJMass_abs->fill(jjmass/GeV, weight);
            _hVis_addJJMass    ->fill(jjmass/GeV, weight);
            _hVis_addJJDR_abs  ->fill(jjdR, weight);
            _hVis_addJJDR      ->fill(jjdR, weight);
          }
        }
      }


      // Same set of plots if there are additional b-jets
      if (addBJets.size() >= 1) {
        const Jet& b1 = addBJets[0];
        _hFull_addBJet1Pt_abs ->fill(b1.pT()/GeV, weight);
        _hFull_addBJet1Pt     ->fill(b1.pT()/GeV, weight);
        _hFull_addBJet1Eta_abs->fill(b1.abseta(), weight);
        _hFull_addBJet1Eta    ->fill(b1.abseta(), weight);
        if (isVisiblePS) {
          _hVis_addBJet1Pt_abs ->fill(b1.pT()/GeV, weight);
          _hVis_addBJet1Pt     ->fill(b1.pT()/GeV, weight);
          _hVis_addBJet1Eta_abs->fill(b1.abseta(), weight);
          _hVis_addBJet1Eta    ->fill(b1.abseta(), weight);
        }

        if (addBJets.size() >= 2) {
          const Jet& b2 = addBJets[1];

          _hFull_addBJet2Pt_abs ->fill(b2.pT()/GeV, weight);
          _hFull_addBJet2Pt     ->fill(b2.pT()/GeV, weight);
          _hFull_addBJet2Eta_abs->fill(b2.abseta(), weight);
          _hFull_addBJet2Eta    ->fill(b2.abseta(), weight);
          if (isVisiblePS) {
            _hVis_addBJet2Pt_abs ->fill(b2.pT()/GeV, weight);
            _hVis_addBJet2Pt     ->fill(b2.pT()/GeV, weight);
            _hVis_addBJet2Eta_abs->fill(b2.abseta(), weight);
            _hVis_addBJet2Eta    ->fill(b2.abseta(), weight);
          }

          const double bbmass = (b1.mom() + b2.mom()).mass();
          const double bbdR = deltaR(b1, b2);
          _hFull_addBBMass_abs->fill(bbmass/GeV, weight);
          _hFull_addBBMass    ->fill(bbmass/GeV, weight);
          _hFull_addBBDR_abs  ->fill(bbdR, weight);
          _hFull_addBBDR      ->fill(bbdR, weight);
          if (isVisiblePS) {
            _hVis_addBBMass_abs->fill(bbmass/GeV, weight);
            _hVis_addBBMass    ->fill(bbmass/GeV, weight);
            _hVis_addBBDR_abs  ->fill(bbdR, weight);
            _hVis_addBBDR      ->fill(bbdR, weight);
          }
        }
      }

    }


    void finalize() {
      const double ttbarXS = !std::isnan(crossSectionPerEvent()) ? crossSection() : 252.89*picobarn;
      if (std::isnan(crossSectionPerEvent()))
        MSG_INFO("No valid cross-section given, using NNLO (arXiv:1303.6254; sqrt(s)=8 TeV, m_t=172.5 GeV): " << ttbarXS/picobarn << " pb");

      normalize({_hVis_nJet30,_hVis_nJet60, _hVis_nJet100,
            _hVis_addJet1Pt, _hVis_addJet1Eta, _hVis_addJet2Pt, _hVis_addJet2Eta,
            _hVis_addJJMass, _hVis_addJJDR, _hVis_addJJHT,
            _hFull_addJet1Pt, _hFull_addJet1Eta, _hFull_addJet2Pt, _hFull_addJet2Eta,
            _hFull_addJJMass, _hFull_addJJDR, _hFull_addJJHT,
            _hVis_addBJet1Pt, _hVis_addBJet1Eta, _hVis_addBJet2Pt, _hVis_addBJet2Eta,
            _hVis_addBBMass, _hVis_addBBDR,
            _hFull_addBJet1Pt, _hFull_addBJet1Eta, _hFull_addBJet2Pt, _hFull_addBJet2Eta,
            _hFull_addBBMass, _hFull_addBBDR});

      const double xsPerWeight = ttbarXS/picobarn / sumOfWeights();
      scale({_hVis_nJet30_abs, _hVis_nJet60_abs, _hVis_nJet100_abs,
            _hVis_addJet1Pt_abs, _hVis_addJet1Eta_abs, _hVis_addJet2Pt_abs, _hVis_addJet2Eta_abs,
            _hVis_addJJMass_abs, _hVis_addJJDR_abs, _hVis_addJJHT_abs,
            _hVis_addBJet1Pt_abs, _hVis_addBJet1Eta_abs, _hVis_addBJet2Pt_abs, _hVis_addBJet2Eta_abs,
            _hVis_addBBMass_abs, _hVis_addBBDR_abs}, xsPerWeight);

      const double sfull = xsPerWeight / 0.0454; //< correct for dilepton branching fraction
      scale({_hFull_addJet1Pt_abs, _hFull_addJet1Eta_abs, _hFull_addJet2Pt_abs, _hFull_addJet2Eta_abs,
            _hFull_addJJMass_abs, _hFull_addJJDR_abs, _hFull_addJJHT_abs,
            _hFull_addBJet1Pt_abs, _hFull_addBJet1Eta_abs, _hFull_addBJet2Pt_abs, _hFull_addBJet2Eta_abs,
            _hFull_addBBMass_abs, _hFull_addBBDR_abs}, sfull);
    }

    //@}


    void fillWithOF(Histo1DPtr h, double x, double w) {
      h->fill(std::min(x, h->xMax()-1e-9), w);
    }


    void fillGapFractions(const Jets& addJets, Profile1DPtr h_gap_addJet1Pt, Profile1DPtr h_gap_addJet2Pt, Profile1DPtr h_gap_addJetHT, double weight) {
      const double j1pt = (addJets.size() > 0) ? addJets[0].pT() : 0;
      for (size_t i = 0; i < h_gap_addJet1Pt->numBins(); ++i) {
        const double binCenter = h_gap_addJet1Pt->bin(i).xMid();
        h_gap_addJet1Pt->fillBin(i, int(j1pt/GeV < binCenter), weight);
      }

      const double j2pt = (addJets.size() > 1) ? addJets[1].pT() : 0;
      for (size_t i = 0; i < h_gap_addJet2Pt->numBins(); ++i) {
        const double binCenter = h_gap_addJet2Pt->bin(i).xMid();
        h_gap_addJet2Pt->fillBin(i, int(j2pt/GeV < binCenter), weight);
      }

      const double ht = sum(addJets, pT, 0.);
      for (size_t i = 0; i < h_gap_addJetHT->numBins(); ++i) {
        const double binCenter = h_gap_addJetHT->bin(i).xMid();
        h_gap_addJetHT->fillBin(i, int(ht/GeV < binCenter) , weight);
      }
    }


    // @name Histogram data members
    //@{

    Histo1DPtr _hVis_nJet30_abs, _hVis_nJet60_abs, _hVis_nJet100_abs;
    Histo1DPtr _hVis_addJet1Pt_abs, _hVis_addJet1Eta_abs, _hVis_addJet2Pt_abs, _hVis_addJet2Eta_abs;
    Histo1DPtr _hVis_addJJMass_abs, _hVis_addJJDR_abs, _hVis_addJJHT_abs;
    Histo1DPtr _hFull_addJet1Pt_abs, _hFull_addJet1Eta_abs, _hFull_addJet2Pt_abs, _hFull_addJet2Eta_abs;
    Histo1DPtr _hFull_addJJMass_abs, _hFull_addJJDR_abs, _hFull_addJJHT_abs;
    Histo1DPtr _hVis_addBJet1Pt_abs, _hVis_addBJet1Eta_abs, _hVis_addBJet2Pt_abs, _hVis_addBJet2Eta_abs;
    Histo1DPtr _hVis_addBBMass_abs, _hVis_addBBDR_abs;
    Histo1DPtr _hFull_addBJet1Pt_abs, _hFull_addBJet1Eta_abs, _hFull_addBJet2Pt_abs, _hFull_addBJet2Eta_abs;
    Histo1DPtr _hFull_addBBMass_abs, _hFull_addBBDR_abs;

    Histo1DPtr _hVis_nJet30, _hVis_nJet60, _hVis_nJet100;
    Histo1DPtr _hVis_addJet1Pt, _hVis_addJet1Eta, _hVis_addJet2Pt, _hVis_addJet2Eta;
    Histo1DPtr _hVis_addJJMass, _hVis_addJJDR, _hVis_addJJHT;
    Histo1DPtr _hFull_addJet1Pt, _hFull_addJet1Eta, _hFull_addJet2Pt, _hFull_addJet2Eta;
    Histo1DPtr _hFull_addJJMass, _hFull_addJJDR, _hFull_addJJHT;
    Histo1DPtr _hVis_addBJet1Pt, _hVis_addBJet1Eta, _hVis_addBJet2Pt, _hVis_addBJet2Eta;
    Histo1DPtr _hVis_addBBMass, _hVis_addBBDR;
    Histo1DPtr _hFull_addBJet1Pt, _hFull_addBJet1Eta, _hFull_addBJet2Pt, _hFull_addBJet2Eta;
    Histo1DPtr _hFull_addBBMass, _hFull_addBBDR;

    Profile1DPtr _h_gap_addJet1Pt, _h_gap_addJet1Pt_eta0, _h_gap_addJet1Pt_eta1, _h_gap_addJet1Pt_eta2;
    Profile1DPtr _h_gap_addJet2Pt, _h_gap_addJet2Pt_eta0, _h_gap_addJet2Pt_eta1, _h_gap_addJet2Pt_eta2;
    Profile1DPtr _h_gap_addJetHT, _h_gap_addJetHT_eta0, _h_gap_addJetHT_eta1, _h_gap_addJetHT_eta2;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1397174);


}
#line 1 "CMS_2016_I1473674.cc"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

  class CMS_2016_I1473674 : public Analysis {
  public:

    // Minimal constructor
    CMS_2016_I1473674() : Analysis("CMS_2016_I1473674") {
    }

    // Set up projections and book histograms
    void init() {
      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);

      // Projection for taus
      TauFinder taus(TauFinder::ANY);
      addProjection(taus, "Tau");
      IdentifiedFinalState nu_taus(fs);
      nu_taus.acceptIdPair(PID::NU_TAU);
      addProjection(nu_taus, "NuTau");

      // Projection for electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      addProjection(electrons, "Electrons");
      DressedLeptons dressed_electrons(photons, electrons, 0.1, Cuts::open(), true, false);
      addProjection(dressed_electrons, "DressedElectrons");

      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      addProjection(muons, "Muons");
      DressedLeptons dressed_muons(photons, muons, 0.1, Cuts::open(), true, false);
      addProjection(dressed_muons, "DressedMuons");

      // Projection for jets
      VetoedFinalState fs_jets(FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV));
      fs_jets.addVetoOnThisFinalState(dressed_muons);
      addProjection(FastJets(fs_jets, FastJets::ANTIKT, 0.5), "Jets");

      // Projections for MET
      addProjection(MissingMomentum(), "MET");

      // Booking of histograms
      _hist_met = bookHisto1D(5, 1, 1);
      _hist_ht  = bookHisto1D(6, 1, 1);
      _hist_st  = bookHisto1D(7, 1, 1);
      _hist_wpt = bookHisto1D(8, 1, 1);
    }


    // per event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // select ttbar -> lepton+jets
      const DressedLeptons& dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");
      if (dressed_electrons.dressedLeptons().size() +
          dressed_muons.dressedLeptons().size() != 1) {
        vetoEvent;
      }

      FourMomentum lepton;
      if (dressed_electrons.dressedLeptons().size() == 1) {
        lepton = dressed_electrons.dressedLeptons()[0].momentum();
      } else {
        lepton = dressed_muons.dressedLeptons()[0].momentum();
      }

      // veto if lepton is tau
      const TauFinder& taus = applyProjection<TauFinder>(event, "Tau");
      const IdentifiedFinalState nu_taus = applyProjection<IdentifiedFinalState>(event, "NuTau");
      foreach (const Particle& tau, taus.taus()) {
        foreach (const Particle& nu, nu_taus.particles()) {
          if (tau.pid() * nu.pid() < 0)
            continue;

          const FourMomentum w_candidate = tau.momentum() + nu.momentum();
          if (abs(w_candidate.mass() - 80.4) > 5.)
            vetoEvent;
        }
      }


      // MET
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      _hist_met->fill(met.visibleMomentum().pT()/GeV, weight);

      // HT and ST
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets jets = jetpro.jetsByPt(20*GeV);

      double ht = 0.0;
      foreach (const Jet& j, jets) {
        if (deltaR(j.momentum(), lepton) > 0.3) {
          ht += j.pT();
        }
      }

      double st = ht + lepton.pT() + met.visibleMomentum().pT();
      _hist_ht->fill(ht/GeV, weight);
      _hist_st->fill(st/GeV, weight);

      // WPT
      FourMomentum w = lepton - met.visibleMomentum();
      _hist_wpt->fill(w.pT()/GeV, weight);
    }

    // scale by 1 over weight
    void finalize() {
      normalize(_hist_met);
      normalize(_hist_ht);
      normalize(_hist_st);
      normalize(_hist_wpt);
    }

  private:
    Histo1DPtr _hist_met, _hist_ht, _hist_st, _hist_wpt;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1473674);
}
#line 1 "CMSTOTEM_2014_I1294140.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  class CMSTOTEM_2014_I1294140 : public Analysis {
  public:

    CMSTOTEM_2014_I1294140()
      : Analysis("CMSTOTEM_2014_I1294140")
    {     }


    void init() {
      ChargedFinalState cfs(-7.0, 7.0, 0.0*GeV);
      declare(cfs, "CFS");

      _Nevt_after_cuts_or = 0;
      _Nevt_after_cuts_and = 0;
      _Nevt_after_cuts_xor = 0;

      if (fuzzyEquals(sqrtS(), 8000*GeV, 1E-3)) {
        _h_dNch_dEta_OR = bookHisto1D(1, 1, 1);
        _h_dNch_dEta_AND = bookHisto1D(2, 1, 1);
        _h_dNch_dEta_XOR = bookHisto1D(3, 1, 1);
      }
    }


    void analyze(const Event& event) {
      // Count forward and backward charged particles
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");
      int count_plus = 0, count_minus = 0;
      foreach (const Particle& p, charged.particles()) {
        if (inRange(p.eta(),  5.3,  6.5)) count_plus++;
        if (inRange(p.eta(), -6.5, -5.3)) count_minus++;
      }

      // Cut combinations
      const bool cutsor  = (count_plus > 0 || count_minus > 0);
      const bool cutsand = (count_plus > 0 && count_minus > 0);
      const bool cutsxor = ( (count_plus > 0 && count_minus == 0) || (count_plus == 0 && count_minus > 0) );

      // Increment counters and fill histos
      const double weight = event.weight();
      if (cutsor)  _Nevt_after_cuts_or  += weight;
      if (cutsand) _Nevt_after_cuts_and += weight;
      if (cutsxor) _Nevt_after_cuts_xor += weight;
      foreach (const Particle& p, charged.particles()) {
        if (cutsor)  _h_dNch_dEta_OR ->fill(p.abseta(), weight);
        if (cutsand) _h_dNch_dEta_AND->fill(p.abseta(), weight);
        if (cutsxor) _h_dNch_dEta_XOR->fill(p.abseta(), weight);
      }

    }


    void finalize() {
      scale(_h_dNch_dEta_OR,  0.5/_Nevt_after_cuts_or);
      scale(_h_dNch_dEta_AND, 0.5/_Nevt_after_cuts_and);
      scale(_h_dNch_dEta_XOR, 0.5/_Nevt_after_cuts_xor);
    }


  private:

    Histo1DPtr _h_dNch_dEta_OR, _h_dNch_dEta_AND, _h_dNch_dEta_XOR;
    double _Nevt_after_cuts_or, _Nevt_after_cuts_and, _Nevt_after_cuts_xor;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMSTOTEM_2014_I1294140);

}
#line 1 "TOTEM_2014_I1328627.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class TOTEM_2014_I1328627 : public Analysis {
  public:


    TOTEM_2014_I1328627()
      : Analysis("TOTEM_2014_I1328627")
    {    }



    void init() {
      ChargedFinalState cfsm(-7.0, -6.0, 0.0*GeV);
      ChargedFinalState cfsp( 3.7,  4.8, 0.0*GeV);
      declare(cfsm, "CFSM");
      declare(cfsp, "CFSP");

      _h_eta = bookHisto1D(1, 1, 1);
      _sumofweights = 0.;
    }


    void analyze(const Event& event) {
      const ChargedFinalState cfsm = apply<ChargedFinalState>(event, "CFSM");
      const ChargedFinalState cfsp = apply<ChargedFinalState>(event, "CFSP");
      if (cfsm.size() == 0 && cfsp.size() == 0) vetoEvent;

      _sumofweights += event.weight();
      foreach (const Particle& p, cfsm.particles() + cfsp.particles()) {
        _h_eta->fill(p.abseta(), event.weight());
      }
    }


    void finalize() {
      scale(_h_eta, 1./_sumofweights);
    }


  private:

    double _sumofweights;
    Histo1DPtr _h_eta;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TOTEM_2014_I1328627);


}
#line 1 "CMS_QCD_10_024.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Particle.hh"

namespace Rivet {


  class CMS_QCD_10_024 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_QCD_10_024() : Analysis("CMS_QCD_10_024"),
		       _weight_pt05_eta08(0.), _weight_pt10_eta08(0.),
		       _weight_pt05_eta24(0.), _weight_pt10_eta24(0.) {  }


    void init() {
      declare(ChargedFinalState(-0.8, 0.8, 0.5*GeV), "CFS_08_05");
      declare(ChargedFinalState(-0.8, 0.8, 1.0*GeV), "CFS_08_10");
      declare(ChargedFinalState(-2.4, 2.4, 0.5*GeV), "CFS_24_05");
      declare(ChargedFinalState(-2.4, 2.4, 1.0*GeV), "CFS_24_10");

      size_t offset = 0;
      if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) offset = 0;
      if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3)) offset = 4;
      _hist_dNch_deta_pt05_eta08 = bookHisto1D(1+offset, 1, 1);
      _hist_dNch_deta_pt10_eta08 = bookHisto1D(2+offset, 1, 1);
      _hist_dNch_deta_pt05_eta24 = bookHisto1D(3+offset, 1, 1);
      _hist_dNch_deta_pt10_eta24 = bookHisto1D(4+offset, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& cfs_08_05 = apply<ChargedFinalState>(event, "CFS_08_05");
      const ChargedFinalState& cfs_08_10 = apply<ChargedFinalState>(event, "CFS_08_10");
      const ChargedFinalState& cfs_24_05 = apply<ChargedFinalState>(event, "CFS_24_05");
      const ChargedFinalState& cfs_24_10 = apply<ChargedFinalState>(event, "CFS_24_10");

      // Plot distributions
      if(!cfs_08_05.particles().empty()) _weight_pt05_eta08 += weight;
      if(!cfs_24_05.particles().empty()) _weight_pt05_eta24 += weight;
      foreach (const Particle& p, cfs_24_05.particles()) {
        _hist_dNch_deta_pt05_eta24->fill(p.eta(), weight);
        if(!cfs_08_05.particles().empty())
	  _hist_dNch_deta_pt05_eta08->fill(p.eta(), weight);
      }
      if(!cfs_08_10.particles().empty()) _weight_pt10_eta08 += weight;
      if(!cfs_24_10.particles().empty()) _weight_pt10_eta24 += weight;
      foreach (const Particle& p, cfs_24_10.particles()) {
        _hist_dNch_deta_pt10_eta24->fill(p.eta(), weight);
	if(!cfs_08_10.particles().empty())
	  _hist_dNch_deta_pt10_eta08->fill(p.eta(), weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_dNch_deta_pt05_eta08,1./_weight_pt05_eta08);
      scale(_hist_dNch_deta_pt10_eta08,1./_weight_pt10_eta08);
      scale(_hist_dNch_deta_pt05_eta24,1./_weight_pt05_eta24);
      scale(_hist_dNch_deta_pt10_eta24,1./_weight_pt10_eta24);
    }


  private:

    Histo1DPtr _hist_dNch_deta_pt05_eta08;
    Histo1DPtr _hist_dNch_deta_pt10_eta08;
    Histo1DPtr _hist_dNch_deta_pt05_eta24;
    Histo1DPtr _hist_dNch_deta_pt10_eta24;
    double _weight_pt05_eta08,_weight_pt10_eta08,_weight_pt05_eta24,_weight_pt10_eta24;
  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_QCD_10_024);

}
#line 1 "CMS_2012_PAS_QCD_11_010.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  class CMS_2012_PAS_QCD_11_010 : public Analysis {
  public:

    CMS_2012_PAS_QCD_11_010()
      : Analysis("CMS_2012_PAS_QCD_11_010")
    {  }

    void init() {
      const FastJets jets(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 0.5*GeV), FastJets::ANTIKT, 0.5);
      declare(jets, "Jets");

      const UnstableFinalState ufs(Cuts::abseta < 2 && Cuts::pT > 0.6*GeV);
      declare(ufs, "UFS");

      _h_nTrans_Lambda     = bookProfile1D(1, 1, 1);
      _h_nTrans_Kaon       = bookProfile1D(2, 1, 1);
      _h_ptsumTrans_Lambda = bookProfile1D(3, 1, 1);
      _h_ptsumTrans_Kaon   = bookProfile1D(4, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(1.0*GeV);
      if (jets.size() < 1) vetoEvent;

      if (fabs(jets[0].eta()) >= 2) { // cuts on leading jets
        vetoEvent;
      }

      FourMomentum p_lead = jets[0].momentum();
      const double pTlead  = p_lead.pT();

      const UnstableFinalState& ufs = apply<UnstableFinalState>(event, "UFS");

      int numTrans_Kaon = 0;
      int numTrans_Lambda = 0;
      double ptSumTrans_Kaon = 0.;
      double ptSumTrans_Lambda = 0.;

      foreach (const Particle& p, ufs.particles()) {
        double dphi = deltaPhi(p, p_lead);
        double pT = p.pT();
        const PdgId id = p.abspid();

        if (dphi > PI/3. && dphi < 2./3.*PI) {
          if (id == 310 && pT > 0.6*GeV) {
            ptSumTrans_Kaon += pT/GeV;
            numTrans_Kaon++;
          }
          else if (id == 3122 && pT > 1.5*GeV) {
            ptSumTrans_Lambda += pT/GeV;
            numTrans_Lambda++;
          }
        }
      }

      _h_nTrans_Kaon->fill(pTlead/GeV, numTrans_Kaon / (8.0 * PI/3.0), weight);
      _h_nTrans_Lambda->fill(pTlead/GeV, numTrans_Lambda / (8.0 * PI/3.0), weight);
      _h_ptsumTrans_Kaon->fill(pTlead/GeV, ptSumTrans_Kaon / (GeV * (8.0 * PI/3.0)), weight);
      _h_ptsumTrans_Lambda->fill(pTlead/GeV, ptSumTrans_Lambda / (GeV * (8.0 * PI/3.0)), weight);
    }


    void finalize() { }

  private:

    Profile1DPtr _h_nTrans_Kaon;
    Profile1DPtr _h_nTrans_Lambda;
    Profile1DPtr _h_ptsumTrans_Kaon;
    Profile1DPtr _h_ptsumTrans_Lambda;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_PAS_QCD_11_010);

}
#line 1 "CMS_2016_PAS_SUS_16_14.cc"
// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedMET.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// @brief CMS 2016 0-lepton SUSY search, from 13/fb PAS note
  class CMS_2016_PAS_SUS_16_14 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_PAS_SUS_16_14);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState calofs(Cuts::abseta < 5.0);
      FastJets fj(calofs, FastJets::ANTIKT, 0.4);
      declare(fj, "TruthJets");
      declare(SmearedJets(fj, JET_SMEAR_CMS_RUN2, [](const Jet& j) {
            if (j.abseta() > 2.5) return 0.;
            return j.bTagged() ? 0.55 : j.cTagged() ? 0.12 : 0.016; }), "Jets");

      FinalState es(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.5);
      declare(es, "TruthElectrons");
      declare(SmearedParticles(es, ELECTRON_EFF_CMS_RUN2, ELECTRON_SMEAR_CMS_RUN2), "Electrons");

      FinalState mus(Cuts::abspid == PID::MUON && Cuts::abseta < 2.4);
      declare(mus, "TruthMuons");
      declare(SmearedParticles(mus, MUON_EFF_CMS_RUN2, MUON_SMEAR_CMS_RUN2), "Muons");

      FinalState isofs(Cuts::abseta < 3.0 && Cuts::abspid != PID::ELECTRON && Cuts::abspid != PID::MUON);
      declare(isofs, "IsoFS");
      FinalState cfs(Cuts::abseta < 2.5 && Cuts::abscharge != 0);
      declare(cfs, "TruthTracks");
      declare(SmearedParticles(cfs, TRK_EFF_CMS_RUN2), "Tracks");

      // Book histograms/counters
      _h_srcounts.resize(160);
      for (size_t ij = 0; ij < 4; ++ij) {
        for (size_t ib = 0; ib < 4; ++ib) {
          for (size_t ih = 0; ih < 10; ++ih) {
            const size_t i = 40*ij + 10*ib + ih;
            _h_srcounts[i] = bookCounter(toString(2*ij+3) + "j-" + toString(ib) + "b-" + toString(ih));
          }
        }
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get jets and require Nj >= 3
      const Jets jets24 = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 2.4);
      if (jets24.size() < 3) vetoEvent;

      // HT cut
      vector<double> jetpts24; transform(jets24, jetpts24, pT);
      const double ht = sum(jetpts24, 0.0);
      if (ht < 300*GeV) vetoEvent;

      // HTmiss cut
      const Jets jets50 = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 5.0);
      const FourMomentum htmissvec = -sum(jets24, mom, FourMomentum());
      const double htmiss = htmissvec.pT();
      if (htmissvec.pT() < 300*GeV) vetoEvent;


      // Get baseline electrons & muons
      Particles elecs = apply<ParticleFinder>(event, "Electrons").particles(Cuts::pT > 10*GeV);
      Particles muons = apply<ParticleFinder>(event, "Muons").particles(Cuts::pT > 10*GeV);

      // Electron/muon isolation
      const Particles calofs = apply<ParticleFinder>(event, "IsoFS").particles();
      ifilter_discard(elecs, [&](const Particle& e) {
          const double R = max(0.05, min(0.2, 10*GeV/e.pT()));
          double ptsum = -e.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,e) < R) ptsum += p.pT();
          return ptsum / e.pT() > 0.1;
        });
      ifilter_discard(muons, [&](const Particle& m) {
          const double R = max(0.05, min(0.2, 10*GeV/m.pT()));
          double ptsum = -m.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,m) < R) ptsum += p.pT();
          return ptsum / m.pT() > 0.2;
        });

      // Veto the event if there are any remaining baseline leptons
      if (!elecs.empty()) vetoEvent;
      if (!muons.empty()) vetoEvent;


      // Get isolated tracks
      Particles trks25 = apply<ParticleFinder>(event, "Tracks").particles();
      ifilter_discard(trks25, [&](const Particle& t) {
          double ptsum = -t.pT();
          for (const Particle& p : trks25)
            if (deltaR(p,t) < 0.3) ptsum += p.pT();
          return ptsum/t.pT() > ((t.abspid() == PID::ELECTRON || t.abspid() == PID::MUON) ? 0.2 : 0.1);
        });
      const Particles trks = filter_select(trks25, Cuts::abseta < 2.4);

      // Isolated track pT, pTmiss and mT cut
      // mT^2 = m1^2 + m2^2 + 2(ET1 ET2 - pT1 . pT2))
      // => mT0^2 = 2(ET1 |pT2| - pT1 . pT2)) for m1, m2 -> 0
      FourMomentum ptmissvec = htmissvec; ///< @todo Can we do better? No e,mu left...
      const double ptmiss = ptmissvec.pT();
      for (const Particle& t : trks) {
        const double ptcut = (t.abspid() == PID::ELECTRON || t.abspid() == PID::MUON) ? 5*GeV : 10*GeV;
        const double mT = sqrt( t.mass2() + 2*(t.Et()*ptmiss - t.pT()*ptmiss*cos(deltaPhi(t,ptmissvec))) );
        if (mT < 100*GeV && t.pT() < ptcut) vetoEvent;
      }

      // Lead jets isolation from Htmiss
      if (deltaPhi(htmissvec, jets24[0]) < 0.5) vetoEvent;
      if (deltaPhi(htmissvec, jets24[1]) < 0.5) vetoEvent;
      if (deltaPhi(htmissvec, jets24[2]) < 0.3) vetoEvent;
      if (jets24.size() >= 4 && deltaPhi(htmissvec, jets24[3]) < 0.3) vetoEvent;


      ////////


      // Calculate a bin index for this event
      // Nj bin
      static const vector<double> njedges = {3., 5., 7., 9.};
      const size_t nj = jets24.size();
      // Nbj bin
      static const vector<double> njbedges = {0., 1., 2., 3.};
      const size_t inj = binIndex(nj, njedges, true);
      size_t nbj = 0;
      for (const Jet& j : jets24)
        if (j.bTagged()) nbj += 1;
      const size_t inbj = binIndex(nbj, njbedges, true);
      // HTmiss vs HT 2D bin
      int iht = 0;
      if (htmiss < 350*GeV) {
        iht = ht < 500 ? 1 : ht < 1000 ? 2 : 3;
      } if (htmiss < 500*GeV && ht > 350*GeV) {
        iht = ht < 500 ? 4 : ht < 1000 ? 5 : 6;
      } if (htmiss < 750*GeV && ht > 500*GeV) {
        iht = ht < 1000 ? 7 : 8;
      } if (ht > 750*GeV) {
        iht = ht < 1500 ? 9 : 10;
      }
      if (iht == 0) vetoEvent;
      iht -= 1; //< change from the paper's indexing scheme to C++ zero-indexed
      // Total bin number
      const size_t ibin = 40*inj + 10*inbj + (size_t)iht;

      // Fill SR counter
      _h_srcounts[ibin]->fill(event.weight());

    }


    /// Normalise counters after the run
    void finalize() {

      const double sf = 12.9*crossSection()/femtobarn/sumOfWeights();
      scale(_h_srcounts, sf);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    vector<CounterPtr> _h_srcounts;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_PAS_SUS_16_14);


}
