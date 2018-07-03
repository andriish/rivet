#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  class ATLAS_2014_I1319490 : public Analysis {
  public:

    ATLAS_2014_I1319490(string name = "ATLAS_2014_I1319490")
      : Analysis(name)
    {
      _mode = 0; // using electron channel for combined data by default

    }


    // Book histograms and initialise projections before the run
    void init() {

      FinalState fs;

      Cut cuts;
      if (_mode == 2) { // muon channel
        cuts = (Cuts::pT > 25.0*GeV) & Cuts::etaIn(-2.4, 2.4);
      } else if (_mode) { // electron channel
        cuts = (Cuts::pT > 25.0*GeV) & ( Cuts::etaIn(-2.47, -1.52) | Cuts::etaIn(-1.37, 1.37) | Cuts::etaIn(1.52, 2.47) );
      } else { // combined data extrapolated to common phase space
        cuts = (Cuts::pT > 25.0*GeV) & Cuts::etaIn(-2.5, 2.5);
      }

      // bosons
      WFinder wfinder(fs, cuts, _mode > 1? PID::MUON : PID::ELECTRON, 40.0*GeV, DBL_MAX, 0.0*GeV, 0.1,
                      WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
      declare(wfinder, "WF");

      // jets
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(getProjection<WFinder>("WF"));
      FastJets jets(jet_fs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      declare(jets, "Jets");

      // book histograms
      book(histos["h_N_incl"]            ,1,1,_mode+1);
      book(histos["h_N"]                 ,4,1,_mode+1);
      book(histos["h_pt_jet1_1jet"]      ,5,1,_mode+1);
      book(histos["h_pt_jet1_1jet_excl"] ,6,1,_mode+1);
      book(histos["h_pt_jet1_2jet"]      ,7,1,_mode+1);
      book(histos["h_pt_jet1_3jet"]      ,8,1,_mode+1);
      book(histos["h_pt_jet2_2jet"]      ,9,1,_mode+1);
      book(histos["h_pt_jet3_3jet"]      ,10,1,_mode+1);
      book(histos["h_pt_jet4_4jet"]      ,11,1,_mode+1);
      book(histos["h_pt_jet5_5jet"]      ,12,1,_mode+1);
      book(histos["h_y_jet1_1jet"]       ,13,1,_mode+1);
      book(histos["h_y_jet2_2jet"]       ,14,1,_mode+1);
      book(histos["h_HT_1jet"]           ,15,1,_mode+1);
      book(histos["h_HT_1jet_excl"]      ,16,1,_mode+1);
      book(histos["h_HT_2jet"]           ,17,1,_mode+1);
      book(histos["h_HT_2jet_excl"]      ,18,1,_mode+1);
      book(histos["h_HT_3jet"]           ,19,1,_mode+1);
      book(histos["h_HT_3jet_excl"]      ,20,1,_mode+1);
      book(histos["h_HT_4jet"]           ,21,1,_mode+1);
      book(histos["h_HT_5jet"]           ,22,1,_mode+1);
      book(histos["h_deltaPhi_jet12"]    ,23,1,_mode+1);
      book(histos["h_deltaRap_jet12"]    ,24,1,_mode+1);
      book(histos["h_deltaR_jet12"]      ,25,1,_mode+1);
      book(histos["h_M_Jet12_2jet"]      ,26,1,_mode+1);
      book(histos["h_y_jet3_3jet"]       ,27,1,_mode+1);
      book(histos["h_y_jet4_4jet"]       ,28,1,_mode+1);
      book(histos["h_y_jet5_5jet"]       ,29,1,_mode+1);
      book(histos["h_ST_1jet"]           ,30,1,_mode+1);
      book(histos["h_ST_2jet"]           ,31,1,_mode+1);
      book(histos["h_ST_2jet_excl"]      ,32,1,_mode+1);
      book(histos["h_ST_3jet"]           ,33,1,_mode+1);
      book(histos["h_ST_3jet_excl"]      ,34,1,_mode+1);
      book(histos["h_ST_4jet"]           ,35,1,_mode+1);
      book(histos["h_ST_5jet"]           ,36,1,_mode+1);
    }


    void fillPlots(const Particle& lepton, const double& missET, Jets& all_jets) {
      // do jet-lepton overlap removal
      Jets jets;
      double ST = 0.0; // scalar pT sum of all selected jets
      for (const Jet &j : all_jets) {
        if (deltaR(j, lepton) > 0.5) {
          jets += j;
          ST += j.pT() / GeV;
        }
      }

      const size_t njets = jets.size();

      const double HT = ST + lepton.pT() / GeV + missET;

      histos["h_N"]->fill(njets + 0.5);
      for (size_t i = 0; i <= njets; ++i) {
        histos["h_N_incl"]->fill(i + 0.5);
      }

      if (njets) {
        const double pT1  = jets[0].pT() / GeV;
        const double rap1 = jets[0].absrap();
        histos["h_pt_jet1_1jet" ]->fill(pT1);
        histos["h_y_jet1_1jet"]->fill(rap1);
        histos["h_HT_1jet"]->fill(HT);
        histos["h_ST_1jet"]->fill(ST);
        if (njets == 1) {
          histos["h_pt_jet1_1jet_excl"]->fill(pT1);
          histos["h_HT_1jet_excl"]->fill(HT);
        } else {
          const double pT2  = jets[1].pT() / GeV;
          const double rap2 = jets[1].absrap();
          const double dR   = deltaR(jets[0], jets[1]);
          const double dRap = deltaRap(jets[0], jets[1]);
          const double dPhi = deltaPhi(jets[0], jets[1]);
          const double mjj  = (jets[0].momentum() + jets[1].momentum()).mass() / GeV;
          histos["h_pt_jet1_2jet"]->fill(pT1);
          histos["h_pt_jet2_2jet"]->fill(pT2);
          histos["h_y_jet2_2jet"]->fill(rap2);
          histos["h_M_Jet12_2jet"]->fill(mjj);
          histos["h_HT_2jet"]->fill(HT);
          histos["h_ST_2jet"]->fill(ST);
          histos["h_deltaPhi_jet12"]->fill(dPhi);
          histos["h_deltaRap_jet12"]->fill(dRap);
          histos["h_deltaR_jet12"]->fill(dR);
          if (njets == 2) {
            histos["h_ST_2jet_excl"]->fill(ST);
            histos["h_HT_2jet_excl"]->fill(HT);
          } else {
            const double pT3  = jets[2].pT() / GeV;
            const double rap3 = jets[2].absrap();
            histos["h_pt_jet1_3jet"]->fill(pT1);
            histos["h_pt_jet3_3jet"]->fill(pT3);
            histos["h_y_jet3_3jet"]->fill(rap3);
            histos["h_HT_3jet"]->fill(HT);
            histos["h_ST_3jet"]->fill(ST);
            if(njets == 3) {
              histos["h_ST_3jet_excl"]->fill(ST);
              histos["h_HT_3jet_excl"]->fill(HT);
            } else {
              const double pT4  = jets[3].pT() / GeV;
              const double rap4 = jets[3].absrap();
              histos["h_pt_jet4_4jet"]->fill(pT4);
              histos["h_y_jet4_4jet"]->fill(rap4);
              histos["h_HT_4jet"]->fill(HT);
              histos["h_ST_4jet"]->fill(ST);
              if (njets > 4) {
                const double pT5  = jets[4].pT() / GeV;
                const double rap5 = jets[4].absrap();
                histos["h_pt_jet5_5jet"]->fill(pT5);
                histos["h_y_jet5_5jet"]->fill(rap5);
                histos["h_HT_5jet"]->fill(HT);
                histos["h_ST_5jet"]->fill(ST);
              }
            }
          }
        }
      }
    }


    // Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve boson candidate
      const WFinder& wf = apply<WFinder>(event, "WF");
      if (wf.empty()) vetoEvent;

      // Retrieve jets
      const JetAlg& jetfs = apply<JetAlg>(event, "Jets");
      Jets all_jets = jetfs.jetsByPt(Cuts::pT > 30.0*GeV && Cuts::absrap < 4.4);

      const Particles& leptons = wf.constituentLeptons();
      const double missET = wf.constituentNeutrino().pT() / GeV;
      if (leptons.size() == 1 && missET > 25.0 && wf.mT() > 40.0*GeV) {
        const Particle& lep = leptons[0];
        fillPlots(lep, missET, all_jets);
      }
    }


    void finalize() {
      const double scalefactor(crossSection() / sumOfWeights());
      /// @todo Update to use C++11 range-for
      for (map<string, Histo1DPtr>::iterator hit = histos.begin(); hit != histos.end(); ++hit) {
        scale(hit->second, scalefactor);
      }
    }


  protected:

    size_t _mode;


  private:

    map<string, Histo1DPtr> histos;

  };


  class ATLAS_2014_I1319490_EL : public ATLAS_2014_I1319490 {
  public:
    ATLAS_2014_I1319490_EL()
      : ATLAS_2014_I1319490("ATLAS_2014_I1319490_EL")
    {
      _mode = 1;
    }
  };


  class ATLAS_2014_I1319490_MU : public ATLAS_2014_I1319490 {
  public:
    ATLAS_2014_I1319490_MU()
      : ATLAS_2014_I1319490("ATLAS_2014_I1319490_MU")
    {
      _mode = 2;
    }
  };


  // The hooks for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1319490);
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1319490_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1319490_MU);

}
