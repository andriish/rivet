// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetYODA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeptonClusters.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Particle.hh"
#include "HepMC/GenEvent.h"

namespace Rivet {

  struct ATLAS_2012_I1094568_plots {
    // Keep track of which veto region this is, to match
    // the autobook-ed histograms
    int region_index;

    // lower rapidity boundary or veto region
    double y_low;
    // upper rapidity boundary or veto region
    double y_high;

    double vetoJetPt_Q0;
    double vetoJetPt_Qsum;

    // Histograms to store the veto jet pT and
    // sum(veto jet pT) histograms.
    Histo1DPtr _h_vetoJetPt_Q0;
    Histo1DPtr _h_vetoJetPt_Qsum;

    // DataPointSets for the gap fractions
    Scatter2DPtr _d_gapFraction_Q0;
    Scatter2DPtr _d_gapFraction_Qsum;
  };


  class ATLAS_2012_I1094568 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1094568() : Analysis("ATLAS_2012_I1094568")
    {}


  public:

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(-4.5, 4.5);
      addProjection(fs, "ALL_FS");

      /// Get electrons from truth record
      IdentifiedFinalState elec_fs(-2.47, 2.47, 25.0*GeV);
      elec_fs.acceptIdPair(ELECTRON);
      addProjection(elec_fs, "ELEC_FS");

      /// Get muons which pass the initial kinematic cuts:
      IdentifiedFinalState muon_fs(-2.5, 2.5, 20.0*GeV);
      muon_fs.acceptIdPair(MUON);
      addProjection(muon_fs, "MUON_FS");

      /// Get all neutrinos. These will not be used to form jets.
      /// We'll use the highest 2 pT neutrinos to calculate the MET
      IdentifiedFinalState neutrino_fs(-4.5, 4.5, 0.0*GeV);
      neutrino_fs.acceptNeutrinos();
      addProjection(neutrino_fs, "NEUTRINO_FS");

      // Final state used as input for jet-finding.
      // We include everything except the muons and neutrinos
      VetoedFinalState jet_input(fs);
      jet_input.vetoNeutrinos();
      jet_input.addVetoPairId(MUON);
      addProjection(jet_input, "JET_INPUT");

      // Get the jets
      FastJets jets(jet_input, FastJets::ANTIKT, 0.4);
      addProjection(jets, "JETS");

      for(int i=0; i<201; ++i) {
        double bin_edge = i*5;
        m_q0BinEdges += bin_edge;
      }

      m_total_weight = 0.0;

      m_plots[0].region_index = 1;
      m_plots[0].y_low = 0.0;
      m_plots[0].y_high = 0.8;
      InitializePlots(m_plots[0]);

      m_plots[1].region_index = 2;
      m_plots[1].y_low = 0.8;
      m_plots[1].y_high = 1.5;
      InitializePlots(m_plots[1]);

      m_plots[2].region_index = 3;
      m_plots[2].y_low = 1.5;
      m_plots[2].y_high = 2.1;
      InitializePlots(m_plots[2]);

      m_plots[3].region_index = 4;
      m_plots[3].y_low = 0.0;
      m_plots[3].y_high = 2.1;
      InitializePlots(m_plots[3]);
    }

    void InitializePlots(ATLAS_2012_I1094568_plots& plots) {
      int q0_index = 1;
      int qsum_index = 2;

      std::stringstream vetoPt_Q0_name;
      vetoPt_Q0_name << "vetoJetPt_Q0_" << plots.region_index;

      std::stringstream vetoPt_Qsum_name;
      vetoPt_Qsum_name << "vetoJetPt_Qsum_" << plots.region_index;

      plots._h_vetoJetPt_Q0   = bookHisto1D(vetoPt_Q0_name.str(), m_q0BinEdges);
      plots._h_vetoJetPt_Qsum = bookHisto1D(vetoPt_Qsum_name.str(), m_q0BinEdges);

      plots._d_gapFraction_Q0   = bookScatter2D(plots.region_index, q0_index, 1);
      plots._d_gapFraction_Qsum = bookScatter2D(plots.region_index, qsum_index, 1);

      plots.vetoJetPt_Q0 = 0.0;
      plots.vetoJetPt_Qsum = 0.0;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      /// Get the various sets of final state particles
      const ParticleVector& elecFS = applyProjection<IdentifiedFinalState>(event, "ELEC_FS").particlesByPt();
      const ParticleVector& muonFS = applyProjection<IdentifiedFinalState>(event, "MUON_FS").particlesByPt();
      const ParticleVector& neutrinoFS = applyProjection<IdentifiedFinalState>(event, "NEUTRINO_FS").particlesByPt();

      // Get all jets with pT > 25 GeV
      const Jets& jets = applyProjection<FastJets>(event, "JETS").jetsByPt(25.0*GeV);

      // Keep any jets that pass the initial rapidity cut
      vector<const Jet*> central_jets;
      double rapMax = 2.4;
      foreach(const Jet& j, jets) {
        double rapidity = fabs(j.momentum().rapidity());
        if(rapidity < rapMax) central_jets.push_back(&j);
      }

      // For each of the jets that pass the rapidity cut, only keep those that are not
      // too close to any leptons
      vector<const Jet*> good_jets;
      foreach(const Jet* j, central_jets) {
        bool goodJet = true;
        foreach(const Particle& e, elecFS) {
          double elec_jet_dR = deltaR(e.momentum(), j->momentum());
          if(elec_jet_dR < 0.4) goodJet = false;
        }
        foreach(const Particle& m, muonFS) {
          double muon_jet_dR = deltaR(m.momentum(), j->momentum());
          if(muon_jet_dR < 0.4) goodJet = false;
        }
        if(goodJet == true) good_jets.push_back(j);
      }

      // Temporary fix to get B-hadrons in evgen files where they don't show up
      // in the UnstableFinalState projection
      // (e.g. mc10_7TeV.105200.T1_McAtNlo_Jimmy.evgen.EVNT.e598/)
      // This will be updated for MC12 to just use UnstableFinalState
      // (Thanks to Steve Bieniek for this!)
      std::vector<HepMC::GenParticle*> B_hadrons;
      std::vector<HepMC::GenParticle*> allParticles = particles(event.genEvent());
      for(unsigned int i = 0; i < allParticles.size(); i++) {
        GenParticle* p = allParticles.at(i);
        if ( !(Rivet::PID::isHadron( p->pdg_id() ) && Rivet::PID::hasBottom( p->pdg_id() )) ) continue;
        if(p->momentum().perp()*GeV < 5) continue;
        B_hadrons.push_back(p);
      }


      // For each of the good jets, check whether any are b-jets
      vector<const Jet*> b_jets;
      foreach(const Jet* j, good_jets) {
        bool isbJet = false;
        foreach(HepMC::GenParticle* b, B_hadrons) {
          FourMomentum hadron = b->momentum();
          double hadron_jet_dR = deltaR(j->momentum(), hadron);
          if(hadron_jet_dR < 0.3) isbJet = true;
        }
        if(isbJet) b_jets.push_back(j);
      }


      // Check the good jets again and keep track of the "additional jets"
      // I.e. those which are not either of the 2 highest pT b-jets
      vector<const Jet*> veto_jets;
      int n_bjets_matched = 0;
      foreach(const Jet* j, good_jets) {
        bool isBJet = false;
        foreach(const Jet* b, b_jets) {
          if(n_bjets_matched == 2) break;
          if(b == j){isBJet = true; ++ n_bjets_matched;}
        }
        if(!isBJet) veto_jets.push_back(j);
      }


      // Get the MET by taking the vector sum of all neutrinos
      double MET = 0;
      FourMomentum p_MET(0., 0., 0., 0.);
      foreach(const Particle& p, neutrinoFS) {
        p_MET = p_MET + p.momentum();
      }
      MET = p_MET.pT();


      // Now we have everything we need to start doing the event selections
      bool passed_ee = false;
      vector<const Jet*> vetoJets_ee;

      // We want exactly 2 electrons...
      if(elecFS.size() == 2) {
        // ... with opposite sign charges.
        if(PID::charge(elecFS.at(0)) != PID::charge(elecFS.at(1))) {
          // Check the MET
          if(MET >= 40*GeV) {
            // Do some dilepton mass cuts
            double dilepton_mass = (elecFS.at(0).momentum() + elecFS.at(1).momentum()).mass();
            if(dilepton_mass >= 15*GeV) {
              if(fabs(dilepton_mass - 91.0*GeV) >= 10.0*GeV) {
                // we need at least 2 b-jets
                if(b_jets.size() > 1) {
                  // This event has passed all the cuts;
                  passed_ee = true;
                }
              }
            }
          }
        }
      }


      bool passed_mumu = false;
      // Now do the same checks for the mumu channel
      vector<const Jet*> vetoJets_mumu;
      // So we now want 2 good muons...
      if(muonFS.size() == 2) {
        // ...with opposite sign charges.
        if(PID::charge(muonFS.at(0)) != PID::charge(muonFS.at(1))) {
          // Check the MET
          if(MET >= 40*GeV) {
            // and do some di-muon mass cuts
            double dilepton_mass = (muonFS.at(0).momentum() + muonFS.at(1).momentum()).mass();
            if(dilepton_mass >= 15*GeV) {
              if(fabs(dilepton_mass - 91.0*GeV) >= 10.0*GeV) {
                // Need atleast 2 b-jets
                if(b_jets.size() > 1) {
                  // This event has passed all mumu-channel cuts
                  passed_mumu = true;
                }
              }
            }
          }
        }
      }


      bool passed_emu = false;
      // Finally, the same again with the emu channel
      vector<const Jet*> vetoJets_emu;
      // We want exactly 1 electron and 1 muon
      if(elecFS.size() == 1 && muonFS.size() == 1) {
        // With opposite sign charges
        if(PID::charge(elecFS.at(0)) != PID::charge(muonFS.at(0))) {
          // Calculate the HT from the scalar sum of the pT of the leptons
          // and all good jets
          double HT = 0;
          HT += fabs(elecFS.at(0).momentum().pT());
          HT += fabs(muonFS.at(0).momentum().pT());
          foreach(const Jet* j, good_jets) {
            HT += fabs(j->momentum().pT());
          }
          // Keep events with HT > 130 GeV
          if(HT > 130.0*GeV) {
            // And again we want 2 or more b-jets
            if(b_jets.size() > 1) {
              passed_emu = true;
            }
          }
        }
      }

      if(passed_ee == true || passed_mumu == true || passed_emu == true) {
        // If the event passes the selection, we use it for all gap fractions
        m_total_weight += weight;

        // Loop over each veto jet
        foreach(const Jet* j, veto_jets) {
          double pt = j->momentum().pT();
          double rapidity = fabs(j->momentum().rapidity());
          // Loop over each region
          for(int i=0; i<4; ++i) {
            // If the jet falls into this region, get its pT and increment sum(pT)
            if( (rapidity > m_plots[i].y_low) && (rapidity < m_plots[i].y_high)) {
              m_plots[i].vetoJetPt_Qsum += pt;

              // If we've already got a veto jet, don't replace it
              if(m_plots[i].vetoJetPt_Q0 == 0.0) m_plots[i].vetoJetPt_Q0 = pt;
            }
          }
        } // end loop over veto jets
        for(int i=0; i<4; ++i) {
          m_plots[i]._h_vetoJetPt_Q0->fill(m_plots[i].vetoJetPt_Q0, weight);
          m_plots[i]._h_vetoJetPt_Qsum->fill(m_plots[i].vetoJetPt_Qsum, weight);
          ClearVetoJetPts(m_plots[i]);
        }
      }
    }

    void ClearVetoJetPts(ATLAS_2012_I1094568_plots& plots) {
      plots.vetoJetPt_Q0 = 0.0;
      plots.vetoJetPt_Qsum = 0.0;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(int i=0; i<4; ++i) {
        // @todo YODA
        //FinalizeGapFraction(m_total_weight, m_plots[i]._d_gapFraction_Q0, m_plots[i]._h_vetoJetPt_Q0, binEdges(i+1, 1, 1));
        //FinalizeGapFraction(m_total_weight, m_plots[i]._d_gapFraction_Qsum, m_plots[i]._h_vetoJetPt_Qsum, binEdges(i+1, 2, 1));
      }
    }

    // @todo YODA
    ////void FinalizeGapFraction(double total_weight, ATLAS_2011_I1094568_plots& plots, int type)
    //void FinalizeGapFraction(double total_weight, Scatter2DPtr gapFraction, Histo1DPtr vetoPt, BinEdges fgap_binEdges) {

    //  // Stores the cumulative frequency of the veto jet pT histogram
    //  double vetoPtWeightSum = 0.0;
    //  
    //  // Keep track of which gap fraction point we're doing
    //  unsigned int fgap_point = 0;
    //  for(unsigned int i=0; i<m_q0BinEdges.size()-2; ++i) {
    //    vetoPtWeightSum += vetoPt->binHeight(i);

    //    // If we've done the last point, stop.
    //    if(fgap_point == fgap_binEdges.size()-1) break;

    //    // Get the x-value of this gap fraction point, from the mid-point of the bin edges
    //    double binCentre = ( fgap_binEdges.at(fgap_point) + fgap_binEdges.at(fgap_point+1) ) / 2;
    //    double errorPlus = fgap_binEdges.at(fgap_point+1) - binCentre;
    //    double errorMinus = binCentre - fgap_binEdges.at(fgap_point);

    //    // If this Q0/Qsum point is not the cut value we need for this gap fraction point, continue
    //    if(m_q0BinEdges.at(i+1) != binCentre) continue;

    //    // Calculate the gap fraction and its uncertainty
    //    double fraction = vetoPtWeightSum/total_weight;
    //    double fraction_error = sqrt(fraction*(1.0-fraction)/total_weight);
    //    if(total_weight == 0.0) fraction = fraction_error = 0.0;

    //    // Set the point
    //    IDataPoint* currentPoint = gapFraction->point(fgap_point);
    //    IMeasurement* xCoord = currentPoint->coordinate(0);
    //    IMeasurement* yCoord = currentPoint->coordinate(1);

    //    xCoord->setValue(binCentre);
    //    xCoord->setErrorPlus(errorPlus);
    //    xCoord->setErrorMinus(errorMinus);
    //    yCoord->setValue(fraction);
    //    yCoord->setErrorPlus(fraction_error);
    //    yCoord->setErrorMinus(fraction_error);

    //    ++fgap_point;
    //  }
    //  tree().rm(tree().findPath(dynamic_cast<AIDA::IManagedObject&>(*vetoPt)));
    //}


  private:

    // define the vetoJet pT binning
    std::vector<double> m_q0BinEdges;
    double m_total_weight;


  private:
    ATLAS_2012_I1094568_plots m_plots[4];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1094568);

}
