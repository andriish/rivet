// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2013_I1243871 : public Analysis {
  public:

    /// Constructor
    ATLAS_2013_I1243871()
      : Analysis("ATLAS_2013_I1243871")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {
      // Set up projections
      const FinalState fs(-4.5, 4.5);
      declare(fs, "ALL_FS");

      /// Get electrons from truth record
      IdentifiedFinalState elec_fs(Cuts::abseta < 2.47 && Cuts::pT > 25*GeV);
      elec_fs.acceptIdPair(PID::ELECTRON);
      declare(elec_fs, "ELEC_FS");

      /// Get muons which pass the initial kinematic cuts:
      IdentifiedFinalState muon_fs(Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      muon_fs.acceptIdPair(PID::MUON);
      declare(muon_fs, "MUON_FS");

      /// Get b-hadrons for tagging
      HeavyHadrons hh(Cuts::pT > 5*GeV);
      declare(hh, "HF_HADRONS");

      // Final state used as input for jet-finding.
      // We include everything except the muons and neutrinos
      VetoedFinalState jet_input(fs);
      jet_input.vetoNeutrinos();
      jet_input.addVetoPairId(PID::MUON);
      declare(jet_input, "JET_INPUT");

      // Get the jets
      FastJets jets(jet_input, FastJets::ANTIKT, 0.4);
      declare(jets, "JETS");

      // Book histograms
      for (size_t d = 0; d < 5; ++d) {
        _p_b_rho[d] = bookProfile1D(d+1, 1, 1);
        _p_l_rho[d] = bookProfile1D(d+1, 2, 1);
        _p_b_Psi[d] = bookProfile1D(d+1, 1, 2);
        _p_l_Psi[d] = bookProfile1D(d+1, 2, 2);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// Get the various sets of final state particles
      const Particles& elecs = apply<IdentifiedFinalState>(event, "ELEC_FS").particlesByPt();
      const Particles& muons = apply<IdentifiedFinalState>(event, "MUON_FS").particlesByPt();

      // Get all jets with pT > 7 GeV (ATLAS standard jet collection)
      const Jets& allJets = apply<FastJets>(event, "JETS").jetsByPt(7*GeV);

      // Keep any jets that pass the pt cut
      const Jets pt_jets = filter_select(allJets, Cuts::pt > 25*GeV && Cuts::abseta < 2.5);

      // Remove jets too close to an electron
      const Jets good_jets = filter_discard(pt_jets, [&](const Jet& j){
          return any(elecs, deltaRLess(j, 0.2));
        });

      // Classify the event type
      const size_t nElec = elecs.size();
      const size_t nMuon = muons.size();
      bool isSemilepton = false, isDilepton = false;
      if (nElec == 1 && nMuon == 0) {
        isSemilepton = true;
      } else if (nElec == 0 && nMuon == 1) {
        isSemilepton = true;
      } else if (nElec == 2 && nMuon == 0) {
        if (charge(elecs[0]) != charge(elecs[1])) isDilepton = true;
      } else if (nElec == 1 && nMuon == 1) {
        if (charge(elecs[0]) != charge(muons[0])) isDilepton = true;
      } else if (nElec == 0 && nMuon == 2) {
        if (charge(muons[0]) != charge(muons[1])) isDilepton = true;
      }
      const bool isGoodEvent = (isSemilepton && good_jets.size() >= 4) || (isDilepton && good_jets.size() >= 2);
      if (!isGoodEvent) vetoEvent;

      // Weakly-decaying b-hadrons for tagging
      const Particles& b_hadrons = apply<HeavyHadrons>(event, "HF_HADRONS").bHadrons();

      // Select b-jets as those containing a b-hadron and not overlapped with any other jet
      Jets b_jets = filter_select(good_jets, [&](const Jet& j){
          if (all(b_hadrons, deltaRGtr(j, 0.3))) return false;
          for (const Jet& k : allJets)
            if (inRange(deltaR(j,k), 0.01, 0.8)) return false;
          return true;
        });
      MSG_DEBUG(b_jets.size() << " b-jets selected");


      // Select light-jets as the pair of non-b-jets with invariant mass closest to the W mass
      /// @todo Use built-in b-tagging (dR < 0.3 defn), avoid HepMC
      const double nominalW = 80.4*GeV;
      double deltaM = 500*GeV;
      const Jet* light1 = NULL; const Jet* light2 = NULL; // NB: const Jets, not const pointers!
      for (const Jet& i : good_jets) {
        const bool isbJet1 = any(b_hadrons, deltaRLess(i, 0.3));
        if (isbJet1) continue;
        for (const Jet& j : good_jets) {
          const bool isbJet2 = any(b_hadrons, deltaRLess(j, 0.3));
          if (isbJet2) continue;
          const double invMass = (i.momentum()+j.momentum()).mass();
          const double dM = fabs(invMass - nominalW);
          if (dM < deltaM) {
            deltaM = dM;
            light1 = &i;
            light2 = &j;
          }
        }
      }

      // Check that both jets are not overlapped, and populate the light jets list
      Jets light_jets;
      const bool hasGoodLight = light1 != NULL && light2 != NULL && light1 != light2;
      if (hasGoodLight) {
        bool isOverlap1 = false, isOverlap2 = false;
        for (const Jet& j : allJets) {
          if (light1 == &j) continue;
          const double dR1j = deltaR(*light1, j);
          if (dR1j < 0.8) { isOverlap1 = true; break; }
        }
        for (const Jet& j : allJets) {
          if (light2 == &j) continue;
          const double dR2j = deltaR(*light2, j);
          if (dR2j < 0.8) { isOverlap2 = true; break; }
        }
        if (!isOverlap1 && !isOverlap2) {
          light_jets = {*light1, *light2};
        }
      }
      MSG_DEBUG(light_jets.size() << " light jets selected");


      // Calculate the jet shapes
      const double binWidth = 0.04; // -> 10 bins from 0.0-0.4
      const vector<double> ptEdges = {{ 30, 40, 50, 70, 100, 150 }};

      // b-jet shapes
      MSG_DEBUG("Filling b-jet shapes");
      for (const Jet& bJet : b_jets) {
        // Work out jet pT bin and skip this jet if out of range
        const double jetPt = bJet.pT();
        MSG_DEBUG("Jet pT = " << jetPt/GeV << " GeV");
        if (!inRange(jetPt/GeV, 30., 150.)) continue;
        /// @todo Use YODA bin index lookup tools
        size_t ipt; for (ipt = 0; ipt < 5; ++ipt) if (inRange(jetPt/GeV, ptEdges[ipt], ptEdges[ipt+1])) break;
        MSG_DEBUG("Jet pT index = " << ipt);

        // Calculate jet shape
        vector<double> rings(10, 0);
        for (const Particle& p : bJet.particles()) {
          const double dR = deltaR(bJet, p);
          const size_t idR = (size_t) floor(dR/binWidth);
          for (size_t i = idR; i < 10; ++i) rings[i] += p.pT();
        }

        // Fill each dR bin of the histos for this jet pT
        for (int iBin = 0; iBin < 10; ++iBin) {
          const double rcenter = 0.02 + iBin*binWidth;
          const double rhoval = (iBin != 0 ? (rings[iBin]-rings[iBin-1]) : rings[iBin]) / binWidth / rings[9];
          const double psival = rings[iBin] / rings[9];
          MSG_DEBUG(rcenter << ", " << rhoval << ", " << psival);
          _p_b_rho[ipt]->fill(rcenter, rhoval, weight);
          _p_b_Psi[ipt]->fill(rcenter, psival, weight);
        }
      }

      // Light jet shapes
      MSG_DEBUG("Filling light jet shapes");
      for (const Jet& lJet : light_jets) {
        // Work out jet pT bin and skip this jet if out of range
        const double jetPt = lJet.pT();
        MSG_DEBUG("Jet pT = " << jetPt/GeV << " GeV");
        if (!inRange(jetPt/GeV, 30., 150.)) continue;
        /// @todo Use YODA bin index lookup tools
        size_t ipt; for (ipt = 0; ipt < 5; ++ipt) if (inRange(jetPt/GeV, ptEdges[ipt], ptEdges[ipt+1])) break;
        MSG_DEBUG("Jet pT index = " << ipt);

        // Calculate jet shape
        vector<double> rings(10, 0);
        for (const Particle& p : lJet.particles()) {
          const double dR = deltaR(lJet, p);
          const size_t idR = (size_t) floor(dR/binWidth);
          for (size_t i = idR; i < 10; ++i) rings[i] += p.pT();
        }

        // Fill each dR bin of the histos for this jet pT
        for (int iBin = 0; iBin < 10; ++iBin) {
          const double rcenter = 0.02 + iBin*binWidth;
          const double rhoval = (iBin != 0 ? (rings[iBin]-rings[iBin-1]) : rings[iBin]) / binWidth / rings[9];
          const double psival = rings[iBin] / rings[9];
          _p_l_rho[ipt]->fill(rcenter, rhoval, weight);
          _p_l_Psi[ipt]->fill(rcenter, psival, weight);
        }
      }

    }


  private:

    Profile1DPtr _p_b_rho[5];
    Profile1DPtr _p_l_rho[5];
    Profile1DPtr _p_b_Psi[5];
    Profile1DPtr _p_l_Psi[5];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1243871);

}
