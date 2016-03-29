// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Measurement of isolated diphoton + X differential cross-sections
  ///
  /// Inclusive isolated gamma gamma cross-sections, differential in M(gg), pT(gg),
  /// dphi(gg), cos(theta*)_CS
  ///
  /// @author Giovanni Marchiori
  ///
  class ATLAS_2012_I1199269 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1199269()
      : Analysis("ATLAS_2012_I1199269"),
        _eta_bins_areaoffset{ 0.0, 1.5, 3.0 }
    {    }


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      addProjection(fs, "FS");

      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      addProjection(fj, "KtJetsD05");

      IdentifiedFinalState photonfs(Cuts::abseta < 2.37 && Cuts::pT > 22*GeV);
      photonfs.acceptId(PID::PHOTON);
      addProjection(photonfs, "Photon");

      _h_M            = bookHisto1D(1, 1, 1);
      _h_pT           = bookHisto1D(2, 1, 1);
      _h_dPhi         = bookHisto1D(3, 1, 1);
      _h_cosThetaStar = bookHisto1D(4, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Require at least 2 photons in final state
      const Particles photons = applyProjection<IdentifiedFinalState>(event, "Photon").particlesByPt();
      if (photons.size() < 2) vetoEvent;

      // Get jets, and corresponding jet areas
      vector<vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);
      const auto clust_seq_area = applyProjection<FastJets>(event, "KtJetsD05").clusterSeqArea();
      for (const Jet& jet : applyProjection<FastJets>(event, "KtJetsD05").jets()) {
        const double area = clust_seq_area->area(jet); // implicit .pseudojet()
        if (area < 1e-3) continue;
        const int ieta = binIndex(jet.abseta(), _eta_bins_areaoffset);
        if (ieta != -1) ptDensities[ieta].push_back(jet.pT()/area);
      }

      // Compute median jet properties over the jets in the event
      vector<double> vptDensity; //, vsigma, vNjets;
      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; ++b) {
        double median = 0.0;
        // , sigma = 0.0;
        // int Njets = 0;
        if (ptDensities[b].size() > 0) {
          std::sort(ptDensities[b].begin(), ptDensities[b].end());
          int nDens = ptDensities[b].size();
          median = (nDens % 2 == 0) ? (ptDensities[b][nDens/2]+ptDensities[b][(nDens-2)/2])/2 : ptDensities[b][(nDens-1)/2];
          // sigma = ptDensities[b][(int)(.15865*nDens)];
          // Njets = nDens;
        }
        vptDensity.push_back(median);
        // vsigma.push_back(sigma);
        // vNjets.push_back(Njets);
      }


      // Loop over photons and fill vector of isolated ones
      Particles isolated_photons;
      for (const Particle& photon : photons) {
        /// Remove photons in ECAL crack region
        if (inRange(photon.abseta(), 1.37, 1.52)) continue;
        const double eta_P = photon.eta();
        const double phi_P = photon.phi();

        // Compute isolation via particles within an R=0.4 cone of the photon
        Particles fs = applyProjection<FinalState>(event, "FS").particles();
        FourMomentum mom_in_EtCone;
        for (const Particle& p : fs) {
          // Reject if not in cone
          if (deltaR(photon.momentum(), p.momentum()) > 0.4) continue;
          // Reject if in the 5x7 cell central core
          if (fabs(eta_P - p.eta()) < 0.025 * 5 * 0.5 &&
              fabs(phi_P - p.phi()) < PI/128. * 7 * 0.5) continue;
          // Sum momentum
          mom_in_EtCone += p.momentum();
        }
        // Now figure out the correction (area*density)
        const double EtCone_area = PI*sqr(0.4) - (7*.025)*(5*PI/128.); // cone area - central core rectangle
        const double correction = vptDensity[binIndex(fabs(eta_P), _eta_bins_areaoffset)] * EtCone_area;

        // Discard the photon if there is more than 4 GeV of cone activity
        // NOTE: Shouldn't need to subtract photon itself (it's in the central core)
        // NOTE: using expected cut at hadron/particle level, not at reco level
        if (mom_in_EtCone.Et() - correction > 4*GeV) continue;
        // Add isolated photon to list
        isolated_photons.push_back(photon);
      }

      // Require at least two isolated photons and select leading pT pair
      if (isolated_photons.size() < 2) vetoEvent;
      sortByPt(isolated_photons);
      FourMomentum y1 = isolated_photons[0].momentum();
      FourMomentum y2 = isolated_photons[1].momentum();

      // Leading photon should have pT > 25 GeV
      if (y1.pT() < 25*GeV) vetoEvent;

      // Require the two photons to be separated by dR > 0.4
      if (deltaR(y1, y2) < 0.4) vetoEvent;

      // Compute diphoton vector and fill histos
      const double weight = event.weight();
      FourMomentum yy = y1 + y2;
      const double costhetayy = 2 * y1.pT() * y2.pT() * sinh(y1.eta() - y2.eta()) / yy.mass() / add_quad(yy.mass(), yy.pT());
      _h_M->fill(yy.mass()/GeV, weight);
      _h_pT->fill(yy.pT()/GeV, weight);
      _h_dPhi->fill(mapAngle0ToPi(y1.phi() - y2.phi()), weight);
      _h_cosThetaStar->fill(costhetayy, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_M, crossSection()/sumOfWeights());
      scale(_h_pT, crossSection()/sumOfWeights());
      scale(_h_dPhi, crossSection()/sumOfWeights());
      scale(_h_cosThetaStar, crossSection()/sumOfWeights());
    }


  private:

    Histo1DPtr _h_M, _h_pT, _h_dPhi, _h_cosThetaStar;

    vector<double> _eta_bins, _eta_bins_areaoffset;

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1199269);

}
