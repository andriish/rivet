// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Measurement of isolated diphoton + X differential cross-sections
  ///
  /// Inclusive isolated gamma gamma cross-sections, differential in M(gg), pT(gg),
  /// dphi(gg)
  ///
  /// @author Giovanni Marchiori
  class ATLAS_2011_S9120807 : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_S9120807()
      : Analysis("ATLAS_2011_S9120807"),
        _eta_bins_areaoffset{0.0, 1.5, 3.0}
    {    }


    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;
      addProjection(fs, "FS");

      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      addProjection(fj, "KtJetsD05");

      IdentifiedFinalState photonfs(Cuts::abseta < 2.37 && Cuts::pT > 16*GeV);
      photonfs.acceptId(PID::PHOTON);
      addProjection(photonfs, "Photon");

      _h_M    = bookHisto1D(1, 1, 1);
      _h_pT   = bookHisto1D(2, 1, 1);
      _h_dPhi = bookHisto1D(3, 1, 1);
    }


    /// @todo Prefer to use Rivet::binIndex()
    size_t getEtaBin(double eta_w) const {
      const double aeta = fabs(eta_w);
      size_t v_iter = 0;
      for (; v_iter+1 < _eta_bins_areaoffset.size(); ++v_iter) {
        if (inRange(aeta, _eta_bins_areaoffset[v_iter], _eta_bins_areaoffset[v_iter+1])) break;
      }
      return v_iter;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Require at least 2 photons in final state
      Particles photons = applyProjection<IdentifiedFinalState>(event, "Photon").particlesByPt();
      if (photons.size() < 2) vetoEvent;

      // Compute the median energy density
      vector<double> _ptDensity;
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);

      const shared_ptr<fastjet::ClusterSequenceArea> clust_seq_area = applyProjection<FastJets>(event, "KtJetsD05").clusterSeqArea();
      for (const Jet& jet : applyProjection<FastJets>(event, "KtJetsD05").jets()) {
        const double area = clust_seq_area->area(jet); // .pseudojet() called implicitly
        /// @todo Should be 1e-4?
        if (area > 10e-4 && jet.abseta() < _eta_bins_areaoffset.back()) {
          ptDensities.at(getEtaBin(jet.abseta())).push_back(jet.pT()/area);
        }
      }

      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; ++b) {
        double median = 0.0;
        if (!ptDensities[b].empty()) {
          std::sort(ptDensities[b].begin(), ptDensities[b].end());
          int nDens = ptDensities[b].size();
          if (nDens % 2 == 0) {
            median = (ptDensities[b][nDens/2] + ptDensities[b][(nDens-2)/2]) / 2;
          } else {
            median = ptDensities[b][(nDens-1)/2];
          }
        }
        _ptDensity.push_back(median);
      }

      // Loop over photons and fill vector of isolated ones
      Particles isolated_photons;
      for (const Particle& photon : photons) {

        // Remove photons in crack
        if (inRange(photon.abseta(), 1.37, 1.52)) continue;

        // Standard ET cone isolation
        const Particles& fs = applyProjection<FinalState>(event, "FS").particles();
        FourMomentum mom_in_EtCone;
        for (const Particle& p : fs) {
          // Check if it's in the cone of .4
          if (deltaR(photon, p) >= 0.4) continue;
          // Veto if it's in the 5x7 central core
          if (fabs(deltaEta(photon, p)) < 0.025*5.0*0.5 &&
              fabs(deltaPhi(photon, p)) < (M_PI/128.)*7.0*0.5) continue;
          // Increment isolation cone ET sum
          mom_in_EtCone += p.momentum();
        }

        // Now figure out the correction (area*density)
        const double EtCone_area = M_PI*.4*.4 - (7.0*.025)*(5.0*M_PI/128.);
        const double correction = _ptDensity[getEtaBin(photon.abseta())]*EtCone_area;

        // Shouldn't need to subtract photon
        // NB. Using expected cut at hadron/particle level, not cut at reco level
        if (mom_in_EtCone.Et() - correction > 4.0*GeV) continue;

        // Add to list of isolated photons
        isolated_photons.push_back(photon);
      }

      // Require at least two isolated photons
      if (isolated_photons.size() < 2) vetoEvent;

      // Select leading pT pair
      std::sort(isolated_photons.begin(), isolated_photons.end(), cmpMomByPt);
      const FourMomentum y1 = isolated_photons[0].momentum();
      const FourMomentum y2 = isolated_photons[1].momentum();

      // Require the two photons to be separated (dR>0.4)
      if (deltaR(y1, y2) < 0.4) vetoEvent;

      const FourMomentum yy = y1 + y2;
      const double Myy = yy.mass()/GeV;
      const double pTyy = yy.pT()/GeV;
      const double dPhiyy = deltaPhi(y1.phi(), y2.phi());

      const double weight = event.weight();
      _h_M->fill(Myy, weight);
      _h_pT->fill(pTyy, weight);
      _h_dPhi->fill(dPhiyy, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_M, crossSection()/sumOfWeights());
      scale(_h_pT, crossSection()/sumOfWeights());
      scale(_h_dPhi, crossSection()/sumOfWeights());
    }


  private:

    Histo1DPtr _h_M, _h_pT, _h_dPhi;
    vector<double> _eta_bins_areaoffset;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_S9120807);


}
