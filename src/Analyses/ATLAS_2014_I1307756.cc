// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2014_I1307756 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1307756()
      : Analysis("ATLAS_2014_I1307756"),
        _eta_bins_areaoffset{0.0, 1.5, 3.0}
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections here
      FinalState fs;
      addProjection(fs, "FS");

      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      addProjection(fj, "KtJetsD05");

      IdentifiedFinalState photonfs(Cuts::abseta < 2.37 && Cuts::pT > 22*GeV);
      photonfs.acceptId(PID::PHOTON);
      addProjection(photonfs, "photons");

      // Initialize event count here:
      _fidWeights = 0.;
    }


    /// Utility to compute ambiant energy density per eta bin
    /// @todo Use bin index lookup util instead...
    int getEtaBin(double eta_w) const {
      double eta = fabs(eta_w);
      int v_iter = 0;
      for (v_iter = 0; v_iter < (int)_eta_bins_areaoffset.size()-1; ++v_iter) {
        if (inRange(eta, _eta_bins_areaoffset[v_iter], _eta_bins_areaoffset[v_iter+1]))
          break;
      }
      return v_iter;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      /// Require at least 2 photons in final state
      Particles photons = applyProjection<IdentifiedFinalState>(event, "photons").particlesByPt();
      if (photons.size() < 2) vetoEvent;

      /// compute the median energy density per eta bin
      vector<double> ptDensity;
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);

      const auto clust_seq_area = applyProjection<FastJets>(event, "KtJetsD05").clusterSeqArea();
      for (const Jet& jet : applyProjection<FastJets>(event, "KtJetsD05").jets()) {
        const double area = clust_seq_area->area(jet);
        if (area > 1e-4 && jet.abseta() < _eta_bins_areaoffset.back()) {
          ptDensities.at(getEtaBin(jet.abseta())).push_back(jet.pT()/area);
        }
      }

      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; ++b) {
        double median = 0.0;
        if (ptDensities[b].size() > 0) {
          std::sort(ptDensities[b].begin(), ptDensities[b].end());
          const int nDens = ptDensities[b].size();
          if (nDens % 2 == 0) {
            median = (ptDensities[b][nDens/2] + ptDensities[b][(nDens-2)/2]) / 2;
          } else {
            median = ptDensities[b][(nDens-1)/2];
          }
        }
        ptDensity.push_back(median);
      }

      // Loop over photons and find isolated ones
      Particles isolated_photons;
      for (const Particle& ph : photons) {
        Particles fs = applyProjection<FinalState>(event, "FS").particles();
        FourMomentum mom_in_EtCone;
        for (const Particle& p : fs) {

          // Reject if the particle is not in DR=0.4 cone
          if (deltaR(ph.momentum(), p.momentum()) > 0.4) continue;

          // Reject if the particle falls in the photon core
          if (fabs(ph.eta() - p.eta()) < 0.025 * 7 * 0.5 &&
              fabs(ph.phi() - p.phi()) < PI/128. * 5 * 0.5) continue;

          // Reject if the particle is a neutrino (muons are kept)
          if (p.isNeutrino()) continue;

          // Sum momenta
          mom_in_EtCone += p.momentum();
        }

        // Subtract the UE correction (area*density)
        double EtCone_area = M_PI*.4*.4 - (7.0*.025)*(5.0*M_PI/128.);
        double correction = ptDensity[getEtaBin(ph.eta())]*EtCone_area;

        // Add isolated photon to list
        if (mom_in_EtCone.Et() - correction > 12*GeV) continue;
        isolated_photons.push_back(ph);
      }

      // Require at least two isolated photons
      if (isolated_photons.size() < 2)  vetoEvent ;

      // Select leading pT pair
      std::sort(isolated_photons.begin(), isolated_photons.end(), cmpMomByPt);
      FourMomentum y1 = isolated_photons[0].momentum();
      FourMomentum y2 = isolated_photons[1].momentum();

      // Compute invariant mass
      FourMomentum yy = y1 + y2;
      double Myy = yy.mass();

      // If Myy >= 110 GeV, apply relative cuts
      if (Myy/GeV >= 110 && (y1.Et()/Myy < 0.4 || y2.Et()/Myy < 0.3) ) vetoEvent;

      // Add to cross-section
      _fidWeights += event.weight();
    }


    /// @todo Add to the YODA output rather than print to log
    void finalize() {

      // Compute selection efficiency & statistical error
      const double eff = _fidWeights/sumOfWeights();
      const double err = sqrt(eff*(1-eff)/numEvents());

      // Compute fiducial cross-section in fb
      const double fidCrossSection = eff * crossSection()/femtobarn;

      // Print out result
      MSG_INFO("==================================================");
      MSG_INFO("==== Total cross-section: " << crossSection()/femtobarn<< " fb");
      MSG_INFO("==== Fiducial cross-section: " << fidCrossSection << " fb");
      MSG_INFO("==================================================");
      MSG_INFO("==== Selection efficiency: " << eff << " +/- " << err << " (statistical error)");
      MSG_INFO("==================================================");
    }

    //@}


  private:

    vector<double> _eta_bins_areaoffset;
    float _fidWeights;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1307756);

}
