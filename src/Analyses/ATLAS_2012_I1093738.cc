// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Measurement of isolated gamma + jet + X differential cross-sections
  ///
  /// Inclusive isolated gamma + jet cross-sections, differential in pT(gamma), for
  /// various photon and jet rapidity configurations.
  ///
  /// @author Giovanni Marchiori
  class ATLAS_2012_I1093738 : public Analysis {
  public:

    // Constructor
    ATLAS_2012_I1093738()
      : Analysis("ATLAS_2012_I1093738"),
        _eta_bins_ph{0.0, 1.37, 1.52, 2.37},
        _eta_bins_jet{0.0, 1.2, 2.8, 4.4}
    {    }


    // Book histograms and initialise projections before the run
    void init() {
      // Final state
      FinalState fs;
      addProjection(fs, "FS");

      // Voronoi eta-phi tessellation with KT jets, for ambient energy density calculation
      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      addProjection(fj, "KtJetsD05");

      // Leading photon
      LeadingParticlesFinalState photonfs(FinalState(-1.37, 1.37, 25.0*GeV));
      photonfs.addParticleId(PID::PHOTON);
      addProjection(photonfs, "LeadingPhoton");

      // FS excluding the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      addProjection(vfs, "JetFS");

      // Jets
      FastJets jetpro(vfs, FastJets::ANTIKT, 0.4);
      jetpro.useInvisibles();
      addProjection(jetpro, "Jets");

      _h_phbarrel_jetcentral_SS = bookHisto1D(1, 1, 1);
      _h_phbarrel_jetmedium_SS  = bookHisto1D(2, 1, 1);
      _h_phbarrel_jetforward_SS = bookHisto1D(3, 1, 1);

      _h_phbarrel_jetcentral_OS = bookHisto1D(4, 1, 1);
      _h_phbarrel_jetmedium_OS  = bookHisto1D(5, 1, 1);
      _h_phbarrel_jetforward_OS = bookHisto1D(6, 1, 1);
    }


    /// @todo Use Rivet::binIndex()
    int getEtaBin(double eta_w, int what) const {
      double eta = fabs(eta_w);
      int v_iter = 0;
      if (what==0) {
        for (v_iter=0; v_iter < (int)_eta_bins_ph.size()-1; v_iter++){
          if (eta >= _eta_bins_ph.at(v_iter) && eta < _eta_bins_ph.at(v_iter+1))
            break;
        }
      }
      else if (what==1) {
        for (v_iter=0; v_iter < (int)_eta_bins_jet.size()-1; v_iter++){
          if (eta >= _eta_bins_jet.at(v_iter) && eta < _eta_bins_jet.at(v_iter+1))
            break;
        }
      }
      else {
        for (v_iter=0; v_iter < (int)_eta_bins_areaoffset.size()-1; v_iter++){
          if (eta >= _eta_bins_areaoffset.at(v_iter) && eta < _eta_bins_areaoffset.at(v_iter+1))
            break;
        }
      }
      return v_iter;
    }


    // Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the photon
      const FinalState& photonfs = applyProjection<FinalState>(event, "LeadingPhoton");
      if (photonfs.particles().size() < 1) vetoEvent;
      const FourMomentum photon = photonfs.particles().front().momentum();

      // Get the jet
      Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(20.0*GeV);
      if (jets.empty()) vetoEvent;
      FourMomentum leadingJet = jets[0].momentum();

      // Require jet separated from photon
      if (deltaR(photon, leadingJet) < 1.0) vetoEvent;

      // Veto if leading jet is outside plotted rapidity regions
      if (leadingJet.absrap() > 4.4) vetoEvent;

      // Compute the median event energy density
      const unsigned int skipnhardjets = 0;
      _ptDensity.clear();
      _sigma.clear();
      _Njets.clear();
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);

      FastJets fastjets = applyProjection<FastJets>(event, "KtJetsD05");
      const shared_ptr<fastjet::ClusterSequenceArea> clust_seq_area = fastjets.clusterSeqArea();
      for (const Jet& jet : fastjets.jets()) {
        const double area = clust_seq_area->area(jet); //< Implicit call to pseudojet()
        /// @todo Should be 1e-4?
        if (area > 10e-4 && jet.abseta() < _eta_bins_areaoffset.back()) {
          ptDensities.at(getEtaBin(jet.abseta(),2)).push_back(jet.pT()/area);
        }
      }

      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; b++) {
        double median = 0.0;
        double sigma = 0.0;
        int Njets = 0;
        if (ptDensities[b].size() > skipnhardjets) {
          std::sort(ptDensities[b].begin(), ptDensities[b].end());
          int nDens = ptDensities[b].size() - skipnhardjets;
          if (nDens % 2 == 0) {
            median = (ptDensities[b][nDens/2]+ptDensities[b][(nDens-2)/2])/2;
          } else {
            median = ptDensities[b][(nDens-1)/2];
          }
          sigma = ptDensities[b][(int)(.15865*nDens)];
          Njets = nDens;
        }
        _ptDensity.push_back(median);
        _sigma.push_back(sigma);
        _Njets.push_back(Njets);
      }

      // Compute photon isolation with a standard ET cone
      const Particles fs = applyProjection<FinalState>(event, "JetFS").particles();
      FourMomentum mom_in_EtCone;
      const double ISO_DR = 0.4;
      const double CLUSTER_ETA_WIDTH = 0.25*5.0;
      const double CLUSTER_PHI_WIDTH = (PI/128.)*7.0;
      for (const Particle& p : fs) {
        // Check if it's in the cone of .4
        if (deltaR(photon, p) >= ISO_DR) continue;
        // Check if it's in the 5x7 central core
        if (fabs(deltaEta(photon, p)) < CLUSTER_ETA_WIDTH*0.5 &&
            fabs(deltaPhi(photon, p)) < CLUSTER_PHI_WIDTH*0.5) continue;
        // Increment sum
        mom_in_EtCone += p.momentum();
      }

      // Figure out the correction (area*density)
      const double EtCone_area = PI*ISO_DR*ISO_DR - CLUSTER_ETA_WIDTH*CLUSTER_PHI_WIDTH;
      const double correction = _ptDensity[getEtaBin(photon.abseta(),2)]*EtCone_area;

      // Require photon to be isolated
      if (mom_in_EtCone.Et()-correction > 4.0*GeV) vetoEvent;

      const int photon_jet_sign = sign( leadingJet.rapidity() * photon.rapidity() );

      // Fill histos
      const double abs_jet_rapidity = fabs(leadingJet.rapidity());
      const double photon_pt = photon.pT()/GeV;
      const double abs_photon_eta = fabs(photon.eta());
      const double weight = event.weight();
      if (abs_photon_eta < 1.37) {
        if (abs_jet_rapidity < 1.2) {
          if (photon_jet_sign >= 1) {
            _h_phbarrel_jetcentral_SS->fill(photon_pt, weight);
          } else {
            _h_phbarrel_jetcentral_OS->fill(photon_pt, weight);
          }
        } else if (abs_jet_rapidity < 2.8) {
          if (photon_jet_sign >= 1) {
            _h_phbarrel_jetmedium_SS->fill(photon_pt, weight);
          } else {
            _h_phbarrel_jetmedium_OS->fill(photon_pt, weight);
          }
        } else if (abs_jet_rapidity < 4.4) {
          if (photon_jet_sign >= 1) {
            _h_phbarrel_jetforward_SS->fill(photon_pt, weight);
          } else {
            _h_phbarrel_jetforward_OS->fill(photon_pt, weight);
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_phbarrel_jetcentral_SS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetcentral_OS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetmedium_SS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetmedium_OS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetforward_SS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetforward_OS, crossSection()/sumOfWeights());
    }


  private:

    Histo1DPtr _h_phbarrel_jetcentral_SS, _h_phbarrel_jetmedium_SS, _h_phbarrel_jetforward_SS;
    Histo1DPtr _h_phbarrel_jetcentral_OS, _h_phbarrel_jetmedium_OS, _h_phbarrel_jetforward_OS;

    vector<float> _eta_bins_ph, _eta_bins_jet, _eta_bins_areaoffset;
    vector<float> _ptDensity, _sigma, _Njets;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1093738);


}
