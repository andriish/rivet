// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include <fastjet/ClusterSequence.hh>


namespace Rivet {


  class ATLAS_2012_I1094564 : public Analysis {
  public:

    ATLAS_2012_I1094564()
      : Analysis("ATLAS_2012_I1094564")
    {}



  public:

    // Returns constituents to make it easier to do the filtering
    PseudoJets splitjet(fastjet::PseudoJet jet,
                        double& last_R,
                        const FastJets& fj,
                        bool& unclustered) const {

      // Build a new cluster sequence just using the constituents of this jet.
      fastjet::ClusterSequence cs(fj.clusterSeq()->constituents(jet), fastjet::JetDefinition(fastjet::cambridge_algorithm, M_PI/2.));

      // Get the jet back again
      vector<fastjet::PseudoJet> remadeJets = cs.inclusive_jets(0.);

      if ( remadeJets.size() != 1 ) return remadeJets;

      fastjet::PseudoJet remadeJet = remadeJets[0];
      fastjet::PseudoJet parent1, parent2;
      unclustered = false;

      while ( cs.has_parents(remadeJet, parent1, parent2) ) {

        if (parent1.squared_distance(parent2) < 0.09) {
          break;
        }

        if (parent1.m2() < parent2.m2()) {
          fastjet::PseudoJet tmp;
          tmp = parent1; parent1 = parent2; parent2 = tmp;
        }

        double ktdist = parent1.kt_distance(parent2);
        double rtycut2 = 0.3*0.3;

        if (parent1.m() < 0.67*remadeJet.m() && ktdist > rtycut2*remadeJet.m2()) {
          unclustered = true;
          break;
        } else {
          remadeJet = parent1;
        }
      }

      last_R = 0.5 * sqrt(parent1.squared_distance(parent2));

      return cs.constituents(remadeJet);
    }


    fastjet::PseudoJet filterjet(PseudoJets jets,
                                 double& stingy_R,
                                 const double def_R) const {
      if (stingy_R == 0.0) {
        stingy_R = def_R;
      }

      stingy_R = def_R < stingy_R ? def_R : stingy_R;

      fastjet::JetDefinition stingy_jet_def(fastjet::cambridge_algorithm, stingy_R);

      fastjet::ClusterSequence scs(jets, stingy_jet_def);
      std::vector<fastjet::PseudoJet> stingy_jets = sorted_by_pt(scs.inclusive_jets(0.));

      fastjet::PseudoJet reconst_jet(0.0, 0.0, 0.0, 0.0);

      for (unsigned isj = 0; isj < std::min((unsigned int) 3, (unsigned int) stingy_jets.size()); ++isj) {
        reconst_jet += stingy_jets[isj];
      }
      return reconst_jet;
    }


    // These are custom functions for n-subjettiness.
    PseudoJets jetGetAxes(unsigned int n_jets,
                          PseudoJets& inputJets,
                          double subR) const {

      //sanity check
      if (inputJets.size() < n_jets) {
        std::cout << "Not enough input particles." << endl;
        return inputJets;
      }

      //get subjets, return
      fastjet::ClusterSequence sub_clust_seq(inputJets, fastjet::JetDefinition(fastjet::kt_algorithm, subR, fastjet::E_scheme, fastjet::Best));

      return sub_clust_seq.exclusive_jets((signed)n_jets);
    }


    double jetTauValue(double beta, double jet_rad,
                       PseudoJets& particles, PseudoJets& axes, double Rcut) const {

      double tauNum = 0.0;
      double tauDen = 0.0;

      if (particles.size() == 0) return 0.0;

      for (unsigned int i = 0; i < particles.size(); i++) {
        // find minimum distance (set R large to begin)
        double minR = 10000.0;
        for (unsigned int j = 0; j < axes.size(); j++) {
          double tempR = sqrt(particles[i].squared_distance(axes[j]));
          if (tempR < minR) minR = tempR;
        }
        if (minR > Rcut) minR = Rcut;
        //calculate nominator and denominator
        tauNum += particles[i].perp() * pow(minR,beta);
        tauDen += particles[i].perp() * pow(jet_rad,beta);
      }

      //return N-subjettiness
      if (tauDen == 0) return 0;
      else return tauNum/tauDen;
    }


    void jetUpdateAxes(double beta,
                       PseudoJets& particles, PseudoJets& axes) const {

      vector<int> belongsto;
      for (unsigned int i = 0; i < particles.size(); i++) {
        // find minimum distance axis
        int assign = 0;
        double minR = 10000.0;
        for (unsigned int j = 0; j < axes.size(); j++) {
          double tempR = sqrt(particles[i].squared_distance(axes[j]));
          if (tempR < minR) {
            minR = tempR;
            assign = j;
          }
        }
        belongsto.push_back(assign);
      }

      // iterative step
      double deltaR2, distphi;
      vector<double> ynom, phinom, den;

      ynom.resize(axes.size());
      phinom.resize(axes.size());
      den.resize(axes.size());

      for (unsigned int i = 0; i < particles.size(); i++) {

        distphi = particles[i].phi() - axes[belongsto[i]].phi();
        deltaR2 = particles[i].squared_distance(axes[belongsto[i]]);

        if (deltaR2 == 0.) continue;
        if (abs(distphi) <= M_PI) phinom.at(belongsto[i]) += particles[i].perp() * particles[i].phi() * pow(deltaR2, (beta-2)/2);
        else if ( distphi > M_PI) phinom.at(belongsto[i]) += particles[i].perp() * (-2 * M_PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);
        else if ( distphi < M_PI) phinom.at(belongsto[i]) += particles[i].perp() * (+2 * M_PI + particles[i].phi()) * pow(deltaR2, (beta-2)/2);
        ynom.at(belongsto[i]) += particles[i].perp() * particles[i].rap() * pow(deltaR2, (beta-2)/2);
        den.at(belongsto[i])  += particles[i].perp() * pow(deltaR2, (beta-2)/2);
      }

      // reset to new axes
      for (unsigned int j = 0; j < axes.size(); j++) {
        if (den[j] == 0.) axes.at(j) = axes[j];
        else {
          double phi_new = fmod( 2*M_PI + (phinom[j] / den[j]), 2*M_PI );
          double pt_new  = axes[j].perp();
          double y_new   = ynom[j] / den[j];
          double px = pt_new * cos(phi_new);
          double py = pt_new * sin(phi_new);
          double pz = pt_new * sinh(y_new);
          axes.at(j).reset(px, py, pz, axes[j].perp()/2);
        }
      }
    }


    void init() {

      /// Projections:

      FinalState fs(-4.5, 4.5, 0.*GeV);
      addProjection(fs, "FS");
      addProjection( FastJets(fs, FastJets::ANTIKT, 1.0), "AKT");
      addProjection( FastJets(fs, FastJets::CAM, 1.2)   , "CA" );

      /// Histograms:
      _h_camass.addHistogram(200, 300, bookHistogram1D(1, 1, 1));
      _h_camass.addHistogram(300, 400, bookHistogram1D(2, 1, 1));
      _h_camass.addHistogram(400, 500, bookHistogram1D(3, 1, 1));
      _h_camass.addHistogram(500, 600, bookHistogram1D(4, 1, 1));

      _h_filtmass.addHistogram(200, 300, bookHistogram1D(5, 1, 1));
      _h_filtmass.addHistogram(300, 400, bookHistogram1D(6, 1, 1));
      _h_filtmass.addHistogram(400, 500, bookHistogram1D(7, 1, 1));
      _h_filtmass.addHistogram(500, 600, bookHistogram1D(8, 1, 1));

      _h_ktmass.addHistogram(200, 300, bookHistogram1D( 9, 1, 1));
      _h_ktmass.addHistogram(300, 400, bookHistogram1D(10, 1, 1));
      _h_ktmass.addHistogram(400, 500, bookHistogram1D(11, 1, 1));
      _h_ktmass.addHistogram(500, 600, bookHistogram1D(12, 1, 1));

      _h_ktd12.addHistogram(200, 300, bookHistogram1D(13, 1, 1));
      _h_ktd12.addHistogram(300, 400, bookHistogram1D(14, 1, 1));
      _h_ktd12.addHistogram(400, 500, bookHistogram1D(15, 1, 1));
      _h_ktd12.addHistogram(500, 600, bookHistogram1D(16, 1, 1));

      _h_ktd23.addHistogram(200, 300, bookHistogram1D(17, 1 ,1));
      _h_ktd23.addHistogram(300, 400, bookHistogram1D(18, 1 ,1));
      _h_ktd23.addHistogram(400, 500, bookHistogram1D(19, 1 ,1));
      _h_ktd23.addHistogram(500, 600, bookHistogram1D(20, 1 ,1));

      _h_cat21.addHistogram(200, 300, bookHistogram1D(21, 1, 1));
      _h_cat21.addHistogram(300, 400, bookHistogram1D(22, 1, 1));
      _h_cat21.addHistogram(400, 500, bookHistogram1D(23, 1, 1));
      _h_cat21.addHistogram(500, 600, bookHistogram1D(24, 1, 1));

      _h_cat32.addHistogram(200, 300, bookHistogram1D(25, 1, 1));
      _h_cat32.addHistogram(300, 400, bookHistogram1D(26, 1, 1));
      _h_cat32.addHistogram(400, 500, bookHistogram1D(27, 1, 1));
      _h_cat32.addHistogram(500, 600, bookHistogram1D(28, 1, 1));

      _h_ktt21.addHistogram(200, 300, bookHistogram1D(29, 1, 1));
      _h_ktt21.addHistogram(300, 400, bookHistogram1D(30, 1, 1));
      _h_ktt21.addHistogram(400, 500, bookHistogram1D(31, 1, 1));
      _h_ktt21.addHistogram(500, 600, bookHistogram1D(32, 1, 1));

      _h_ktt32.addHistogram(200, 300, bookHistogram1D(33, 1, 1));
      _h_ktt32.addHistogram(300, 400, bookHistogram1D(34, 1, 1));
      _h_ktt32.addHistogram(400, 500, bookHistogram1D(35, 1, 1));
      _h_ktt32.addHistogram(500, 600, bookHistogram1D(36, 1, 1));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      using namespace fastjet;

      // Get anti-kt jets with p_T > 200 GeV, check abs(y) < 2
      // && fill in mass histograms
      const FastJets& ktfj = applyProjection<FastJets>(event, "AKT");
      PseudoJets ktjets = ktfj.pseudoJetsByPt(200*GeV);
      foreach (const PseudoJet ajet, ktjets) {
        if (abs(ajet.rap()) < 2) {
          _h_ktmass.fill(ajet.perp(), ajet.m(), weight);
        }
      }

      // Same as above but C/A jets
      const FastJets& cafj = applyProjection<FastJets>(event, "CA");
      PseudoJets cajets = cafj.pseudoJetsByPt(200*GeV);
      foreach (const PseudoJet ajet, cajets) {
        if (abs(ajet.rap()) < 2) {
          _h_camass.fill(ajet.perp(), ajet.m(), weight);
        }
      }

      // Split and filter.
      // Only do this to C/A jets in this analysis.
      foreach (const PseudoJet pjet, cajets) {
        if ( pjet.perp() > 600 || abs(pjet.rap()) > 2) continue;
        double dR = 0;
        bool unclustered = false;
        PseudoJets split_jets = splitjet(pjet, dR, cafj, unclustered);
        if ( (dR < 0.15) || (unclustered == false) ) continue;
        PseudoJet filt_jet = filterjet(split_jets, dR, 0.3);
        _h_filtmass.fill(filt_jet.perp(), filt_jet.m(), weight);
      }

      // Use the two last stages of clustering to get sqrt(d_12) and sqrt(d_23).
      // Only use anti-kt jets in this analysis.
      foreach (const PseudoJet pjet, ktjets) {
        if (pjet.perp() > 600 || abs(pjet.rap()) > 2) continue;
        ClusterSequence subjet_cseq(ktfj.clusterSeq()->constituents(pjet), JetDefinition(kt_algorithm, M_PI/2.));
        double d_12 = subjet_cseq.exclusive_dmerge(1) * M_PI*M_PI/4.;
        double d_23 = subjet_cseq.exclusive_dmerge(2) * M_PI*M_PI/4.;
        _h_ktd12.fill(pjet.perp(), sqrt(d_12), weight);
        _h_ktd23.fill(pjet.perp(), sqrt(d_23), weight);
      }

      // N-subjettiness, use beta = 1 (no rationale given).
      // Uses the functions defined above.
      // C/A jets first, anti-kt after.
      double beta = 1.;
      //Rcut is used for particles that are very far from the closest axis. At 10
      //is has no impact on the outcome of the calculation
      double Rcut = 10.;
      foreach (const PseudoJet pjet, cajets) {

        if (pjet.perp() > 600 || abs(pjet.rap()) > 2) continue;

        PseudoJets constituents = cafj.clusterSeq()->constituents(pjet);
        if (constituents.size() < 3) continue;

        PseudoJets axis1 = jetGetAxes(1, constituents, M_PI/2.0);
        PseudoJets axis2 = jetGetAxes(2, constituents, M_PI/2.0);
        PseudoJets axis3 = jetGetAxes(3, constituents, M_PI/2.0);

        double radius = 1.2;
        double tau1 = jetTauValue(beta, radius, constituents, axis1, Rcut);
        double tau2 = jetTauValue(beta, radius, constituents, axis2, Rcut);
        double tau3 = jetTauValue(beta, radius, constituents, axis3, Rcut);

        if (tau1 == 0 || tau2 == 0) continue;

        _h_cat21.fill(pjet.perp(), tau2/tau1, weight);
        _h_cat32.fill(pjet.perp(), tau3/tau2, weight);
      }

      foreach (const PseudoJet pjet, ktjets) {

        if (pjet.perp() > 600 || abs(pjet.rap()) > 2) continue;

        PseudoJets constituents = ktfj.clusterSeq()->constituents(pjet);
        if (constituents.size() < 3) continue;

        PseudoJets axis1 = jetGetAxes(1, constituents, M_PI/2.0);
        PseudoJets axis2 = jetGetAxes(2, constituents, M_PI/2.0);
        PseudoJets axis3 = jetGetAxes(3, constituents, M_PI/2.0);
        double radius = 1.0;
        double tau1 = jetTauValue(beta, radius, constituents, axis1, Rcut);
        double tau2 = jetTauValue(beta, radius, constituents, axis2, Rcut);
        double tau3 = jetTauValue(beta, radius, constituents, axis3, Rcut);
        if (tau1 == 0 || tau2 == 0) continue;

        _h_ktt21.fill(pjet.perp(), tau2/tau1, weight);
        _h_ktt32.fill(pjet.perp(), tau3/tau2, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      foreach (AIDA::IHistogram1D* hist, _h_camass.getHistograms()) {
        normalize(hist);
      }
      foreach (AIDA::IHistogram1D* hist, _h_filtmass.getHistograms()) {
        normalize(hist);
      }
      foreach (AIDA::IHistogram1D* hist, _h_ktmass.getHistograms()) {
        normalize(hist);
      }
      foreach (AIDA::IHistogram1D* hist, _h_ktd12.getHistograms()) {
        normalize(hist);
      }
      foreach (AIDA::IHistogram1D* hist, _h_ktd23.getHistograms()) {
        normalize(hist);
      }
      foreach (AIDA::IHistogram1D* hist, _h_cat21.getHistograms()) {
        normalize(hist);
      }
      foreach (AIDA::IHistogram1D* hist, _h_cat32.getHistograms()) {
        normalize(hist);
      }
      foreach (AIDA::IHistogram1D* hist, _h_ktt21.getHistograms()) {
        normalize(hist);
      }
      foreach (AIDA::IHistogram1D* hist, _h_ktt32.getHistograms()) {
        normalize(hist);
      }

    }

  private:

    BinnedHistogram<double> _h_camass;
    BinnedHistogram<double> _h_filtmass;
    BinnedHistogram<double> _h_ktmass;
    BinnedHistogram<double> _h_ktd12;
    BinnedHistogram<double> _h_ktd23;
    BinnedHistogram<double> _h_cat21;
    BinnedHistogram<double> _h_cat32;
    BinnedHistogram<double> _h_ktt21;
    BinnedHistogram<double> _h_ktt32;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1094564);

}
