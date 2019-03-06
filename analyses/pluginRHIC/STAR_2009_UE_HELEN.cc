// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/SISConePlugin.hh"

namespace Rivet {


  /// @brief STAR underlying event
  /// @author Hendrik Hoeth
  class STAR_2009_UE_HELEN : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(STAR_2009_UE_HELEN);


    /// @name Analysis methods
    //@{

    void init() {

      const Cut c = Cuts::abseta < 1.0 && Cuts::pT > 0.2*GeV;

      // Charged and neutral final states
      const ChargedFinalState cfs(c);
      declare(cfs, "CFS");
      const NeutralFinalState nfs(c);
      declare(nfs, "NFS");

      // STAR can't see neutrons and K^0_L
      VetoedFinalState vfs(nfs);
      vfs.vetoNeutrinos();
      vfs.addVetoPairId(PID::K0L);
      vfs.addVetoPairId(PID::NEUTRON);
      declare(vfs, "VFS");

      // Jets are reconstructed from charged and neutral particles,
      // and the cuts are different (pT vs. ET), so we need to merge them.
      const MergedFinalState jfs(cfs, vfs);
      declare(jfs, "JFS");

      // SISCone, R = 0.7, overlap_threshold = 0.75
      declare(FastJets(jfs, FastJets::SISCONE, 0.7), "AllJets");

      // Book histograms
      book(_hist_pmaxnchg   , 1, 1, 1);
      book(_hist_pminnchg   , 2, 1, 1);
      book(_hist_anchg      , 3, 1, 1);
    }


    // Do the analysis
    void analyze(const Event& e) {
      const FinalState& cfs = apply<ChargedFinalState>(e, "CFS");
      if (cfs.particles().empty()) vetoEvent;

      const Jets& alljets = apply<FastJets>(e, "AllJets").jetsByPt();
      MSG_DEBUG("Total jet multiplicity = " << alljets.size());

      // The jet acceptance region is |eta|<(1-R)=0.3  (with R = jet radius)
      // Jets also must have a neutral energy fraction of < 0.7
      Jets jets;
      for (const Jet& jet : alljets) {
        if (jet.abseta() < 0.3 && jet.neutralEnergy()/jet.totalEnergy() < 0.7) jets.push_back(jet);
      }

      // This analysis requires a di-jet like event.
      // WARNING: There is more data in preparation, some of which does _not_ have this constraint!
      if (jets.size() != 2) vetoEvent;

      // The di-jet constraints in this analysis are:
      // - 2 and only 2 jets in the acceptance region
      // - delta(Phi) between the jets is > 150 degrees
      // - Pt_awayjet/Pt_towards_jet > 0.7
      if (deltaPhi(jets[0].phi(), jets[1].phi()) <= 5*PI/6 || jets[1].pT()/jets[0].pT() <= 0.7) vetoEvent;

      // Now lets start ...
      const double jetphi = jets[0].phi();
      const double jetpT  = jets[0].pT();


      // Calculate all the charged stuff
      size_t numTrans1(0), numTrans2(0), numAway(0);
      for (const Particle& p : cfs.particles()) {
        const double dPhi = deltaPhi(p.phi(), jetphi);
        const double pT = p.pT();
        const double phi = p.phi();
        double rotatedphi = phi - jetphi;
        while (rotatedphi < 0) rotatedphi += 2*PI;

        // WARNING: Hack to correct for the STAR tracking efficiency, in lieu of corrected data.
        if (rand01() > 0.87834 - exp(-1.48994-0.788432*pT)) continue;

        if (dPhi < PI/3.0) {
          // toward
        } else if (dPhi < 2*PI/3.0) {
          (rotatedphi <= PI ? numTrans1 : numTrans2) += 1;
        }
        else {
          numAway += 1;
        }
      }

      // Fill the histograms
      _hist_pmaxnchg->fill(jetpT, (numTrans1>numTrans2 ? numTrans1 : numTrans2)/(2*PI/3));
      _hist_pminnchg->fill(jetpT, (numTrans1<numTrans2 ? numTrans1 : numTrans2)/(2*PI/3));
      _hist_anchg->fill(jetpT, numAway/(PI*0.7*0.7)); // jet area = pi*R^2

    }

    //@}


    // Plots
    Profile1DPtr _hist_pmaxnchg, _hist_pminnchg, _hist_anchg;

  };


  DECLARE_RIVET_PLUGIN(STAR_2009_UE_HELEN);

}
