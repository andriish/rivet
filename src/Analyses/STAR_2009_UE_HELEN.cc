// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/SISConePlugin.hh"

namespace Rivet {


  /* STAR underlying event
   * @author Hendrik Hoeth
   */
  class STAR_2009_UE_HELEN : public Analysis {
  public:

    /// Constructor
    STAR_2009_UE_HELEN()
      : Analysis("STAR_2009_UE_HELEN")
    {
      setBeams(PROTON, PROTON);
    }


    /// @name Analysis methods
    //@{

    void init() {
      // Charged final state, |eta|<1, pT>0.2GeV
      const ChargedFinalState cfs(-1.0, 1.0, 0.2*GeV);
      addProjection(cfs, "CFS");

      // Neutral final state, |eta|<1, ET>0.2GeV (needed for the jets)
      const NeutralFinalState nfs(-1.0, 1.0, 0.2*GeV);
      addProjection(nfs, "NFS");

      // STAR can't see neutrons and K^0_L
      VetoedFinalState vfs(nfs);
      vfs.vetoNeutrinos();
      vfs.addVetoPairId(K0L);
      vfs.addVetoPairId(NEUTRON);
      addProjection(vfs, "VFS");

      // Jets are reconstructed from charged and neutral particles,
      // and the cuts are different (pT vs. ET), so we need to merge them.
      const MergedFinalState jfs(cfs, vfs);
      addProjection(jfs, "JFS");

      // SISCone, R = 0.7, overlap_threshold = 0.75
      addProjection(FastJets(jfs, FastJets::SISCONE, 0.7), "AllJets");

      // Book histograms
      _hist_pmaxnchg   = bookProfile1D( 1, 1, 1);
      _hist_pminnchg   = bookProfile1D( 2, 1, 1);
      _hist_anchg      = bookProfile1D( 3, 1, 1);
    }


    // Do the analysis
    void analyze(const Event& e) {
      const FinalState& cfs = applyProjection<ChargedFinalState>(e, "CFS");
      if (cfs.particles().size() < 1) {
        getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
        vetoEvent;
      }

      const Jets& alljets = applyProjection<FastJets>(e, "AllJets").jetsByPt();
      getLog() << Log::DEBUG << "Total jet multiplicity = " << alljets.size() << endl;

      Jets jets;
      foreach (const Jet jet, alljets) {
        if (jet.neutralEnergy() < 0.7)
          jets.push_back(jet);
      }

      // We require the leading jet to be within |eta|<(1-R)=0.3
      if (jets.size() < 1 || fabs(jets[0].momentum().eta()) >= 0.3) {
        getLog() << Log::DEBUG << "Failed jet cut" << endl;
        vetoEvent;
      }

      const double jetphi = jets[0].momentum().phi();
      const double jetpT  = jets[0].momentum().pT();

      // Get the event weight
      const double weight = e.weight();

      size_t numOverall(0),     numToward(0),     numTrans1(0),     numTrans2(0),     numAway(0)  ;
      double ptSumOverall(0.0), ptSumToward(0.0), ptSumTrans1(0.0), ptSumTrans2(0.0), ptSumAway(0.0);
      //double EtSumOverall(0.0), EtSumToward(0.0), EtSumTrans1(0.0), EtSumTrans2(0.0), EtSumAway(0.0);
      double ptMaxOverall(0.0), ptMaxToward(0.0), ptMaxTrans1(0.0), ptMaxTrans2(0.0), ptMaxAway(0.0);

      // Calculate all the charged stuff
      foreach (const Particle& p, cfs.particles()) {
        const double dPhi = deltaPhi(p.momentum().phi(), jetphi);
        const double pT = p.momentum().pT();
        const double phi = p.momentum().phi();
        double rotatedphi = phi - jetphi;
        while (rotatedphi < 0) rotatedphi += 2*PI;

        ptSumOverall += pT;
        ++numOverall;
        if (pT > ptMaxOverall) {
          ptMaxOverall = pT;
        }

        if (dPhi < PI/3.0) {
          ptSumToward += pT;
          ++numToward;
          if (pT > ptMaxToward)
            ptMaxToward = pT;
        }
        else if (dPhi < 2*PI/3.0) {
          if (rotatedphi <= PI) {
            ptSumTrans1 += pT;
            ++numTrans1;
            if (pT > ptMaxTrans1)
              ptMaxTrans1 = pT;
          }
          else {
            ptSumTrans2 += pT;
            ++numTrans2;
            if (pT > ptMaxTrans2)
              ptMaxTrans2 = pT;
          }
        }
        else {
          ptSumAway += pT;
          ++numAway;
          if (pT > ptMaxAway)
            ptMaxAway = pT;
        }
      } // end charged particle loop



      // Fill the histograms
      // @TODO: TAKE OUT THE 0.8 SCALE FACTOR WHEN WE GET THE FINAL DATA!!!!!!
      const double efficiency = 0.8;
      _hist_pmaxnchg->fill(jetpT, efficiency*(numTrans1>numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);
      _hist_pminnchg->fill(jetpT, efficiency*(numTrans1<numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);
      _hist_anchg->fill(jetpT, efficiency*numAway/(PI*0.7*0.7), weight); // jet area = pi*R^2

    }


    void finalize() {
      //
    }

    //@}


  private:

    AIDA::IProfile1D *_hist_pmaxnchg;
    AIDA::IProfile1D *_hist_pminnchg;
    AIDA::IProfile1D *_hist_anchg;

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<STAR_2009_UE_HELEN> plugin_STAR_2009_UE_HELEN;

}
