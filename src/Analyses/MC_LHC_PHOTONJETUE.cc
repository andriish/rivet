// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /* Underlying event in jet + isolated photon events
   * @author Andy Buckley
   */ 
  class MC_LHC_PHOTONJETUE : public Analysis {
  public:
    
    /// Constructor
    MC_LHC_PHOTONJETUE()
      : Analysis("MC_LHC_PHOTONJETUE")
    { 
      setBeams(PROTON, PROTON);
    }
    
    
    /// @name Analysis methods
    //@{
    
    // Book histograms and projections
    void init() {
      // Final state for the jet finding
      const FinalState fsj(-4.0, 4.0, 0.1*GeV);
      addProjection(fsj, "FSJ");
      addProjection(FastJets(fsj, FastJets::ANTIKT, 0.7), "Jets");
      
      // Charged final state for the distributions
      const ChargedFinalState cfs(-2.0, 2.0, 0.2*GeV);
      addProjection(cfs, "Tracks");

      // Photons (for isolation)
      const FinalState fsp(-2.0, 2.0, 20.0*GeV);
      IdentifiedFinalState photonfs(fsp);
      photonfs.acceptId(PHOTON);
      addProjection(photonfs, "Photons");

      // Histograms
      _hist_jetgamma_dR   = bookHistogram1D("gammajet-dR", 52, 0.0, 5.2);
      _hist_jetgamma_dphi = bookHistogram1D("gammajet-dphi", 50, 0.0, PI);
      const double MAXPT1 = 50.0;
      _hist_pnchg      = bookProfile1D("trans-nchg",     50, 0.0, MAXPT1);
      _hist_pmaxnchg   = bookProfile1D("trans-maxnchg",  50, 0.0, MAXPT1);
      _hist_pminnchg   = bookProfile1D("trans-minnchg",  50, 0.0, MAXPT1);
      _hist_pcptsum    = bookProfile1D("trans-ptsum",    50, 0.0, MAXPT1);
      _hist_pmaxcptsum = bookProfile1D("trans-maxptsum", 50, 0.0, MAXPT1);
      _hist_pmincptsum = bookProfile1D("trans-minptsum", 50, 0.0, MAXPT1);
      _hist_pcptave    = bookProfile1D("trans-ptavg",    50, 0.0, MAXPT1);
    }


    // Do the analysis
    void analyze(const Event& evt) {

      // Get jets
      const Jets jets = applyProjection<FastJets>(evt, "Jets").jetsByPt();
      getLog() << Log::DEBUG << "Jet multiplicity = " << jets.size() << endl;
      if (jets.size() < 1) {
        getLog() << Log::DEBUG << "No jets found" << endl;
        vetoEvent;
      }

      // Get leading jet properties
      const FourMomentum pjet = jets.front().momentum();
      const double jeteta = pjet.eta();
      const double jetphi = pjet.phi();
      const double jetpT  = pjet.pT();
      getLog() << Log::DEBUG << "Leading jet: pT = " << jetpT/GeV << " GeV"
               << ", eta = " << jeteta
               << ", phi = " << jetphi << endl;

      // Require the leading jet to be within |eta| < 2
      if (fabs(jeteta) > 2) {
        getLog() << Log::DEBUG << "Failed leading jet eta cut" << endl;
        vetoEvent;
      }

      // Get the leading photon
      const FinalState& photonfs = applyProjection<FinalState>(evt, "Photons");
      if (photonfs.size() < 1) {
        getLog() << Log::DEBUG << "No hard photons found" << endl;
        vetoEvent;
      }
      const FourMomentum pgamma = photonfs.particlesByPt().front().momentum();      

      // Check that leading photon is isolated from jets
      bool isolated = true;
      foreach (const Jet& j, jets) {
        if (deltaR(j.momentum(), pgamma) < 0.1) {
          isolated = false;
          break;
        }
      }
      if (!isolated) {
        getLog() << Log::DEBUG << "Leading photon is not isolated from jets" << endl;
        vetoEvent;
      }

      // Get leading photon properties
      const double gammaeta = pgamma.eta();
      const double gammaphi = pgamma.phi();
      const double gammapT  = pgamma.pT();
      getLog() << Log::DEBUG << "Leading photon: pT = " << gammapT/GeV << " GeV"
               << ", eta = " << gammaeta
               << ", phi = " << gammaphi << endl;

      // Get the event weight
      const double weight = evt.weight();

      // Fill jet1-photon separation histos
      _hist_jetgamma_dR->fill(deltaR(pgamma, pjet), weight);
      _hist_jetgamma_dphi->fill(deltaPhi(gammaphi, jetphi), weight);

      /// @todo Cut on back-to-backness of jet-photon?

      /// @todo Plot evolution of UE as a function of jet-photon angle
      /// @todo Plot evolution of UE as a function of photon pT

      // Get the final states to work with for filling the distributions
      const FinalState& cfs = applyProjection<ChargedFinalState>(evt, "Tracks");

      size_t   numOverall(0),     numToward(0),     numTrans1(0),     numTrans2(0),     numAway(0);
      double ptSumOverall(0.0), ptSumToward(0.0), ptSumTrans1(0.0), ptSumTrans2(0.0), ptSumAway(0.0);
      double ptMaxOverall(0.0), ptMaxToward(0.0), ptMaxTrans1(0.0), ptMaxTrans2(0.0), ptMaxAway(0.0);

      // Calculate all the charged stuff
      foreach (const Particle& p, cfs.particles()) {
        const double dPhi = deltaPhi(p.momentum().phi(), jetphi);
        const double pT = p.momentum().pT();
        const double phi = p.momentum().azimuthalAngle();
        const double rotatedphi = phi - jetphi;

        ptSumOverall += pT;
        ++numOverall;
        if (pT > ptMaxOverall) ptMaxOverall = pT;

        if (dPhi < PI/3.0) {
          ptSumToward += pT;
          ++numToward;
          if (pT > ptMaxToward) ptMaxToward = pT;
        }
        else if (dPhi < 2*PI/3.0) {
          if (rotatedphi <= PI) {
            ptSumTrans1 += pT;
            ++numTrans1;
            if (pT > ptMaxTrans1) {
              ptMaxTrans1 = pT;
            } else {
              ptSumTrans2 += pT;
              ++numTrans2;
              if (pT > ptMaxTrans2) ptMaxTrans2 = pT;
            }
          }
        }
        else {
          ptSumAway += pT;
          ++numAway;
          if (pT > ptMaxAway) ptMaxAway = pT;
        }
      }
      
      
      // Fill the histograms
      _hist_pnchg->fill(jetpT, (numTrans1+numTrans2)/(4*PI/3), weight);
      _hist_pmaxnchg->fill(jetpT, (numTrans1>numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);
      _hist_pminnchg->fill(jetpT, (numTrans1<numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);      
      _hist_pcptsum->fill(jetpT, (ptSumTrans1+ptSumTrans2)/(4*PI/3), weight);
      _hist_pmaxcptsum->fill(jetpT, (ptSumTrans1>ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3), weight);
      _hist_pmincptsum->fill(jetpT, (ptSumTrans1<ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3), weight);
      if ((numTrans1+numTrans2) > 0) {
        _hist_pcptave->fill(jetpT, (ptSumTrans1+ptSumTrans2)/(numTrans1+numTrans2), weight);
      }
    }
    
    
    void finalize() {  
      //
    }
    
    
  private:

    AIDA::IHistogram1D* _hist_jetgamma_dR;
    AIDA::IHistogram1D* _hist_jetgamma_dphi;
    
    AIDA::IProfile1D *_hist_pnchg;
    AIDA::IProfile1D *_hist_pmaxnchg;
    AIDA::IProfile1D *_hist_pminnchg;
    AIDA::IProfile1D *_hist_pcptsum;
    AIDA::IProfile1D *_hist_pmaxcptsum;
    AIDA::IProfile1D *_hist_pmincptsum;
    AIDA::IProfile1D *_hist_pcptave;  
    
  };
  
  
  
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<MC_LHC_PHOTONJETUE> plugin_MC_LHC_PHOTONJETUE;
  
}
