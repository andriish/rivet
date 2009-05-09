// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_LEADINGJETS.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  CDF_2008_LEADINGJETS::CDF_2008_LEADINGJETS()
  { 
    setBeams(PROTON, ANTIPROTON);
    
    // Final state for the jet finding
    const FinalState fsj(-4.0, 4.0, 0.0*GeV);
    addProjection(fsj, "FSJ");
    addProjection(FastJets(fsj, FastJets::CDFMIDPOINT, 0.7), "MidpointJets");
    
    // Final state for the sum(ET) distributions
    const FinalState fs(-1.0, 1.0, 0.0*GeV);
    addProjection(fs, "FS");
    
    // Charged final state for the distributions
    const ChargedFinalState cfs(-1.0, 1.0, 0.5*GeV);
    addProjection(cfs, "CFS");
  }


  // Book histograms
  void CDF_2008_LEADINGJETS::init() {
    _hist_pnchg      = bookProfile1D( 1, 1, 1, "Transverse region charged particle density",
        "$p_T(\\text{leading jet})$ / GeV", "$\\langle N_\\text{ch} \\rangle / \\text{d}\\eta\\,\\text{d}\\phi$");
    _hist_pmaxnchg   = bookProfile1D( 2, 1, 1, "TransMAX region charged particle density",
        "$p_T(\\text{leading jet})$ / GeV", "$\\langle N_\\text{ch} \\rangle / \\text{d}\\eta\\,\\text{d}\\phi$");
    _hist_pminnchg   = bookProfile1D( 3, 1, 1, "TransMIN region charged particle density",
        "$p_T(\\text{leading jet})$ / GeV", "$\\langle N_\\text{ch} \\rangle / \\text{d}\\eta\\,\\text{d}\\phi$");
    _hist_pdifnchg   = bookProfile1D( 4, 1, 1, "TransDIF region charged particle density",
        "$p_T(\\text{leading jet})$ / GeV", "$\\langle N_\\text{ch} \\rangle / \\text{d}\\eta\\,\\text{d}\\phi$");
    _hist_pcptsum    = bookProfile1D( 5, 1, 1, "Transverse region charged pT sum density",
        "$p_T(\\text{leading jet})$ / GeV", "$\\langle \\sum p_T^\\text{track} \\rangle / \\text{d}\\eta\\,\\text{d}\\phi$");
    _hist_pmaxcptsum = bookProfile1D( 6, 1, 1, "TransMAX region charged pT sum density",
        "$p_T(\\text{leading jet})$ / GeV", "$\\langle \\sum p_T^\\text{track} \\rangle / \\text{d}\\eta\\,\\text{d}\\phi$");
    _hist_pmincptsum = bookProfile1D( 7, 1, 1, "TransMIN region charged pT sum density",
        "$p_T(\\text{leading jet})$ / GeV", "$\\langle \\sum p_T^\\text{track} \\rangle / \\text{d}\\eta\\,\\text{d}\\phi$");
    _hist_pdifcptsum = bookProfile1D( 8, 1, 1, "TransDIF region charged pT sum density",
        "$p_T(\\text{leading jet})$ / GeV", "$\\langle \\sum p_T^\\text{track} \\rangle / \\text{d}\\eta\\,\\text{d}\\phi$");
    _hist_pcptave    = bookProfile1D( 9, 1, 1, "Transverse region charged pT average",
        "$p_T(\\text{leading jet})$ / GeV", "$\\langle p_T^\\text{track} \\rangle$");
    //_hist_onchg   = bookProfile1D( 1, 1, 1, "Overall number of charged particles");
    //_hist_ocptsum = bookProfile1D( 2, 1, 1, "Overall charged $p_\\perp$ sum");
    //_hist_oetsum  = bookProfile1D( 3, 1, 1, "Overall $E_\\perp$ sum");
  }


  // Do the analysis
  void CDF_2008_LEADINGJETS::analyze(const Event& e) {

    const FinalState& fsj = applyProjection<FinalState>(e, "FSJ");
    if (fsj.particles().size() < 1) {
      getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
      vetoEvent(e);
    }

    const FastJets& jetpro = applyProjection<FastJets>(e, "MidpointJets");
    const PseudoJets& jets = jetpro.pseudoJetsByPt();
    getLog() << Log::DEBUG << "Jet multiplicity = " << jets.size() << endl;

    // We require the leading jet to be within |eta|<2
    if (jets.size() < 1 || fabs(jets[0].eta()) >= 2) {
      getLog() << Log::DEBUG << "Failed jet cut" << endl;
      vetoEvent(e);
    }
    
    const double jetphi = jets[0].phi();
    const double jetpT  = jets[0].perp();
    getLog() << Log::DEBUG << "Leading jet: pT = " << jetpT
             << ", eta = " << jets[0].eta()
             << ", phi = " << jetphi << endl;

    // Get the event weight
    const double weight = e.weight();

    // Get the final states to work with for filling the distributions
    const FinalState& cfs = applyProjection<ChargedFinalState>(e, "CFS");

    size_t   numOverall(0),     numToward(0),     numTrans1(0),     numTrans2(0),     numAway(0)  ;
    double ptSumOverall(0.0), ptSumToward(0.0), ptSumTrans1(0.0), ptSumTrans2(0.0), ptSumAway(0.0);
    //double EtSumOverall(0.0), EtSumToward(0.0), EtSumTrans1(0.0), EtSumTrans2(0.0), EtSumAway(0.0);
    double ptMaxOverall(0.0), ptMaxToward(0.0), ptMaxTrans1(0.0), ptMaxTrans2(0.0), ptMaxAway(0.0);

    // Calculate all the charged stuff
    foreach (const Particle& p, cfs.particles()) {
      const double deltaPhi = delta_phi(p.momentum().azimuthalAngle(), jetphi);
      const double pT = p.momentum().pT();
      const double phi = p.momentum().azimuthalAngle();

      // Jets come with phi in [0 .. 2*Pi]. Particles come in [-Pi .. Pi].
      // Lovely, isn't it?
      /// @todo Use fastjet PseudoJet::phi_pi_pi (or whatever it's called)
      double rotatedphi = phi - jetphi;
      while (rotatedphi < 0) rotatedphi += 2*PI;

      ptSumOverall += pT;
      ++numOverall;
      if (pT > ptMaxOverall)
        ptMaxOverall = pT;

      if (deltaPhi < PI/3.0) {
        ptSumToward += pT;
        ++numToward;
        if (pT > ptMaxToward)
          ptMaxToward = pT;
      }
      else if (deltaPhi < 2*PI/3.0) {
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


#if 0   
    /// @todo Enable this part when we have the numbers from Rick Field

    // And now the same business for all particles (including neutrals)
    foreach (const Particle& p, fs.particles()) {
      const double deltaPhi = delta_phi(p.momentum().azimuthalAngle(), jetphi);
      const double ET = p.momentum().Et();
      const double phi = p.momentum().azimuthalAngle();

      // Jets come with phi in [0 .. 2*Pi]. Particles come in [-Pi .. Pi].
      // Lovely, isn't it?
      /// @todo Use FastJet methos to get the appropriate phi mapping
      double rotatedphi = phi - jetphi;
      while (rotatedphi < 0) rotatedphi += 2*PI;

      EtSumOverall += ET;

      if (deltaPhi < PI/3.0) {
        EtSumToward += ET;
      }
      else if (deltaPhi < 2*PI/3.0) {
        if (rotatedphi <= PI) {
          EtSumTrans1 += ET;
        }
        else {
          EtSumTrans2 += ET;
        }
      }
      else {
        EtSumAway += ET;
      }
    } // end all particle loop
#endif


    // Fill the histograms
    //_hist_tnchg->fill(jetpT, numToward/(4*PI/3), weight);
    _hist_pnchg->fill(jetpT, (numTrans1+numTrans2)/(4*PI/3), weight);
    _hist_pmaxnchg->fill(jetpT, (numTrans1>numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);
    _hist_pminnchg->fill(jetpT, (numTrans1<numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);
    _hist_pdifnchg->fill(jetpT, abs(numTrans1-numTrans2)/(2*PI/3), weight);
    //_hist_anchg->fill(jetpT, numAway/(4*PI/3), weight);

    //_hist_tcptsum->fill(jetpT, ptSumToward/(4*PI/3), weight);
    _hist_pcptsum->fill(jetpT, (ptSumTrans1+ptSumTrans2)/(4*PI/3), weight);
    _hist_pmaxcptsum->fill(jetpT, (ptSumTrans1>ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3), weight);
    _hist_pmincptsum->fill(jetpT, (ptSumTrans1<ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3), weight);
    _hist_pdifcptsum->fill(jetpT, fabs(ptSumTrans1-ptSumTrans2)/(2*PI/3), weight);
    //_hist_acptsum->fill(jetpT, ptSumAway/(4*PI/3), weight);

    //if (numToward > 0) {
    //  _hist_tcptave->fill(jetpT, ptSumToward/numToward, weight);
    //  _hist_tcptmax->fill(jetpT, ptMaxToward, weight);
    //}
    if ((numTrans1+numTrans2) > 0) {
      _hist_pcptave->fill(jetpT, (ptSumTrans1+ptSumTrans2)/(numTrans1+numTrans2), weight);
      //_hist_pcptmax->fill(jetpT, (ptMaxTrans1 > ptMaxTrans2 ? ptMaxTrans1 : ptMaxTrans2), weight);
    }
    //if (numAway > 0) {
    //  _hist_acptave->fill(jetpT, ptSumAway/numAway, weight);
    //  _hist_acptmax->fill(jetpT, ptMaxAway, weight);
    //}
  }


  void CDF_2008_LEADINGJETS::finalize() {  }


}
