// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/MC_LHC_LEADINGJETS.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  MC_LHC_LEADINGJETS::MC_LHC_LEADINGJETS()
  { 
    setBeams(PROTON, PROTON);
    
    // Final state for the jet finding
    const FinalState fsj(-4.0, 4.0, 0.0*GeV);
    addProjection(fsj, "FSJ");
    addProjection(FastJets(fsj, FastJets::KT, 0.7), "Jets");
        
    // Charged final state for the distributions
    const ChargedFinalState cfs(-1.0, 1.0, 0.5*GeV);
    addProjection(cfs, "CFS");
  }
  

  
  // Book histograms
  void MC_LHC_LEADINGJETS::init() {
    const string xlabel = "Leading jet $p_\\perp$ / GeV";
    const string pdensityTeX = "$\\mathrm{d}{N_\\text{ch}}/\\mathrm{d}{\\phi}$";
    const string ptsumdensityTeX = "$\\mathrm{d}{\\sum p_\\perp^\\text{sum}}/\\mathrm{d}{\\phi}$";
    const string ptavgTeX = "$\\langle p_\\perp \\rangle$";
    const double maxpt1 = 500.0/GeV;
    _hist_pnchg      = bookProfile1D(histoPath("trans-nchg"),
                                     "Transverse region charged particle density", 
                                     xlabel, pdensityTeX, 50, 0.0, maxpt1);
    _hist_pmaxnchg   = bookProfile1D(histoPath("trans-maxnchg"),
                                     "TransMAX region charged particle density",
                                     xlabel, pdensityTeX, 50, 0.0, maxpt1);
    _hist_pminnchg   = bookProfile1D(histoPath("trans-minnchg"),
                                     "TransMIN region charged particle density",
                                     xlabel, pdensityTeX, 50, 0.0, maxpt1);
    _hist_pcptsum    = bookProfile1D(histoPath("trans-ptsum"),
                                     "Transverse region charged pT sum density",
                                     xlabel, ptsumdensityTeX, 50, 0.0, maxpt1);
    _hist_pmaxcptsum = bookProfile1D(histoPath("trans-maxptsum"),
                                     "TransMAX region charged pT sum density",
                                     xlabel, ptsumdensityTeX, 50, 0.0, maxpt1);
    _hist_pmincptsum = bookProfile1D(histoPath("trans-minptsum"),
                                     "TransMIN region charged pT sum density",
                                     xlabel, ptsumdensityTeX, 50, 0.0, maxpt1);
    _hist_pcptave    = bookProfile1D(histoPath("trans-ptavg"),
                                     "Transverse region charged pT average",
                                     xlabel, ptavgTeX, 50, 0.0, maxpt1);
  }


  // Do the analysis
  void MC_LHC_LEADINGJETS::analyze(const Event& e) {

    const FinalState& fsj = applyProjection<FinalState>(e, "FSJ");
    if (fsj.particles().empty()) {
      getLog() << Log::DEBUG << "Failed multiplicity cut" << endl;
      vetoEvent(e);
    }

    const FastJets& jetpro = applyProjection<FastJets>(e, "Jets");
    const Jets jets = jetpro.jetsByPt();
    getLog() << Log::DEBUG << "Jet multiplicity = " << jets.size() << endl;

    // We require the leading jet to be within |eta| < 2
    if (jets.size() < 1 || fabs(jets[0].momentum().pseudorapidity()) > 2) {
      getLog() << Log::DEBUG << "Failed jet cut" << endl;
      vetoEvent(e);
    }
    
    const double jetphi = jets[0].momentum().phi();
    const double jetpT  = jets[0].momentum().pT();
    getLog() << Log::DEBUG << "Leading jet: pT = " << jetpT
             << ", eta = " << jets[0].momentum().pseudorapidity()
             << ", phi = " << jetphi << endl;

    // Get the event weight
    const double weight = e.weight();

    // Get the final states to work with for filling the distributions
    const FinalState& cfs = applyProjection<ChargedFinalState>(e, "CFS");

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
    //_hist_tnchg->fill(jetpT, numToward/(4*PI/3), weight);
    _hist_pnchg->fill(jetpT, (numTrans1+numTrans2)/(4*PI/3), weight);
    _hist_pmaxnchg->fill(jetpT, (numTrans1>numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);
    _hist_pminnchg->fill(jetpT, (numTrans1<numTrans2 ? numTrans1 : numTrans2)/(2*PI/3), weight);
    //_hist_pdifnchg->fill(jetpT, abs(numTrans1-numTrans2)/(2*PI/3), weight);
    //_hist_anchg->fill(jetpT, numAway/(4*PI/3), weight);

    //_hist_tcptsum->fill(jetpT, ptSumToward/(4*PI/3), weight);
    _hist_pcptsum->fill(jetpT, (ptSumTrans1+ptSumTrans2)/(4*PI/3), weight);
    _hist_pmaxcptsum->fill(jetpT, (ptSumTrans1>ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3), weight);
    _hist_pmincptsum->fill(jetpT, (ptSumTrans1<ptSumTrans2 ? ptSumTrans1 : ptSumTrans2)/(2*PI/3), weight);
    //_hist_pdifcptsum->fill(jetpT, fabs(ptSumTrans1-ptSumTrans2)/(2*PI/3), weight);
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


  void MC_LHC_LEADINGJETS::finalize() {  
    //
  }


}
