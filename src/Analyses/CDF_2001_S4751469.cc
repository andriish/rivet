// -*- C++ -*-

// Field & Stuart underlying event analysis at CDF.
// Phys.Rev.D65:092002,2002 - no arXiv code.
// FNAL-PUB 01/211-E

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2001_S4751469.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  // Book histograms
  void CDF_2001_S4751469::init() {
    _numTowardMB = bookProfile1D(3, 1, 1, "Num (toward) for min-bias");
    _numTransMB = bookProfile1D(3, 1, 2, "Num (transverse) for min-bias");
    _numAwayMB = bookProfile1D(3, 1, 3, "Num (away) for min-bias");
    _numTowardJ20 = bookProfile1D(4, 1, 1, "Num (toward) for JET20");
    _numTransJ20 = bookProfile1D(4, 1, 2, "Num (transverse) for JET20");
    _numAwayJ20 = bookProfile1D(4, 1, 3, "Num (away) for JET20");

    _ptsumTowardMB = bookProfile1D(5, 1, 1, "pT sum (toward) for min-bias");
    _ptsumTransMB = bookProfile1D(5, 1, 2, "pT sum (transverse) for min-bias");
    _ptsumAwayMB = bookProfile1D(5, 1, 3, "pT sum (away) for min-bias");
    _ptsumTowardJ20 = bookProfile1D(6, 1, 1, "pT sum (toward) for JET20");
    _ptsumTransJ20 = bookProfile1D(6, 1, 2, "pT sum (transverse) for JET20");
    _ptsumAwayJ20 = bookProfile1D(6, 1, 3, "pT sum (away) for JET20");

    _ptTrans2 = bookHistogram1D(7, 1, 1, "pT distribution (transverse, pT1 > 2 GeV)");
    _ptTrans5 = bookHistogram1D(7, 1, 2, "pT distribution (transverse, pT1 > 5 GeV)");
    _ptTrans30 = bookHistogram1D(7, 1, 3, "pT distribution (transverse, pT1 > 30 GeV)");
  }


  // Do the analysis
  void CDF_2001_S4751469::analyze(const Event& event) {
    Log log = getLog();

    // Analyse, with pT > 0.5 GeV AND |eta| < 1
    const TrackJet& tj = event.applyProjection(_trackjetproj);

    // Get jets, sorted by pT
    const Jets jets = tj.getJets();
    if (jets.empty()) { vetoEvent(event); }

    Jet leadingJet = jets[0];
    const double phiLead = leadingJet.getPtWeightedPhi();
    const double ptLead = leadingJet.getPtSum();

    // Cut on highest pT jet: combined 0.5 GeV < pT(lead) < 50 GeV
    if (ptLead/GeV < 0.5) vetoEvent(event);
    if (ptLead/GeV > 50.0) vetoEvent(event);

    // Get the event weight
    const double weight = event.weight();

    // Run over tracks
    double ptSumToward(0.0), ptSumAway(0.0), ptSumTrans(0.0);
    size_t numToward(0), numTrans(0), numAway(0);
    for (Jets::const_iterator j = jets.begin(); j != jets.end(); ++j) {
      for (Jet::const_iterator p = j->begin(); p != j->end(); ++p) {
        // Calculate Delta(phi) from leading jet
        const double deltaPhi = delta_phi(p->azimuthalAngle(), phiLead);
        
        // Get pT sum and multiplicity values for each region 
        // (each is 1 number for each region per event)
        /// @todo Include event weight factor?
        if (deltaPhi < PI/3.0) {
          ptSumToward += p->pT();
          ++numToward;
        } else if (deltaPhi < 2*PI/3.0) {
          ptSumTrans += p->pT();
          ++numTrans;
          
          // Fill transverse pT distributions
          _totalNumTrans += weight;
          if (ptLead/GeV > 2)  _ptTrans2->fill(p->pT()/GeV, weight);
          if (ptLead/GeV > 5)  _ptTrans5->fill(p->pT()/GeV, weight);
          if (ptLead/GeV > 30) _ptTrans30->fill(p->pT()/GeV, weight);

        } else {
          ptSumAway += p->pT();
          ++numAway;
        }
        
      }
    }
    
    // Log some event details
    log << Log::DEBUG 
        << "pT [lead; twd, away, trans] = ["
        << ptLead << "; " 
        << ptSumToward << ", " 
        << ptSumAway << ", " 
        << ptSumTrans << "]" 
        << endl;

    // Update the pT profile histograms
    _ptsumTowardMB->fill(ptLead/GeV, ptSumToward/GeV, weight);
    _ptsumTowardJ20->fill(ptLead/GeV, ptSumToward/GeV, weight);

    _ptsumTransMB->fill(ptLead/GeV, ptSumTrans/GeV, weight);
    _ptsumTransJ20->fill(ptLead/GeV, ptSumTrans/GeV, weight);

    _ptsumAwayMB->fill(ptLead/GeV, ptSumAway/GeV, weight);
    _ptsumAwayJ20->fill(ptLead/GeV, ptSumAway/GeV, weight);

    // Log some event details
    log << Log::DEBUG 
        << "N [twd, away, trans] = ["
        << numToward << ", " 
        << numTrans << ", " 
        << numAway << "]" 
        << endl;

    // Update the N_jet profile histograms
    _numTowardMB->fill(ptLead/GeV, numToward, weight);
    _numTowardJ20->fill(ptLead/GeV, numToward, weight);

    _numTransMB->fill(ptLead/GeV, numTrans, weight);
    _numTransJ20->fill(ptLead/GeV, numTrans, weight);

    _numAwayMB->fill(ptLead/GeV, numAway, weight);
    _numAwayJ20->fill(ptLead/GeV, numAway, weight);
  }


  void CDF_2001_S4751469::finalize() { 
    const double avgNumTrans = _totalNumTrans / sumOfWeights();
    normalize(_ptTrans2, avgNumTrans);
    normalize(_ptTrans5, avgNumTrans);
    normalize(_ptTrans30, avgNumTrans);
  }


}
