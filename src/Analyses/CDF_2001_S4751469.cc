// -*- C++ -*-

// Field & Stuart underlying event analysis at CDF.
// Phys.Rev.D65:092002,2002 // no hep-ex code
// FNAL-PUB 01/211-E

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2001_S4751469.hh"
using namespace Rivet;

#include "Rivet/RivetAIDA.hh"
using namespace AIDA;


/////////////////////////////////////////////////


// Book histograms
void CDF_2001_S4751469::init() {
  _ptsumToward = bookProfile1D(1, 1, 1, "pT sum toward total");
  _ptsumTrans = bookProfile1D(1, 1, 2, "pT sum transverse total");
  _ptsumAway = bookProfile1D(1, 1, 3, "pT sum away total");
  /// @todo Get proper HepData entries for other plots in this paper
  _numToward = bookProfile1D("NToward", "Num toward", 50, 0.0, 50);
  _numTrans = bookProfile1D("NTrans", "Num transverse", 50, 0.0, 50);
  _numAway = bookProfile1D("NAway", "Num away", 50, 0.0, 50);
  // _numToward = bookProfile1D(2, 1, 1, "Num toward");
  // _numTrans = bookProfile1D(2, 1, 2, "Num transverse");
  // _numAway = bookProfile1D(2, 1, 3, "Num away");
}


// Do the analysis
void CDF_2001_S4751469::analyze(const Event& event) {
  Log log = getLog();

  // Analyse, with pT > 0.5 GeV AND |eta| < 1
  /// @todo Detector correction: randomly discard 8% of charged particles... how does this fit in with projections?
  const TrackJet& tj = event.applyProjection(_trackjetproj);

  // Get jets, sorted by pT
  const Jets jets = tj.getJets();
  if (jets.empty()) { return; }

  Jet leadingJet = jets[0];
  const double phiLead = leadingJet.getPtWeightedPhi();
  const double ptLead = leadingJet.getPtSum();

  // Cut on highest pT jet: combined 0.5 GeV < pT(lead) < 50 GeV
  if (ptLead/GeV < 0.5) return;
  if (ptLead/GeV > 50.0) return;

  // Run over tracks
  double ptSumToward(0.0), ptSumAway(0.0), ptSumTrans(0.0);
  size_t numToward(0), numTrans(0), numAway(0);
  for (Jets::const_iterator j = jets.begin(); j != jets.end(); ++j) {
    for (Jet::const_iterator p = j->begin(); p != j->end(); ++p) {
      // Calculate Delta(phi) from leading jet
      const double deltaPhi = delta_phi(p->azimuthalAngle(), phiLead);

      // Get pT sum and multiplicity values for each region 
      // (each is 1 number for each region per event)
      if (deltaPhi < PI/3.0) {
        ptSumToward += p->pT();
        ++numToward;
      } else if (deltaPhi < 2*PI/3.0) {
        ptSumTrans += p->pT();
        ++numTrans;
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
  _ptsumToward->fill(ptLead/GeV, ptSumToward/GeV, event.weight());
  _ptsumTrans->fill(ptLead/GeV, ptSumTrans/GeV, event.weight());
  _ptsumAway->fill(ptLead/GeV, ptSumAway/GeV, event.weight());

  // Log some event details
  log << Log::DEBUG 
      << "N [twd, away, trans] = ["
      << numToward << ", " 
      << numTrans << ", " 
      << numAway << "]" 
      << endl;

  // Update the N_jet profile histograms
  _numToward->fill(ptLead, numToward, event.weight());
  _numTrans->fill(ptLead, numTrans, event.weight());
  _numAway->fill(ptLead, numAway, event.weight());
}


void CDF_2001_S4751469::finalize() { }
