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

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////


// Book histograms
void CDF_2001_S4751469::init() {

  _dataToward = bookProfile1D(1, 1, 1, "pT sum toward total");
  _dataTrans = bookProfile1D(1, 1, 2, "pT sum transverse total");
  _dataAway = bookProfile1D(1, 1, 3, "pT sum away total");

}


// Do the analysis
void CDF_2001_S4751469::analyze(const Event& event) {
  Log log = getLog();

  // Analyse, with pT > 0.5 GeV AND |eta| < 1
  const TrackJet& tj = event.applyProjection(_trackjetproj);

  // Get jets, sorted by pT
  const TrackJet::Jets jets = tj.getJets();
  if (jets.size()==0) { return; }

  TrackJet::Jet leadingJet = jets[0];
  const double phiLead = leadingJet.getPtWeightedPhi();
  const double ptLead = leadingJet.getPtSum();

  // Cut on highest pT jet: combined 0.5 GeV < pT(lead) < 50 GeV
  if (ptLead < 0.5) return;
  if (ptLead > 50.0) return;
  //const size_t nBin = size_t(floor(ptLead));
  
  // Run over tracks in non-leading jets
  double ptSumToward(0.0), ptSumAway(0.0), ptSumTrans(0.0);
  for (TrackJet::Jets::const_iterator j = jets.begin()+1; j != jets.end(); ++j) {
    for (TrackJet::Jet::const_iterator p = j->begin(); p != j->end(); ++p) {
      // Calculate delta phi from leading jet
      double deltaPhi = fabs(p->azimuthalAngle() - phiLead);
      if (deltaPhi > PI) deltaPhi -= PI;
      assert(deltaPhi >= 0);
      assert(deltaPhi <= PI);

      // Get a pT sum value for each region (1 number for each region per event)
      if (deltaPhi < PI/3.0) {
        ptSumToward += pT(*p);
      } else if (deltaPhi < 2*PI/3.0) {
        ptSumTrans += pT(*p);
      } else {
        ptSumAway += pT(*p);
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

  // Update the proto-profile histograms
  _dataToward->fill(ptLead, ptSumToward, event.weight());
  _dataTrans->fill(ptLead, ptSumTrans, event.weight());
  _dataAway->fill(ptLead, ptSumAway, event.weight());
  //}

}


// Create the profile histograms
void CDF_2001_S4751469::finalize() {

}
