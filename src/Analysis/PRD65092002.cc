// -*- C++ -*-
// Field & Stuart underlying event analysis at CDF.
// Phys.Rev.D65:092002,2002 // no hep-ex code
// FNAL-PUB 01/211-E

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/PRD65092002.hh"
using namespace Rivet;

#include "AIDA/IHistogram1D.h"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////






// Book histograms
void PRD65092002::init() {
  // Book mini histos (for the profile histo effect) and storage histos
  const size_t numBins = 50;
  _dataToward.reserve(numBins);
  _dataAway.reserve(numBins);
  _dataTrans.reserve(numBins);
  _histToward = bookHistogram1D("PtSumToward", "pT sum toward total", numBins, 0.0, 50.0);
  _histTrans = bookHistogram1D("PtSumTransverse", "pT sum transverse total", numBins, 0.0, 50.0);
  _histAway = bookHistogram1D("PtSumAway", "pT sum away total", numBins, 0.0, 50.0);
}


// Do the analysis
void PRD65092002::analyze(const Event& event) {
  Log log = getLog();

  // Analyse, with pT > 0.5 GeV AND |eta| < 1
  const TrackJet& tj = event.applyProjection(_trackjetproj);

  // Get jets, sorted by pT
  const TrackJet::Jets jets = tj.getJets();
  TrackJet::Jet leadingJet = jets[0];
  const double phiLead = leadingJet.getPtWeightedPhi();
  const double ptLead = leadingJet.getPtSum();

  // Cut on highest pT jet: combined 0.5 GeV < pT(lead) < 50 GeV
  if (ptLead < 0.5) return;
  if (ptLead < 0.5) return;
  const size_t nBin = size_t(floor(ptLead/50.0));
  
  // Run over tracks in non-leading jets
  double ptSumToward(0.0), ptSumAway(0.0), ptSumTrans(0.0);
  for (TrackJet::Jets::const_iterator j = jets.begin()+1; j != jets.end(); ++j) {
    for (TrackJet::Jet::const_iterator p = j->begin(); p != j->end(); ++p) {
      // Calculate delta phi from leading jet
      double deltaPhi = fabs(p->phi() - phiLead);
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

  // Update the proto-profile histograms
  _dataToward[nBin] += ptSumToward;
  _dataAway[nBin] += ptSumAway;
  _dataTrans[nBin] += ptSumTrans;
}


// Create the profile histograms
void PRD65092002::finalize() { 
  for (size_t bin = 0; bin < 50; ++bin) {
    const double leadPt = double(bin) + 0.5;
    /// @todo Should really use proper profile histograms here.
    /// @todo Should also compute the error, using var = <pt^2> - <pt>^2
    _histToward->fill(leadPt, _dataToward[bin].sumPt/double(_dataToward[bin].numEntries));
    _histAway->fill(leadPt, _dataAway[bin].sumPt/double(_dataAway[bin].numEntries));
    _histTrans->fill(leadPt, _dataTrans[bin].sumPt/double(_dataTrans[bin].numEntries));
  }
}
