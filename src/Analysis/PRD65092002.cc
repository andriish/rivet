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
  _histToward = bookHistogram1D("PtSumToward", "pT sum toward total", 50, 0.0, 50.0);
  _histTrans = bookHistogram1D("PtSumTransverse", "pT sum transverse total", 50, 0.0, 50.0);
  _histAway = bookHistogram1D("PtSumAway", "pT sum away total", 50, 0.0, 50.0);
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

  // Cut on highest pT jet: combined 0.5 GeV < pT(lead) < 50 GeV
  if (leadingJet.getPtSum() < 0.5) return;
  if (leadingJet.getPtSum() < 0.5) return;
  
  // Run over tracks in non-leading jets
  for (TrackJet::Jets::const_iterator j = jets.begin()+1; j != jets.end(); ++j) {
    for (TrackJet::Jet::const_iterator p = j->begin(); p != j->end(); ++p) {
      // Calculate delta phi from leading jet
      const double deltaPhi = fabs(p->phi() - phiLead);
      assert(deltaPhi >= 0);
      assert(deltaPhi <= PI);
      /// @todo Is this really right? Shouldn't the phi values be calculated relative to a weighted vector?
      
      // Histogram the number of particles and the pT sum in this region
      if (deltaPhi < PI/3.0) {
        _histToward->fill(pT(*p));
      } else if (deltaPhi < 2*PI/3.0) {
        _histTrans->fill(pT(*p));
      } else {
        _histAway->fill(pT(*p));
      }

    }
  }
}


// Finalize
void PRD65092002::finalize() { }
