// -*- C++ -*-

// Field & Stuart underlying event analysis at CDF.
// Phys.Rev.D65:092002,2002 // no hep-ex code
// FNAL-PUB 01/211-E

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analysis/PRD65092002.hh"
using namespace Rivet;

#include "Rivet/RivetAIDA.hh"
using namespace AIDA;

#include "HepPDT/ParticleID.hh"
using namespace HepMC;


/////////////////////////////////////////////////


// Book histograms
void PRD65092002::init() {
  // Book mini histos (for the profile histo effect) and storage histos
  _dataToward.reserve(_numBins);
  _dataAway.reserve(_numBins);
  _dataTrans.reserve(_numBins);

  // Data point sets
  /// @todo Should really use proper profile histograms.
  _dpsToward = bookDataPointSet("PtSumToward", "pT sum toward total");
  _dpsTrans = bookDataPointSet("PtSumTransverse", "pT sum transverse total");
  _dpsAway = bookDataPointSet("PtSumAway", "pT sum away total");

  // Book mini-histos and initialise the DPS x-coords
  for (size_t bin = 0; bin < _numBins; ++bin) {
    _dataToward[bin] = MiniHisto();
    _dataAway[bin] = MiniHisto();
    _dataTrans[bin] = MiniHisto();

    const double binwidth = 50.0/_numBins;
    const double bincentre = (bin + 0.5) * binwidth;
    IMeasurement* meas;

    _dpsToward->addPoint();
    meas = _dpsToward->point(bin)->coordinate(0);
    meas->setValue(bincentre);
    meas->setErrorPlus(binwidth/2.0);
    meas->setErrorMinus(binwidth/2.0);

    _dpsTrans->addPoint();
    meas = _dpsTrans->point(bin)->coordinate(0);
    meas->setValue(bincentre);
    meas->setErrorPlus(binwidth/2.0);
    meas->setErrorMinus(binwidth/2.0);

    _dpsAway->addPoint();
    meas = _dpsAway->point(bin)->coordinate(0);
    meas->setValue(bincentre);
    meas->setErrorPlus(binwidth/2.0);
    meas->setErrorMinus(binwidth/2.0);
  }

}


// Do the analysis
void PRD65092002::analyze(const Event& event) {
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
  const size_t nBin = size_t(floor(ptLead));
  
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

  // Log some event details
  log << Log::DEBUG 
      << "pT [lead; twd, away, trans] = ["
      << ptLead << "; " 
      << ptSumToward << ", " 
      << ptSumAway << ", " 
      << ptSumTrans << "]" 
      << ", nbin = " << nBin
      << endl;

  // Update the proto-profile histograms
  if (nBin>_numBins) {
    log << Log::ERROR << "nBin out of range: " << nBin << endl;
  } else {
    _dataToward[nBin] += ptSumToward;
    _dataAway[nBin] += ptSumAway;
    _dataTrans[nBin] += ptSumTrans;
  }
}


// Create the profile histograms
void PRD65092002::finalize() {
  for (size_t bin = 0; bin < _numBins; ++bin) {

    const double nToward = _dataToward[bin].numEntries;
    if (nToward) {
      const double avgPt = _dataToward[bin].sumPt/nToward;
      const double avgPt2 = _dataToward[bin].sumPtSq/nToward;
      const double err = sqrt(avgPt2 - avgPt*avgPt);
      IDataPoint* pt = _dpsToward->point(bin);
      IMeasurement* meas = pt->coordinate(1);
      meas->setValue(avgPt);
      meas->setErrorPlus(err);
      meas->setErrorMinus(err);
    }

    const double nTrans = _dataTrans[bin].numEntries;
    if (nTrans) {
      const double avgPt = _dataTrans[bin].sumPt/nTrans;
      const double avgPt2 = _dataTrans[bin].sumPtSq/nTrans;
      const double err = sqrt(avgPt2 - avgPt*avgPt);
      IMeasurement* meas = _dpsTrans->point(bin)->coordinate(1);
      meas->setValue(avgPt);
      meas->setErrorPlus(err);
      meas->setErrorMinus(err);
    }

    const double nAway = _dataAway[bin].numEntries;
    if (nAway) {
      const double avgPt = _dataAway[bin].sumPt/nAway;
      const double avgPt2 = _dataAway[bin].sumPtSq/nAway;
      const double err = sqrt(avgPt2 - avgPt*avgPt);
      IMeasurement* meas = _dpsAway->point(bin)->coordinate(1);
      meas->setValue(avgPt);
      meas->setErrorPlus(err);
      meas->setErrorMinus(err);
    }

  }

}
