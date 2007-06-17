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
  for (size_t i = 0; i < _numBins; ++i) {
    _dataToward[i] = MiniHisto();
    _dataAway[i] = MiniHisto();
    _dataTrans[i] = MiniHisto();
  }
  _histToward = bookHistogram1D("PtSumToward", "pT sum toward total", _numBins, 0.0, 50.0);
  _histTrans = bookHistogram1D("PtSumTransverse", "pT sum transverse total", _numBins, 0.0, 50.0);
  _histAway = bookHistogram1D("PtSumAway", "pT sum away total", _numBins, 0.0, 50.0);

  // Data point sets
  _dpsToward = bookDataPointSet("PtSumTowardDPS", "pT sum toward total");
  _dpsTrans = bookDataPointSet("PtSumTransverseDPS", "pT sum transverse total");
  _dpsAway = bookDataPointSet("PtSumAwayDPS", "pT sum away total");
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
  Log& log = getLog();

  vector<double> valsToward, errsToward;
  vector<double> valsTrans, errsTrans;
  vector<double> valsAway, errsAway;
  vector<double> xvalsTw, xerrsTw, xvalsTr, xerrsTr, xvalsA, xerrsA;

  for (size_t bin = 0; bin < _numBins; ++bin) {
    const double binwidth = 50.0/_numBins;
    const double bincentre = (bin + 0.5) * binwidth;
    const double leadPt = double(bin) + 0.5;
    /// @todo Should really use proper profile histograms here.

    const double nToward = _dataToward[bin].numEntries;
    if (nToward) {
      const double avgPtToward = _dataToward[bin].sumPt/nToward;
      const double avgPt2Toward = _dataToward[bin].sumPtSq/nToward;
      _histToward->fill(leadPt, avgPtToward);
      xvalsTw.push_back(bincentre);
      xerrsTw.push_back(binwidth/2.0);
      _dpsToward->addPoint();
      valsToward.push_back(avgPtToward);
      errsToward.push_back(sqrt( avgPt2Toward - avgPtToward*avgPtToward ));
    }

    const double nTrans = _dataTrans[bin].numEntries;
    if (nTrans) {
      const double avgPtTrans = _dataTrans[bin].sumPt/nTrans;
      const double avgPt2Trans = _dataTrans[bin].sumPtSq/nTrans;
      _histTrans->fill(leadPt, avgPtTrans);
      xvalsTr.push_back(bincentre);
      xerrsTr.push_back(binwidth/2.0);
      _dpsTrans->addPoint();
      valsTrans.push_back(avgPtTrans);
      errsTrans.push_back(sqrt( avgPt2Trans - avgPtTrans*avgPtTrans ));
    }

    const double nAway = _dataAway[bin].numEntries;
    if (nAway) {
      const double avgPtAway = _dataAway[bin].sumPt/nAway;
      const double avgPt2Away = _dataAway[bin].sumPtSq/nAway;
      _histAway->fill(leadPt, avgPtAway);
      xvalsA.push_back(bincentre);
      xerrsA.push_back(binwidth/2.0);
      _dpsAway->addPoint();
      valsAway.push_back(avgPtAway);
      errsAway.push_back(sqrt( avgPt2Away - avgPtAway*avgPtAway ));
    }
  }

  log << Log::INFO << "Vals toward = " << valsToward << endl;
  log << Log::INFO << "Vals transverse = " << valsTrans << endl;
  log << Log::INFO << "Vals away = " << valsAway << endl;

  // Set DPS x-coordinate values and errors
  _dpsToward->setCoordinate(0, xvalsTw, xerrsTw);
  _dpsTrans->setCoordinate(0, xvalsTr, xerrsTr);
  _dpsAway->setCoordinate(0, xvalsA, xerrsA);

  // Set DPS y-coordinate values and errors
  const bool st1 = _dpsToward->setCoordinate(1, valsToward, errsToward);
  const bool st2 = _dpsTrans->setCoordinate(1, valsTrans, errsTrans);
  const bool st3 = _dpsAway->setCoordinate(1, valsAway, errsAway);
  log << Log::INFO << "DPS statuses = " << st1 << ", "<< st2 << ", " << st3 << endl;
}
