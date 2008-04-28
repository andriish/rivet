// -*- C++ -*-
// Underlying event analysis at CDF.

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2004_S5839831.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  // Book histograms
  void CDF_2004_S5839831::init() {
    _pt90Max1800 = bookProfile1D(2, 1, 1, "pTmax vs ET at sqrt{s} = 1800 GeV");
    _pt90Min1800 = bookProfile1D(2, 1, 2, "pTmin vs ET at sqrt{s} = 1800 GeV");
    _pt90Diff1800 = bookProfile1D(2, 1, 3, "pTdiff vs ET at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et40 = bookHistogram1D(3, 1, 1, "pT distribution in MAX+MIN transverse cones for 40 < ET < 80 GeV at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et80 = bookHistogram1D(3, 1, 2, "pT distribution in MAX+MIN transverse cones for 80 < ET < 120 GeV at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et120 = bookHistogram1D(3, 1, 3, "pT distribution in MAX+MIN transverse cones for 120 < ET < 160 GeV at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et160 = bookHistogram1D(3, 1, 4, "pT distribution in MAX+MIN transverse cones for 160 < ET < 200 GeV at sqrt{s} = 1800 GeV");
    _pt90Dbn1800Et200 = bookHistogram1D(3, 1, 5, "pT distribution in MAX+MIN transverse cones for 200 < ET < 270 GeV at sqrt{s} = 1800 GeV");
    _num90Max1800 = bookProfile1D(4, 1, 1, "Nmax vs ET at sqrt{s} = 1800 GeV");
    _num90Min1800 = bookProfile1D(4, 1, 2, "Nmin vs ET at sqrt{s} = 1800 GeV");    
    _numTracksDbn1800 = bookHistogram1D(5, 1, 1, "Track multiplicity distribution at sqrt{s} = 1800 GeV");
    _ptDbn1800 = bookHistogram1D(6, 1, 1, "pT distribution at sqrt{s} = 1800 GeV");
    _pTSum1800_2Jet = bookProfile1D(7, 1, 1, "pTsum vs ET (for removal of 2 jets) at sqrt{s} = 1800 GeV");
    _pTSum1800_3Jet = bookProfile1D(7, 1, 2, "pTsum vs ET (for removal of 3 jets) at sqrt{s} = 1800 GeV");            
    _pt90Max630 = bookProfile1D(8, 1, 1, "pTmax vs ET at sqrt{s} = 630 GeV");
    _pt90Min630 = bookProfile1D(8, 1, 2, "pTmin vs ET at sqrt{s} = 630 GeV");
    /// @todo Problem with HepData indexes?
    //_pt90Diff630 = bookProfile1D(8, 1, 3, "pTdiff vs ET at sqrt{s} = 630 GeV");
    _pTSum630_2Jet = bookProfile1D(9, 1, 1, "pTsum vs ET (for removal of 2 jets) at sqrt{s} = 630 GeV");
    _pTSum630_3Jet = bookProfile1D(9, 1, 2, "pTsum vs ET (for removal of 3 jets) at sqrt{s} = 630 GeV");
  }


  bool cmpJetsByEt(const Jet& a, const Jet& b) {
    /// @todo Use Et instead... definition?
    return a.getPtSum() < b.getPtSum();
  }


  // Do the analysis
  void CDF_2004_S5839831::analyze(const Event& event) {
    const ParticleVector tracks = applyProjection<FinalState>(event, "FS").particles();
    if (tracks.empty()) vetoEvent(event);
    /// @todo Make FastJets handle "no particles" events nicely
    const FastJets& jetproj = applyProjection<FastJets>(event, "Jets");
    vector<Jet> jets = jetproj.getJets();
    if (jets.empty()) vetoEvent(event);
    sort(jets.begin(), jets.end(), cmpJetsByEt);

    /// @todo Ensure there is only one well-defined primary vertex, i.e. no pileup. 
    /// Should be automatic as long as generator is run sensibly.

    /// @todo Need to implement track pT > 0.4 GeV/c and within 5.0/0.5cm of pp vertex (nominal?)

    // NB. Charged track reconstruction efficiency has already been corrected.

    // Leading jet must be in central |eta| < 0.5 region.
    const Jet leadingjet = jets.front();
    const double etaLead = leadingjet.vector().pseudorapidity();
    if (fabs(etaLead) > 0.5) vetoEvent(event);
    
    // Get Et of the leading jet: used to bin histograms.
    const double ETlead = leadingjet.getEtSum();

    // Get azimuthal angle of leading jet, to determine transverse regions
    const double phiLead = leadingjet.vector().azimuthalAngle();
    const double phiTransPlus = mapAngleMPiToPi(phiLead + PI/2.0);
    const double phiTransMinus = mapAngleMPiToPi(phiLead - PI/2.0);

    // Get the event weight.
    const double weight = event.weight();

    // Run over all charged tracks
    double ptMax(0), ptMin(0);

    for (ParticleVector::const_iterator t = tracks.begin(); t != tracks.end(); ++t) {
      double ptPlus(0), ptMinus(0);
      FourMomentum trackMom = t->getMomentum();
      const double pt = trackMom.pT();

      // Plot total pT distribution for min bias
      _ptDbn1800->fill(pt, weight);

      // Find if track mom is in either transverse cone
      if (deltaR(trackMom, etaLead, phiTransPlus) < 0.7) {
        ptPlus = pt;
      } else if (deltaR(trackMom, etaLead, phiTransMinus) < 0.7) {
        ptMinus = pt;
      }

      // Assign pT_{min,max} from pT_{plus,minus}
      ptMin = min(ptPlus, ptMinus);
      ptMax = max(ptPlus, ptMinus);

      // Increment total num particles in transverse regions
      if (ptMin || ptMax) {
        _totalNumTrans += weight;
      }
    }
    // Fill track multiplicity histo
    _numTracksDbn1800->fill(tracks.size(), weight);

    // AIDA::IProfile1D *_pt90Max1800,  *_pt90Min1800,  *_pt90Diff1800;
    // AIDA::IProfile1D *_pt90Max630,   *_pt90Min630,   *_pt90Diff630;
    // AIDA::IProfile1D *_num90Max1800,  *_num90Min1800;
    // AIDA::IProfile1D *_pTSum1800_2Jet, *_pTSum1800_3Jet;
    // AIDA::IProfile1D *_pTSum630_2Jet, *_pTSum630_3Jet;

    const double ptTransTotal = ptMax + ptMin;
    if (inRange(ETlead/GeV, 40, 80)) {
      _pt90Dbn1800Et40->fill(ptTransTotal/GeV, weight);
    } else if (inRange(ETlead/GeV, 80, 120)) {
      _pt90Dbn1800Et80->fill(ptTransTotal/GeV, weight);
    } else if (inRange(ETlead/GeV, 120, 160)) {
      _pt90Dbn1800Et120->fill(ptTransTotal/GeV, weight);
    } else if (inRange(ETlead/GeV, 160, 200)) {
      _pt90Dbn1800Et160->fill(ptTransTotal/GeV, weight);
    } else if (inRange(ETlead/GeV, 200, 270)) {
      _pt90Dbn1800Et200->fill(ptTransTotal/GeV, weight);
    }

  }


  void CDF_2004_S5839831::finalize() { 
    // const double avgNumTrans = _totalNumTrans / sumOfWeights();
    // normalize(_ptTrans2, avgNumTrans);
    // normalize(_ptTrans5, avgNumTrans);
    // normalize(_ptTrans30, avgNumTrans);
  }


}
