// -*- C++ -*-
// "Acosta" underlying event analysis at CDF, inc. "Swiss Cheese"

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2004_S5839831.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  // Book histograms
  void CDF_2004_S5839831::init() {
    _pt90MaxAvg1800     = bookProfile1D(   1, 1, 1, "Average $p_T^\\text{max}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 1800 GeV");
    _pt90MinAvg1800     = bookProfile1D(   1, 1, 2, "Average $p_T^\\text{min}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 1800 GeV");
    _pt90Max1800        = bookProfile1D(   2, 1, 1, "$p_T^\\text{max}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 1800 GeV");
    _pt90Min1800        = bookProfile1D(   2, 1, 2, "$p_T^\\text{min}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 1800 GeV");
    _pt90Diff1800       = bookProfile1D(   2, 1, 3, "$p_T^\\text{diff}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 1800 GeV");
    _pt90Dbn1800Et40    = bookHistogram1D( 3, 1, 1, "$p_T$ distribution in MAX+MIN transverse cones for $40 < E_T^\\text{lead} < 80$ GeV at $\\sqrt{s}$ = 1800 GeV");
    _pt90Dbn1800Et80    = bookHistogram1D( 3, 1, 2, "$p_T$ distribution in MAX+MIN transverse cones for $80 < E_T^\\text{lead} < 120$ GeV at $\\sqrt{s}$ = 1800 GeV");
    _pt90Dbn1800Et120   = bookHistogram1D( 3, 1, 3, "$p_T$ distribution in MAX+MIN transverse cones for $120 < E_T^\\text{lead} < 160$ GeV at $\\sqrt{s}$ = 1800 GeV");
    _pt90Dbn1800Et160   = bookHistogram1D( 3, 1, 4, "$p_T$ distribution in MAX+MIN transverse cones for $160 < E_T^\\text{lead} < 200$ GeV at $\\sqrt{s}$ = 1800 GeV");
    _pt90Dbn1800Et200   = bookHistogram1D( 3, 1, 5, "$p_T$ distribution in MAX+MIN transverse cones for $200 < E_T^\\text{lead} < 270$ GeV at $\\sqrt{s}$ = 1800 GeV");
    _num90Max1800       = bookProfile1D(   4, 1, 1, "$N_\\text{max}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 1800 GeV");
    _num90Min1800       = bookProfile1D(   4, 1, 2, "$N_\\text{min}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 1800 GeV");
    _numTracksDbn1800MB = bookHistogram1D( 5, 1, 1, "Min bias track multiplicity distribution at $\\sqrt{s}$ = 1800 GeV");
    _ptDbn1800MB        = bookHistogram1D( 6, 1, 1, "Min bias $p_T$ distribution at $\\sqrt{s}$ = 1800 GeV");
    _pTSum1800_2Jet     = bookProfile1D(   7, 1, 1, "Swiss Cheese $p_T^\\text{sum}$ vs $E_T^\\text{lead}$ (for removal of 2 jets) at $\\sqrt{s}$ = 1800 GeV");
    _pTSum1800_3Jet     = bookProfile1D(   7, 1, 2, "Swiss Cheese $p_T^\\text{sum}$ vs $E_T^\\text{lead}$ (for removal of 3 jets) at $\\sqrt{s}$ = 1800 GeV");            
    _pt90Max630         = bookProfile1D(   8, 1, 1, "$p_T^\\text{max}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 630 GeV");
    _pt90Min630         = bookProfile1D(   8, 1, 2, "$p_T^\\text{min}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 630 GeV");
    _pt90Diff630        = bookProfile1D(   8, 1, 3, "$p_T^\\text{diff}$ vs $E_T^\\text{lead}$ at $\\sqrt{s}$ = 630 GeV");
    _pTSum630_2Jet      = bookProfile1D(   9, 1, 1, "Swiss Cheese $p_T^\\text{sum}$ vs $E_T^\\text{lead}$ (for removal of 2 jets) at $\\sqrt{s}$ = 630 GeV");
    _pTSum630_3Jet      = bookProfile1D(   9, 1, 2, "Swiss Cheese $p_T^\\text{sum}$ vs $E_T^\\text{lead}$ (for removal of 3 jets) at $\\sqrt{s}$ = 630 GeV");
    _numTracksDbn630MB  = bookHistogram1D(10, 1, 1, "Min bias track multiplicity distribution at $\\sqrt{s}$ = 630 GeV");
    _ptDbn630MB         = bookHistogram1D(11, 1, 1, "Min bias $p_T$ distribution at $\\sqrt{s}$ = 630 GeV");

    // Random number generators
    //_rngEtaMB = UniformRealRNG(RngBase(42u), UniformRealDist<>(-0.5, 0.5));
    //_rngPhiMB = UniformRealRNG(RngBase(47u), UniformRealDist<>(0, TWOPI));
  }



  // // Note that sorting is inverted, so that highest Et is at the front of the list
  // bool cmpJetsByEt(const Jet& a, const Jet& b) {
  //   return a.EtSum() > b.EtSum();
  // }
  // // Note that sorting is inverted, so that highest E is at the front of the list
  // bool cmpJetsByE(const Jet& a, const Jet& b) {
  //   return a.momentum().E() > b.momentum().E();
  // }



  const CDF_2004_S5839831::ConesInfo 
  CDF_2004_S5839831::calcTransCones(const double etaLead, const double phiLead, 
                                    const ParticleVector& tracks) {
    const double phiTransPlus = mapAngle0To2Pi(phiLead + PI/2.0);
    const double phiTransMinus = mapAngle0To2Pi(phiLead - PI/2.0);
    getLog() << Log::DEBUG << "phi_lead = " << phiLead 
             << " -> trans = (" << phiTransPlus 
             << ", " << phiTransMinus << ")" << endl;

    unsigned int numPlus(0), numMinus(0);
    double ptPlus(0), ptMinus(0);
    // Run over all charged tracks
    foreach (const Particle& t, tracks) {
      FourMomentum trackMom = t.momentum();
      const double pt = trackMom.pT();

      // Find if track mom is in either transverse cone
      if (deltaR(trackMom, etaLead, phiTransPlus) < 0.7) {
        ptPlus += pt;
        numPlus += 1;
      } else if (deltaR(trackMom, etaLead, phiTransMinus) < 0.7) {
        ptMinus += pt;
        numMinus += 1;
      }
    }

    ConesInfo rtn;
    // Assign N_{min,max} from N_{plus,minus}
    rtn.numMax = (ptPlus >= ptMinus) ? numPlus : numMinus;
    rtn.numMin = (ptPlus >= ptMinus) ? numMinus : numPlus;
    // Assign pT_{min,max} from pT_{plus,minus}
    rtn.ptMax = (ptPlus >= ptMinus) ? ptPlus : ptMinus;
    rtn.ptMin = (ptPlus >= ptMinus) ? ptMinus : ptPlus;
    rtn.ptDiff = fabs(rtn.ptMax - rtn.ptMin);

    getLog() << Log::DEBUG << "Min cone has " << rtn.numMin << " tracks -> " 
             << "pT_min = " << rtn.ptMin/GeV << " GeV" << endl;
    getLog() << Log::DEBUG << "Max cone has " << rtn.numMax << " tracks -> " 
             << "pT_max = " << rtn.ptMax/GeV << " GeV" << endl;

    return rtn;
  }



  const CDF_2004_S5839831::ConesInfo 
  CDF_2004_S5839831::calcTransCones(const FourMomentum& leadvec, 
                                    const ParticleVector& tracks) {
    const double etaLead = leadvec.pseudorapidity();
    const double phiLead = leadvec.azimuthalAngle();
    return calcTransCones(etaLead, phiLead, tracks);
  }



  vector<Jet> CDF_2004_S5839831::sortjets(vector<Jet>& jets) {
    sort(jets.begin(), jets.end(), cmpJetsByE); // was ...ByEt
    if ( getLog().isActive(Log::DEBUG) ) {
      getLog() << Log::DEBUG << "Jet Et sums = [" << endl;
      foreach (Jet& j, jets) {
        getLog() << Log::DEBUG << "  " << j.EtSum() << endl;
      }
      getLog() << Log::DEBUG << "]" << endl;
    }
    return jets;
  }



  // Do the analysis
  void CDF_2004_S5839831::analyze(const Event& event) {
    const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
    // Get the event weight
    const double weight = event.weight();


    {
      const ParticleVector tracks = applyProjection<FinalState>(event, "FS").particles();
      vector<Jet> jets = applyProjection<JetAlg>(event, "Jets").jets();
      if (!jets.empty()) {
        // Leading jet must be in central |eta| < 0.5 region
        sortjets(jets);
        const Jet leadingjet = jets.front();
        const double etaLead = leadingjet.momentum().pseudorapidity();
        if (fabs(etaLead) > 0.5) {
          getLog() << Log::DEBUG << "Leading jet eta = " << etaLead << " not in |eta| < 0.5" << endl;
        } else {

          // Get Et of the leading jet: used to bin histograms
          const double ETlead = leadingjet.EtSum();
          getLog() << Log::DEBUG << "Leading Et = " << ETlead/GeV << " GeV" << endl;

          // Multiplicity & pT distributions for sqrt(s) = 630 GeV, 1800 GeV
          const ConesInfo cones = calcTransCones(leadingjet.momentum(), tracks);
          if (fuzzyEquals(sqrtS/GeV, 630)) {
            _pt90Max630->fill(ETlead/GeV, cones.ptMax/GeV, weight);
            _pt90Min630->fill(ETlead/GeV, cones.ptMin/GeV, weight);
            _pt90Diff630->fill(ETlead/GeV, cones.ptDiff/GeV, weight);
          } else if (fuzzyEquals(sqrtS/GeV, 1800)) {
            _num90Max1800->fill(ETlead/GeV, cones.numMax, weight);
            _num90Min1800->fill(ETlead/GeV, cones.numMin, weight);
            _pt90Max1800->fill(ETlead/GeV, cones.ptMax/GeV, weight);
            _pt90Min1800->fill(ETlead/GeV, cones.ptMin/GeV, weight);
            _pt90Diff1800->fill(ETlead/GeV, cones.ptDiff/GeV, weight);
            _pt90MaxAvg1800->fill(ETlead/GeV, cones.ptMax/GeV, weight); // /numMax
            _pt90MinAvg1800->fill(ETlead/GeV, cones.ptMin/GeV, weight); // /numMin
            //
            const double ptTransTotal = cones.ptMax + cones.ptMin;
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

        }
      }
    }


    // Fill min bias total track multiplicity histos
    {
      const ParticleVector mbtracks = applyProjection<FinalState>(event, "MBFS").particles();
      if (fuzzyEquals(sqrtS/GeV, 1800)) {
        _numTracksDbn1800MB->fill(mbtracks.size(), weight);
      } else if (fuzzyEquals(sqrtS/GeV, 630)) {
        _numTracksDbn630MB->fill(mbtracks.size(), weight);
      }
      // Run over all charged tracks
      foreach (const Particle& t, mbtracks) {
        FourMomentum trackMom = t.momentum();
        const double pt = trackMom.pT();
        // Plot total pT distribution for min bias
        if (fuzzyEquals(sqrtS/GeV, 1800)) {
          _ptDbn1800MB->fill(pt/GeV, weight);
        } else if (fuzzyEquals(sqrtS/GeV, 630)) {
          _ptDbn630MB->fill(pt/GeV, weight);
        }
      }
    }



    // Construct "Swiss Cheese" pT distributions, with pT contributions from
    // tracks within R = 0.7 of the 1st, 2nd (and 3rd) jets being ignored. A
    // different set of charged tracks, with |eta| < 1.0, is used here, and all
    // the removed jets must have Et > 5 GeV.
    {
      const ParticleVector cheesetracks = applyProjection<FinalState>(event, "CheeseFS").particles();
      vector<Jet> cheesejets = applyProjection<JetAlg>(event, "Jets").jets();
      if (cheesejets.empty()) {
        getLog() << Log::DEBUG << "No 'cheese' jets found in event" << endl;
        return;
      }
      sortjets(cheesejets);
      if (cheesejets.size() > 1 &&
          fabs(cheesejets[0].momentum().pseudorapidity()) <= 0.5 &&
          cheesejets[0].momentum().Et()/GeV > 5.0 &&
          cheesejets[1].momentum().Et()/GeV > 5.0) {

        const double cheeseETlead = cheesejets[0].momentum().Et();

        const double eta1 = cheesejets[0].momentum().pseudorapidity();
        const double phi1 = cheesejets[0].momentum().azimuthalAngle();
        const double eta2 = cheesejets[1].momentum().pseudorapidity();
        const double phi2 = cheesejets[1].momentum().azimuthalAngle();

        double ptSumSub2(0), ptSumSub3(0);
        foreach (const Particle& t, cheesetracks) {
          FourMomentum trackMom = t.momentum();
          const double pt = trackMom.pT();

          // Subtracting 2 leading jets
          const double deltaR1 = deltaR(trackMom, eta1, phi1);
          const double deltaR2 = deltaR(trackMom, eta2, phi2);
          getLog() << Log::TRACE << "Track vs jet(1): "
                   << "|(" << trackMom.pseudorapidity() << ", " << trackMom.azimuthalAngle() << ") - "
                   << "|(" << eta1 << ", " << phi1 << ")| = " << deltaR1 << endl;
          getLog() << Log::TRACE << "Track vs jet(2): "
                   << "|(" << trackMom.pseudorapidity() << ", " << trackMom.azimuthalAngle() << ") - "
                   << "|(" << eta2 << ", " << phi2 << ")| = " << deltaR2 << endl;
          if (deltaR1 > 0.7 && deltaR2 > 0.7) {
            ptSumSub2 += pt;

            // Subtracting 3rd leading jet
            if (cheesejets.size() > 2 && 
                cheesejets[2].momentum().Et()/GeV > 5.0) {
              const double eta3 = cheesejets[2].momentum().pseudorapidity();
              const double phi3 = cheesejets[2].momentum().azimuthalAngle();
              const double deltaR3 = deltaR(trackMom, eta3, phi3);
              getLog() << Log::TRACE << "Track vs jet(3): "
                       << "|(" << trackMom.pseudorapidity() << ", " << trackMom.azimuthalAngle() << ") - "
                       << "|(" << eta3 << ", " << phi3 << ")| = " << deltaR3 << endl;
              if (deltaR3 > 0.7) {
                ptSumSub3 += pt;
              }
            }
          }
        }
      
        // Swiss Cheese sub 2,3 jets distributions for sqrt(s) = 630 GeV, 1800 GeV
        if (fuzzyEquals(sqrtS/GeV, 630)) {
          _pTSum630_2Jet->fill(cheeseETlead/GeV, ptSumSub2/GeV, weight);
          _pTSum630_3Jet->fill(cheeseETlead/GeV, ptSumSub3/GeV, weight);
        } else if (fuzzyEquals(sqrtS/GeV, 1800)) {
          _pTSum1800_2Jet->fill(cheeseETlead/GeV, ptSumSub2/GeV, weight);
          _pTSum1800_3Jet->fill(cheeseETlead/GeV, ptSumSub3/GeV, weight);
        }

      }
    }

    
  }
  
  
  void CDF_2004_S5839831::finalize() { 
    // Normalize to actual number of entries in pT dbn histos
    normalize(_pt90Dbn1800Et40,  1656.75);
    normalize(_pt90Dbn1800Et80,  4657.5);
    normalize(_pt90Dbn1800Et120, 5395.5);
    normalize(_pt90Dbn1800Et160, 7248.75);
    normalize(_pt90Dbn1800Et200, 2442.0);
    // and for min bias distributions:
    normalize(_numTracksDbn1800MB, 309718.25);
    normalize(_numTracksDbn630MB, 1101024.0);
    normalize(_ptDbn1800MB, 33600.0);
    normalize(_ptDbn630MB, 105088.0);
  }


}
