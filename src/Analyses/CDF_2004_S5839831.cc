// -*- C++ -*-
// "Acosta" underlying event analysis at CDF, inc. "Swiss Cheese"

#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2004_S5839831.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  CDF_2004_S5839831::CDF_2004_S5839831() {
    setBeams(PROTON, ANTIPROTON);
    addProjection(Beam(), "Beam");
    // NB. Charged track reconstruction efficiency has already been corrected in the data.
    const ChargedFinalState fs(-1.2, 1.2, 0.4*GeV);
    addProjection(fs, "FS");
    addProjection(FastJets(fs, FastJets::CDFJETCLU, 0.7), "Jets");
    // Restrict tracks to |eta| < 0.7 for the min bias part.
    const ChargedFinalState mbfs(-0.7, 0.7, 0.4*GeV);
    addProjection(mbfs, "MBFS");
    // Restrict tracks to |eta| < 1 for the Swiss-Cheese part.
    const ChargedFinalState cheesefs(-1.0, 1.0, 0.4*GeV);
    addProjection(cheesefs, "CheeseFS");
    addProjection(FastJets(cheesefs, FastJets::CDFJETCLU, 0.7), "CheeseJets");
  }



  // Book histograms
  void CDF_2004_S5839831::init() {
    getLog() << Log::WARN << "\n"
             << "***************************************************\n"
             << "This analysis is not considered reliable enough for\n"
             << "inclusion in MC tuning studies: be careful! Expert\n"
             << "help with ensuring that this analysis matches the\n"
             << "experiment's implementation would be very welcome!\n"
             << "***************************************************" 
             << endl;      

    const string ptmax = "$p_\\perp^\\text{max} \\rangle$";
    const string ptmin = "$p_\\perp^\\text{min}$";
    const string ptdiff = "$p_\\perp^\\text{diff}$";
    const string ptsum = "$p_\\perp^\\text{sum}$";
    const string ptmaxmean = "$\\langle p_\\perp^\\text{max} \\rangle$";
    const string ptminmean = "$\\langle p_\\perp^\\text{min} \\rangle$";
    const string et1 = "$E_\\perp^\\text{lead}$";
    string xlabel = et1 + " / GeV";
    string ylabel = "";

    _pt90MaxAvg1800 = 
      bookProfile1D(1, 1, 1, 
                    ptmaxmean + " vs. " + et1 + " at $\\sqrt{s}$ = 1800 GeV",
                    xlabel, ptmaxmean + " / GeV");
    _pt90MinAvg1800 = 
      bookProfile1D(1, 1, 2, 
                    ptminmean + " vs. " + et1 + " at $\\sqrt{s}$ = 1800 GeV",
                    xlabel, ptminmean + " / GeV"); 
    _pt90Max1800 = 
      bookProfile1D(2, 1, 1, 
                    ptmax + " vs. " + et1 + " at $\\sqrt{s}$ = 1800 GeV",
                    xlabel, ptmax + " / GeV");
    _pt90Min1800 = 
      bookProfile1D(2, 1, 2, 
                    ptmin + " vs. " + et1 + " at $\\sqrt{s}$ = 1800 GeV",
                    xlabel, ptmin + " / GeV");
    _pt90Diff1800 =
      bookProfile1D(2, 1, 3, 
                    ptdiff + " vs. " + et1 + " at $\\sqrt{s}$ = 1800 GeV",
                    xlabel, ptdiff + " / GeV");
    _num90Max1800 = 
      bookProfile1D(4, 1, 1, 
                    "$N_\\text{max}$ vs. " + et1 + " at $\\sqrt{s}$ = 1800 GeV",
                    xlabel, "$N_\\text{max}$");
    _num90Min1800 = 
      bookProfile1D(4, 1, 2, 
                    "$N_\\text{min}$ vs. " + et1 + " at $\\sqrt{s}$ = 1800 GeV",
                    xlabel, "$N_\\text{min}$");
    _pTSum1800_2Jet = 
      bookProfile1D(7, 1, 1, 
                    "Swiss Cheese " + ptsum + " vs. " + et1 + 
                    " (for removal of 2 jets) at $\\sqrt{s}$ = 1800 GeV",
                    xlabel, ptsum + " / GeV (2 jets removed)");
    _pTSum1800_3Jet = 
      bookProfile1D(7, 1, 2, 
                    "Swiss Cheese " + ptsum + " vs. " + et1 + 
                    " (for removal of 3 jets) at $\\sqrt{s}$ = 1800 GeV",
                    xlabel, ptsum + " / GeV (3 jets removed)");
    _pt90Max630 = 
      bookProfile1D(8, 1, 1,
                    ptmax + " vs. " + et1 + " at $\\sqrt{s}$ = 630 GeV",
                    xlabel, ptmax + " / GeV"); 
    _pt90Min630 = 
      bookProfile1D(8, 1, 2,
                    ptmin + " vs. " + et1 + " at $\\sqrt{s}$ = 630 GeV",
                    xlabel, ptmin + " / GeV"); 
    _pt90Diff630 =
      bookProfile1D(8, 1, 3,
                    ptdiff + " vs. " + et1 + " at $\\sqrt{s}$ = 630 GeV",
                    xlabel, ptdiff + " / GeV"); 
    _pTSum630_2Jet =
      bookProfile1D(9, 1, 1,
                    "Swiss Cheese " + ptsum + " vs. " + et1 + 
                    " (for removal of 2 jets) at $\\sqrt{s}$ = 630 GeV",
                    xlabel, ptsum + " / GeV (2 jets removed)");
    _pTSum630_3Jet =
      bookProfile1D(9, 1, 2,
                    "Swiss Cheese " + ptsum + " vs. " + et1 + 
                    " (for removal of 3 jets) at $\\sqrt{s}$ = 630 GeV",
                    xlabel, ptsum + " / GeV (3 jets removed)"); 
    

    string basetitle = "$p_\\perp$ distribution in MAX+MIN transverse cones for ";
    xlabel = "$p_\\perp / GeV";
    ylabel = "$\\mathrm{d}{\\sigma}/\\mathrm{d}{p_\\perp}$";
    /// @todo Check this normalisation defn (num-tracks x xsec?.)
    _pt90Dbn1800Et40 = 
      bookHistogram1D(3, 1, 1,
                      basetitle + "$40 < E_\\perp^\\text{lead} < 80$ GeV at $\\sqrt{s}$ = 1800 GeV",
                      xlabel, ylabel);
    /// @todo Check this normalisation defn (num-tracks x xsec?.)
    _pt90Dbn1800Et80 = 
      bookHistogram1D(3, 1, 2, 
                      basetitle + "$80 < E_\\perp^\\text{lead} < 120$ GeV at $\\sqrt{s}$ = 1800 GeV",
                      xlabel, ylabel);
    /// @todo Check this normalisation defn (num-tracks x xsec?.)
    _pt90Dbn1800Et120 = 
      bookHistogram1D(3, 1, 3, 
                      basetitle + "$120 < E_\\perp^\\text{lead} < 160$ GeV at $\\sqrt{s}$ = 1800 GeV",
                      xlabel, ylabel);
    /// @todo Check this normalisation defn (num-tracks x xsec?.)
    _pt90Dbn1800Et160 =
      bookHistogram1D(3, 1, 4, 
                      basetitle + "$160 < E_\\perp^\\text{lead} < 200$ GeV at $\\sqrt{s}$ = 1800 GeV",
                      xlabel, ylabel);
    /// @todo Check this normalisation defn (num-tracks x xsec?.)
    _pt90Dbn1800Et200 =
      bookHistogram1D(3, 1, 5, 
                      basetitle + "$200 < E_\\perp^\\text{lead} < 270$ GeV at $\\sqrt{s}$ = 1800 GeV",
                      xlabel, ylabel);
    /// @todo Check this normalisation defn (num-tracks x xsec?.)
    _ptDbn1800MB = 
      bookHistogram1D(6, 1, 1, 
                      "Min bias $p_\\perp$ distribution at $\\sqrt{s}$ = 1800 GeV",
                      xlabel, ylabel);


    xlabel = "$N_\\text{ch}$";
    ylabel = "$\\mathrm{d}{\\sigma}/\\mathrm{d}{N_\\text{ch}}$";
    /// @todo Check this normalisation defn.
    _numTracksDbn1800MB = 
      bookHistogram1D(5, 1, 1,
                      "Min bias track multiplicity distribution at $\\sqrt{s}$ = 1800 GeV",
                      xlabel, ylabel);
    /// @todo Check this normalisation defn.
    _numTracksDbn630MB = 
      bookHistogram1D(10, 1, 1, 
                      "Min bias track multiplicity distribution at $\\sqrt{s}$ = 630 GeV",
                      xlabel, ylabel);
    /// @todo Check this normalisation defn.
    _ptDbn630MB = 
      bookHistogram1D(11, 1, 1, 
                      "Min bias $p_\\perp$ distribution at $\\sqrt{s}$ = 630 GeV",
                      xlabel, ylabel);
  }



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
      getLog() << Log::DEBUG << "Running max/min analysis" << endl;
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
      getLog() << Log::DEBUG << "Running min bias multiplicity analysis" << endl;
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
      getLog() << Log::DEBUG << "Running Swiss Cheese analysis" << endl;
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
    /// @todo Check this normalisation defn.
    normalize(_pt90Dbn1800Et40,  1656.75);
    /// @todo Check this normalisation defn.
    normalize(_pt90Dbn1800Et80,  4657.5);
    /// @todo Check this normalisation defn.
    normalize(_pt90Dbn1800Et120, 5395.5);
    /// @todo Check this normalisation defn.
    normalize(_pt90Dbn1800Et160, 7248.75);
    /// @todo Check this normalisation defn.
    normalize(_pt90Dbn1800Et200, 2442.0);

    // and for min bias distributions:
    /// @todo Check this normalisation defn.
    normalize(_numTracksDbn1800MB, 309718.25);
    /// @todo Check this normalisation defn.
    normalize(_numTracksDbn630MB, 1101024.0);
    /// @todo Check this normalisation defn.
    normalize(_ptDbn1800MB, 33600.0);
    /// @todo Check this normalisation defn.
    normalize(_ptDbn630MB, 105088.0);
  }


}
