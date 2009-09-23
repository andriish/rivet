// -*- C++ -*-

// Field & Stuart underlying event analysis at CDF.
// Phys.Rev.D65:092002,2002 - no arXiv code.
// FNAL-PUB 01/211-E

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/LossyFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /* @brief "Field-Stuart" CDF Run I underlying event analysis
   * @author Andy Buckley
   * 
   * 
   * @par Run conditions
   * 
   * @arg \f$ \sqrt{s} = \f$ 1800 GeV
   * @arg Run with generic QCD events.
   * @arg Several \f$ p_\perp^\text{min} \f$ cutoffs are probably required to fill the profile histograms:
   *   @arg \f$ p_\perp^\text{min} = \f$ 0 (min bias), 10, 20 GeV
   * 
   */ 
  class CDF_2001_S4751469 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.5 \f$ GeV. Use a lossy charged FS projection, which
    /// randomly discards 8% of charged particles, as a kind of hacky detector 
    /// correction.
    CDF_2001_S4751469()
      : Analysis("CDF_2001_S4751469"),
        _totalNumTrans2(0), _totalNumTrans5(0), _totalNumTrans30(0),
        _sumWeightsPtLead2(0),_sumWeightsPtLead5(0), _sumWeightsPtLead30(0)
    {
      setBeams(PROTON, ANTIPROTON);
      const ChargedFinalState cfs(-1.0, 1.0, 0.5*GeV);
      const LossyFinalState lfs(cfs, 0.08);
      addProjection(lfs, "FS");
      addProjection(FastJets(lfs, FastJets::TRACKJET, 0.7), "TrackJet");
    }
    
    
    /// @name Analysis methods
    //@{
    
    // Book histograms
    void init() {
      // These histos bin N, pt in dphi
      _hist_num_dphi_2  = bookProfile1D("dummy_num_2", 50, 0, 180);
      _hist_num_dphi_5  = bookProfile1D("dummy_num_5", 50, 0, 180);
      _hist_num_dphi_30 = bookProfile1D("dummy_num_30", 50, 0, 180);
      _hist_pt_dphi_2   = bookProfile1D("dummy_pt_2", 50, 0, 180);
      _hist_pt_dphi_5   = bookProfile1D("dummy_pt_5", 50, 0, 180);
      _hist_pt_dphi_30  = bookProfile1D("dummy_pt_30", 50, 0, 180);
      
      _numvsDeltaPhi2 =  bookProfile1D(1, 1, 1);  
      _numvsDeltaPhi5 =  bookProfile1D(1, 1, 2);  
      _numvsDeltaPhi30 = bookProfile1D(1, 1, 3);  
      _pTvsDeltaPhi2 =   bookProfile1D(2, 1, 1);  
      _pTvsDeltaPhi5 =   bookProfile1D(2, 1, 2);  
      _pTvsDeltaPhi30 =  bookProfile1D(2, 1, 3);
      
      _numTowardMB = bookProfile1D(3, 1, 1);
      _numTransMB = bookProfile1D(3, 1, 2);
      _numAwayMB = bookProfile1D(3, 1, 3);
      _numTowardJ20 = bookProfile1D(4, 1, 1);
      _numTransJ20 = bookProfile1D(4, 1, 2);
      _numAwayJ20 = bookProfile1D(4, 1, 3);
      
      _ptsumTowardMB = bookProfile1D(5, 1, 1);
      _ptsumTransMB = bookProfile1D(5, 1, 2);
      _ptsumAwayMB = bookProfile1D(5, 1, 3);
      _ptsumTowardJ20 = bookProfile1D(6, 1, 1);
      _ptsumTransJ20 = bookProfile1D(6, 1, 2);
      _ptsumAwayJ20 = bookProfile1D(6, 1, 3);
      
      _ptTrans2 = bookHistogram1D(7, 1, 1);
      _ptTrans5 = bookHistogram1D(7, 1, 2);
      _ptTrans30 = bookHistogram1D(7, 1, 3);
    }
    

    /// Do the analysis
    void analyze(const Event& event) {
      
      // Analyse, with pT > 0.5 GeV AND |eta| < 1
      const JetAlg& tj = applyProjection<JetAlg>(event, "TrackJet");
      
      // Get jets, sorted by pT
      const Jets jets = tj.jetsByPt();
      if (jets.empty()) { 
        vetoEvent; 
      }

      Jet leadingJet = jets.front();
      const double phiLead = leadingJet.ptWeightedPhi();
      const double ptLead = leadingJet.ptSum();
      
      // Cut on highest pT jet: combined 0.5 GeV < pT(lead) < 50 GeV
      if (ptLead/GeV < 0.5) vetoEvent;
      if (ptLead/GeV > 50.0) vetoEvent;

      // Get the event weight
      const double weight = event.weight();

      // Count sum of all event weights in three pTlead regions
      if (ptLead/GeV > 2.0) {
        _sumWeightsPtLead2 += weight;
      }
      if (ptLead/GeV > 5.0) { 
        _sumWeightsPtLead5 += weight;
      }
      if (ptLead/GeV > 30.0) {
        _sumWeightsPtLead30 += weight;
      }
      
      // Run over tracks
      double ptSumToward(0.0), ptSumAway(0.0), ptSumTrans(0.0);
      size_t numToward(0), numTrans(0), numAway(0);
      // Reset the histos that bin N, pt in dphi
      _hist_num_dphi_2->reset();
      _hist_num_dphi_5->reset();
      _hist_num_dphi_30->reset();
      _hist_pt_dphi_2->reset();
      _hist_pt_dphi_5->reset();
      _hist_pt_dphi_30->reset();

      foreach (const Jet& j, jets) {
        foreach (const FourMomentum& p, j) {
          // Calculate Delta(phi) from leading jet
          const double dPhi = deltaPhi(p.azimuthalAngle(), phiLead);
          
          // Get pT sum and multiplicity values for each region 
          // (each is 1 number for each region per event)
          /// @todo Include event weight factor?
          if (dPhi < PI/3.0) {
            ptSumToward += p.pT();
            ++numToward;
            
          } else if (dPhi < 2*PI/3.0) {
            ptSumTrans += p.pT();
            ++numTrans;
            // Fill transverse pT distributions
            if (ptLead/GeV > 2.0) {
              _ptTrans2->fill(p.pT()/GeV, weight);
              _totalNumTrans2 += weight;
            }
            if (ptLead/GeV > 5.0) { 
              _ptTrans5->fill(p.pT()/GeV, weight);
              _totalNumTrans5 += weight;
            }
            if (ptLead/GeV > 30.0) {
              _ptTrans30->fill(p.pT()/GeV, weight);
              _totalNumTrans30 += weight;
            }
            
          } else {
            ptSumAway += p.pT();
            ++numAway;
          }
          
          // Fill histos to bin pt, N in dphi
          // dphi in degrees
          const double dPhideg = 360*dPhi/TWOPI;
          //
          if (ptLead/GeV > 2.0) {
            _hist_num_dphi_2->fill(dPhideg, 1, weight         );
            _hist_pt_dphi_2->fill (dPhideg, p.pT()/GeV, weight);

          }
          if (ptLead/GeV > 5.0) { 
            _hist_num_dphi_5->fill(dPhideg, 1, weight         );
            _hist_pt_dphi_5->fill (dPhideg, p.pT()/GeV, weight);

          }
          if (ptLead/GeV > 30.0) {
            _hist_num_dphi_30->fill(dPhideg, 1, weight         );
            _hist_pt_dphi_30->fill (dPhideg, p.pT()/GeV, weight);
          }





        }
      }
      
      // Log some event details
      getLog() << Log::DEBUG 
               << "pT [lead; twd, away, trans] = ["
               << ptLead << "; " 
               << ptSumToward << ", " 
               << ptSumAway << ", " 
               << ptSumTrans << "]" 
               << endl;
      
      // Update the pT profile histograms
      //
      // N, sumpt vs. dphi first
      // TODO: normalisation
      for (int i= 0; i < 50; i++) {
        // pT > 2 GeV  
        _numvsDeltaPhi2->fill(_hist_num_dphi_2->binMean(i), _hist_num_dphi_2->binHeight(i) * _hist_num_dphi_2->axis().binWidth(i));
        _pTvsDeltaPhi2->fill(_hist_pt_dphi_2->binMean(i), _hist_pt_dphi_2->binHeight(i)*_hist_pt_dphi_2->axis().binWidth(i));
        
        // pT > 5 GeV  
        _numvsDeltaPhi5->fill(_hist_num_dphi_5->binMean(i),_hist_num_dphi_5->binHeight(i)*_hist_num_dphi_5->axis().binWidth(i));
        _pTvsDeltaPhi5->fill(_hist_pt_dphi_5->binMean(i), _hist_pt_dphi_5->binHeight(i)*_hist_pt_dphi_5->axis().binWidth(i));
        
        // pT > 30 GeV  
        _numvsDeltaPhi30->fill(_hist_num_dphi_30->binMean(i),_hist_num_dphi_30->binHeight(i)*_hist_num_dphi_30->axis().binWidth(i));
        _pTvsDeltaPhi30->fill(_hist_pt_dphi_30->binMean(i), _hist_pt_dphi_30->binHeight(i)*_hist_pt_dphi_30->axis().binWidth(i));
      }
      
      _ptsumTowardMB->fill(ptLead/GeV, ptSumToward/GeV, weight);
      _ptsumTowardJ20->fill(ptLead/GeV, ptSumToward/GeV, weight);
      
      _ptsumTransMB->fill(ptLead/GeV, ptSumTrans/GeV, weight);
      _ptsumTransJ20->fill(ptLead/GeV, ptSumTrans/GeV, weight);
      
      _ptsumAwayMB->fill(ptLead/GeV, ptSumAway/GeV, weight);
      _ptsumAwayJ20->fill(ptLead/GeV, ptSumAway/GeV, weight);
      
      // Log some event details
      getLog() << Log::DEBUG 
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
    
    
    /// Normalize histos
    void finalize() {
      normalize(_ptTrans2, _totalNumTrans2 / _sumWeightsPtLead2);
      normalize(_ptTrans5, _totalNumTrans5 / _sumWeightsPtLead5);
      normalize(_ptTrans30, _totalNumTrans30 / _sumWeightsPtLead30);
    }
    
    //@}


  private:

    /// Sum total number of charged particles in the trans region, in 3 \f$ p_\perp^\text{lead} \f$ bins.
    double _totalNumTrans2, _totalNumTrans5, _totalNumTrans30;

    /// Sum the total number of events in 3 \f$ p_\perp^\text{lead} \f$ bins.
    double _sumWeightsPtLead2,_sumWeightsPtLead5, _sumWeightsPtLead30;


    /// @name Histogram collections
    //@{
    // These histos (binned in dphi) are filled per event and then reset
    // TODO: use LWH
    AIDA::IProfile1D *_hist_num_dphi_2, *_hist_num_dphi_5, *_hist_num_dphi_30;
    AIDA::IProfile1D *_hist_pt_dphi_2, *_hist_pt_dphi_5, *_hist_pt_dphi_30;

    // The sumpt vs. dphi and Nch vs. dphi histos
    AIDA::IProfile1D *_numvsDeltaPhi2, *_numvsDeltaPhi5, *_numvsDeltaPhi30;
    AIDA::IProfile1D *_pTvsDeltaPhi2, *_pTvsDeltaPhi5, *_pTvsDeltaPhi30;


    /// Profile histograms, binned in the \f$ p_T \f$ of the leading jet, for
    /// the \f$ p_T \f$ sum in the toward, transverse and away regions.
    AIDA::IProfile1D *_ptsumTowardMB,  *_ptsumTransMB,  *_ptsumAwayMB;
    AIDA::IProfile1D *_ptsumTowardJ20, *_ptsumTransJ20, *_ptsumAwayJ20;

    /// Profile histograms, binned in the \f$ p_T \f$ of the leading jet, for
    /// the number of charged particles per jet in the toward, transverse and
    /// away regions.
    AIDA::IProfile1D *_numTowardMB,  *_numTransMB,  *_numAwayMB;
    AIDA::IProfile1D *_numTowardJ20, *_numTransJ20, *_numAwayJ20;

    /// Histogram of \f$ p_T \f$ distribution for 3 different \f$ p_{T1} \f$ IR cutoffs.
    AIDA::IHistogram1D *_ptTrans2, *_ptTrans5, *_ptTrans30;
    //@}

  };

    

  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_2001_S4751469> plugin_CDF_2001_S4751469;

}
