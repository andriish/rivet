// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"

namespace Rivet {


  class UA1_1990_S2044935 : public Analysis {
  public:

    /// Constructor
    UA1_1990_S2044935() : Analysis("UA1_1990_S2044935") {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
      _sumwTrig = 0;
      _sumwTrig08 = 0;
      _sumwTrig40 = 0;
      _sumwTrig80 = 0;
    }
    

    /// @name Analysis methods
    //@{

    /// Book projections and histograms
    void init() { 
      addProjection(ChargedFinalState(-5.5, 5.5), "TriggerFS");
      addProjection(ChargedFinalState(-2.5, 2.5), "TrackFS");
      addProjection(Beam(), "Beam");
      const FinalState calofs(-6.0, 6.0);
      addProjection(TotalVisibleMomentum(calofs), "Mom");

      _hist_Nch200 = bookHistogram1D(1,1,1);
      _hist_Nch500 = bookHistogram1D(1,1,2);
      _hist_Nch900 = bookHistogram1D(1,1,3);
      _hist_Esigd3p200 = bookHistogram1D(2,1,1);
      _hist_Esigd3p500 = bookHistogram1D(2,1,2);
      _hist_Esigd3p900 = bookHistogram1D(2,1,3);
      _hist_Esigd3p08 = bookHistogram1D(3,1,1);
      _hist_Esigd3p40 = bookHistogram1D(4,1,1);
      _hist_Esigd3p80 = bookHistogram1D(5,1,1);
      _hist_Et200 = bookHistogram1D(9,1,1);
      _hist_Et500 = bookHistogram1D(10,1,1);
      _hist_Et900 = bookHistogram1D(11,1,1);
      _hist_Pt63 = bookProfile1D(8,1,1);
      _hist_Pt200 = bookProfile1D(6,1,1);
      _hist_Pt900 = bookProfile1D(7,1,1);
      _hist_Etavg200 = bookProfile1D(12,1,1); 
      _hist_Etavg500 = bookProfile1D(12,1,2); 
      _hist_Etavg900 = bookProfile1D(12,1,3);
    }
    

    void analyze(const Event& event) {
      // Trigger
      const FinalState& trigfs = applyProjection<FinalState>(event, "TriggerFS");
      unsigned int n_minus(0), n_plus(0);
      foreach (const Particle& p, trigfs.particles()) {
        const double eta = p.momentum().eta();
        if (inRange(eta, -5.5, -1.5)) n_minus++;
        else if (inRange(eta, 1.5, 5.5)) n_plus++;
      }
      getLog() << Log::DEBUG << "Trigger -: " << n_minus << ", Trigger +: " << n_plus << endl;
      if (n_plus == 0 || n_minus == 0) vetoEvent;      
      const double weight = event.weight();
      _sumwTrig += weight;

      // Use good central detector tracks
      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
      const FinalState& cfs = applyProjection<FinalState>(event, "TrackFS");
      const double Et = applyProjection<TotalVisibleMomentum>(event, "Mom").scalarET();
      const unsigned int nch = cfs.size();

      // Event level histos
      if (fuzzyEquals(sqrtS/GeV, 200, 1E-4)) {
        _hist_Nch200->fill(nch, weight);
        _hist_Et200->fill(Et/GeV, weight);
      } else if (fuzzyEquals(sqrtS/GeV, 500)) {
        _hist_Nch500->fill(nch, weight);
        _hist_Et500->fill(Et/GeV, weight);
      }	else if (fuzzyEquals(sqrtS/GeV, 900)) {
        _hist_Nch900->fill(nch, weight);
        _hist_Et900->fill(Et/GeV, weight);
      }

      // Particle/track level histos
      const double deta = 2 * 5.0;
      const double dphi = TWOPI;
      const double dnch_deta = nch/5.0; //< @todo No factor of 2 for both sides?
      foreach (const Particle& p, cfs.particles()) {
        const double pt = p.momentum().pT();
        const double scaled_weight = weight/(deta*dphi*pt);
        if (fuzzyEquals(sqrtS/GeV, 63, 1E-3)) {
          _hist_Pt63->fill(nch, pt/GeV, weight);
        } else if (fuzzyEquals(sqrtS/GeV, 200, 1E-4)) {
          _hist_Pt200->fill(nch, pt/GeV, weight);
          _hist_Etavg200->fill(nch, Et/GeV, weight);
          _hist_Esigd3p200->fill(pt/GeV, scaled_weight);
        } else if (fuzzyEquals(sqrtS/GeV, 500)) {
          _hist_Etavg500->fill(nch, Et/GeV, weight);
          _hist_Esigd3p500->fill(pt/GeV, scaled_weight);
        } else if (fuzzyEquals(sqrtS/GeV, 900)) {
          _hist_Pt900->fill(nch, p.momentum().pT()/GeV, weight);
          _hist_Etavg900->fill(nch, Et/GeV, weight);
          _hist_Esigd3p900->fill(pt/GeV, scaled_weight);
          // Also fill for specific dNch/deta ranges for 900 GeV
          if (inRange(dnch_deta, 0.8, 4)) {
            _sumwTrig08 += weight;
            _hist_Esigd3p08->fill(pt/GeV, scaled_weight);
          } else if (dnch_deta > 4 && dnch_deta <= 8) {
            _sumwTrig40 += weight;
            _hist_Esigd3p40->fill(pt/GeV, scaled_weight);
          } else if(dnch_deta > 8) {
            _sumwTrig80 += weight;
            _hist_Esigd3p80->fill(pt/GeV, scaled_weight);
          }
        } 
      }
      
    }
    
    
    void finalize() {
      const double xsec = crossSection();
      if (_sumwTrig > 0) {
        /// @todo Normalisation seems low by factor of ~2:
        normalize(_hist_Nch200, xsec/millibarn * _sumwTrig/sumOfWeights());
        normalize(_hist_Nch500, xsec/millibarn * _sumwTrig/sumOfWeights());
        normalize(_hist_Nch900, xsec/millibarn * _sumwTrig/sumOfWeights());
        //
        scale(_hist_Esigd3p200, xsec/millibarn * 1/_sumwTrig);
        scale(_hist_Esigd3p500, xsec/millibarn * 1/_sumwTrig);
        scale(_hist_Esigd3p900, xsec/millibarn * 1/_sumwTrig);
        //
        if (_sumwTrig08 > 0) scale(_hist_Esigd3p08, xsec/microbarn * 1/_sumwTrig08);
        if (_sumwTrig40 > 0) scale(_hist_Esigd3p40, xsec/microbarn * 1/_sumwTrig40);
        if (_sumwTrig80 > 0) scale(_hist_Esigd3p80, xsec/microbarn * 1/_sumwTrig80);
        //
        normalize(_hist_Et200, xsec/millibarn * _sumwTrig/sumOfWeights());
        normalize(_hist_Et500, xsec/millibarn * _sumwTrig/sumOfWeights());
        normalize(_hist_Et900, xsec/millibarn * _sumwTrig/sumOfWeights());
      }
    }
    
    //@}

    
  private:

    /// Weight counters
    double _sumwTrig, _sumwTrig08, _sumwTrig40, _sumwTrig80;
    
    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _hist_Nch200;
    AIDA::IHistogram1D* _hist_Nch500;
    AIDA::IHistogram1D* _hist_Nch900;
    AIDA::IHistogram1D* _hist_Esigd3p200;
    AIDA::IHistogram1D* _hist_Esigd3p500;
    AIDA::IHistogram1D* _hist_Esigd3p900;
    AIDA::IHistogram1D* _hist_Esigd3p08;
    AIDA::IHistogram1D* _hist_Esigd3p40;
    AIDA::IHistogram1D* _hist_Esigd3p80;
    AIDA::IProfile1D* _hist_Pt63;
    AIDA::IProfile1D* _hist_Pt200;
    AIDA::IProfile1D* _hist_Pt900;
    AIDA::IProfile1D* _hist_Etavg200;
    AIDA::IProfile1D* _hist_Etavg500;
    AIDA::IProfile1D* _hist_Etavg900;
    AIDA::IHistogram1D* _hist_Et200;
    AIDA::IHistogram1D* _hist_Et500;
    AIDA::IHistogram1D* _hist_Et900;
    //@}
    
  };
  
  
  
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<UA1_1990_S2044935> plugin_UA1_1990_S2044935;
  
}
