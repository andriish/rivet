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

      _hist_sigma200 = bookHistogram1D(1,1,1);
      _hist_sigma500 = bookHistogram1D(1,1,2);
      _hist_sigma900 = bookHistogram1D(1,1,3);
      _hist_Esigma200 = bookHistogram1D(2,1,1);
      _hist_Esigma500 = bookHistogram1D(2,1,2);
      _hist_Esigma900 = bookHistogram1D(2,1,3);
      _hist_Esigmapoint8 = bookHistogram1D(3,1,1);
      _hist_Esigma4 = bookHistogram1D(4,1,1);
      _hist_Esigma8 = bookHistogram1D(5,1,1);
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

      // Use good central detector tracks
      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
      const FinalState& cfs = applyProjection<FinalState>(event, "TrackFS");
      const double nch = cfs.size();
      if (fuzzyEquals(sqrtS/GeV, 200, 1E-4)) {
        _hist_sigma200->fill(nch, weight);
      } else if (fuzzyEquals(sqrtS/GeV, 500)) {
        _hist_sigma500->fill(nch, weight);
      }	else if (fuzzyEquals(sqrtS/GeV, 900)) {
        _hist_sigma900->fill(nch, weight);
      }
      foreach (const Particle& p, cfs.particles()) {
        /// @todo Check d3p weight factors
        const double pt = p.momentum().pT();
        const double scaled_weight = weight/(2*5*PI*pt);
        if (fuzzyEquals(sqrtS/GeV, 200, 1E-4)) {
          _hist_Esigma200->fill(pt/GeV, scaled_weight);
        }
        if (fuzzyEquals(sqrtS/GeV, 500)) {
          _hist_Esigma500->fill(pt/GeV, scaled_weight);
        }
        if (fuzzyEquals(sqrtS/GeV, 900)) {
          _hist_Esigma900->fill(pt/GeV, scaled_weight);
          // Also fill dNc/deta for 900 GeV
          const double dnch_deta = nch/5.0;
          if (inRange(dnch_deta, 0.8, 4)) {
            _hist_Esigmapoint8->fill(pt/GeV, scaled_weight);
          } else if (dnch_deta > 4 && dnch_deta <= 8) {
            _hist_Esigma4->fill(pt/GeV, scaled_weight);
          } else if(dnch_deta > 8) {
            _hist_Esigma8->fill(pt/GeV, scaled_weight);
          }
        }                
      }
      
      const double Et = applyProjection<TotalVisibleMomentum>(event, "Mom").scalarET();
      if (fuzzyEquals(sqrtS/GeV, 200, 1E-4)) {
        _hist_Et200->fill(Et/GeV, weight);
      } else if (fuzzyEquals(sqrtS/GeV, 500)) {
        _hist_Et500->fill(Et/GeV, weight);
      } else if (fuzzyEquals(sqrtS/GeV, 900)) {
        _hist_Et900->fill(Et/GeV, weight);
      }

      foreach (const Particle& p, cfs.particles()) {
        if (fuzzyEquals(sqrtS/GeV, 63, 1E-3)) {
          _hist_Pt63->fill(nch, p.momentum().pT()/GeV, weight);
        } else if (fuzzyEquals(sqrtS/GeV, 200, 1E-4)) {
          _hist_Pt200->fill(nch, p.momentum().pT()/GeV, weight);
          _hist_Etavg200->fill(nch, Et/GeV, weight);
        } else if (fuzzyEquals(sqrtS/GeV, 500)) {
          _hist_Etavg500->fill(nch, Et/GeV, weight);
        } else if (fuzzyEquals(sqrtS/GeV, 900)) {
          _hist_Pt900->fill(nch, p.momentum().pT()/GeV, weight);
          _hist_Etavg900->fill(nch, Et/GeV, weight);
        }
      }
    }
    
    
    void finalize() {
      ///@todo: get the total cross-sections from the generator
      ///@todo: check if the scaling for Esigmpoint8, Esigma4 and Esigma8 are correct.
      normalize(_hist_sigma200, 27.9);
      normalize(_hist_sigma500, 31.5);
      normalize(_hist_sigma900, 34.4);
      scale(_hist_Esigma200, 27.9/sumOfWeights());
      scale(_hist_Esigma500, 31.5/sumOfWeights());
      scale(_hist_Esigma900, 34.4/sumOfWeights());
      scale(_hist_Esigmapoint8, 34400./sumOfWeights());
      scale(_hist_Esigma4, 3440./sumOfWeights());
      scale(_hist_Esigma8, 344./sumOfWeights());
      normalize(_hist_Et200, 27.9);
      normalize(_hist_Et500, 31.5);
      normalize(_hist_Et900, 34.4);
    }
    
    //@}

    
  private:
    
    /// @name Histogram collections
    //@{
    AIDA::IHistogram1D* _hist_sigma200;
    AIDA::IHistogram1D* _hist_sigma500;
    AIDA::IHistogram1D* _hist_sigma900;
    AIDA::IHistogram1D* _hist_Esigma200;
    AIDA::IHistogram1D* _hist_Esigma500;
    AIDA::IHistogram1D* _hist_Esigma900;
    AIDA::IHistogram1D* _hist_Esigmapoint8;
    AIDA::IHistogram1D* _hist_Esigma4;
    AIDA::IHistogram1D* _hist_Esigma8;
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
