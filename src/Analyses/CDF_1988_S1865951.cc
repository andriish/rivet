#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Analyses/CDF_1988_S1865951.hh"

namespace Rivet {

    /// Default constructor
    CDF_1988_S1865951::CDF_1988_S1865951() 
    : Analysis("CDF_1988_S1865951") {
      const ChargedFinalState cfs(-1.,1.);
      addProjection(cfs, "CFS");
      addProjection(TotalVisibleMomentum(cfs), "Mom");
      addProjection(Beam(), "Beam");
    }

    /// @name Analysis methods
    //@{
    void CDF_1988_S1865951::init() { 
      _hist_pt630 = bookHistogram1D(1, 1, 1);
      _hist_pt1800 = bookHistogram1D(2, 1, 1);
    }
    
    void CDF_1988_S1865951::analyze(const Event& event) {
      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
      const FinalState& fs = applyProjection<ChargedFinalState>(event, "CFS");
      const double weight = event.weight();

      foreach (Particle p, fs.particles())
      {
        double pt = p.momentum().pT();
        if (fuzzyEquals(sqrtS, 630/GeV))
        {
          _hist_pt630->fill(pt, weight/(10.*M_PI*pt));
        }
        else if (fuzzyEquals(sqrtS, 1800/GeV))
        {
          _hist_pt1800->fill(pt, weight/(10.*M_PI*pt));
        }
      
      }
    }
    
    void CDF_1988_S1865951::finalize() {
    ///@todo Total cross section hard-coded, needs a way to pass variable from pythia.
      scale(_hist_pt630, 32.6/sumOfWeights());
      scale(_hist_pt1800, 38.5/sumOfWeights());
    }
    //@}
    


  }
