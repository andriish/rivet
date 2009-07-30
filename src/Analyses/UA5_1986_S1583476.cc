#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Rivet.hh"
#include "Rivet/Event.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Analyses/UA5_1986_S1583476.hh"

namespace Rivet {
    
    /// Default constructor
    UA5_1986_S1583476::UA5_1986_S1583476()
    : Analysis("UA5_1986_S1583476") {
      const FinalState fs;
      addProjection(fs, "FS");
      addProjection(Beam(), "Beam");
    }

    /// @name Analysis methods
    //@{
    void UA5_1986_S1583476::init() { 
      _hist_eta200 = bookHistogram1D("d01-x01-y01", 20, 0., 5.);
      _hist_eta900 = bookHistogram1D("d01-x01-y03", 30, 0., 10.);
    }
    
    void UA5_1986_S1583476::analyze(const Event& event) {
      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
      const double weight = event.weight();
      foreach (const Particle& p, fs.particles())
      {
        if (fuzzyEquals(sqrtS, 200, 1E-4))
        {
        _hist_eta200->fill(p.momentum().pseudorapidity(), weight/27.9);
        }
        if (fuzzyEquals(sqrtS, 900))
        {
        _hist_eta900->fill(p.momentum().pseudorapidity(), weight/34.4);
        }
      }
    }
    
    void UA5_1986_S1583476::finalize() {
//Normalize using the range of pseudorapidity.
      normalize(_hist_eta200, 10);
      normalize(_hist_eta900, 18);
    }
    //@}

  }
