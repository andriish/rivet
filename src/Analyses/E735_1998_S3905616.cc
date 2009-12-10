// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class E735_1998_S3905616 : public Analysis {
  public:
 
    /// Constructor
    E735_1998_S3905616() : Analysis("E735_1998_S3905616") {
      setBeams(PROTON, ANTIPROTON);
    }
 

    /// @name Analysis methods
    //@{
 
    void init() {
      // Projection
      // const ChargedFinalState cfs();
      addProjection(ChargedFinalState(), "FS");

      // Histo
      _hist_multiplicity = bookHistogram1D(1, 1, 1);
    }


    void analyze(const Event& event) {
      const ChargedFinalState& fs = applyProjection<ChargedFinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      const double weight = event.weight();
      _hist_multiplicity->fill(numParticles, weight);
    }
 
 
    void finalize() {
      normalize(_hist_multiplicity);
    }
 
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_hist_multiplicity;
    //@}
 
  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<E735_1998_S3905616> plugin_E735_1998_S3905616;

}
