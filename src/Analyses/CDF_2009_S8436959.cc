// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  class CDF_2009_S8436959 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2009_S8436959()
      : Analysis("CDF_2009_S8436959") 
    {
      setBeams(PROTON, ANTIPROTON);
      setNeedsCrossSection(true);
    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;
      addProjection(fs, "FS");

      IdentifiedFinalState ifs(-1.0, 1.0, 30.0*GeV);
      ifs.acceptId(PHOTON);
      addProjection(ifs, "IFS");

      _h_Et_photon = bookHistogram1D(1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      ParticleVector photons;
      ParticleVector fs = applyProjection<IdentifiedFinalState>(event, "FS").particles();
      foreach (const Particle& photon, applyProjection<IdentifiedFinalState>(event, "IFS").particles()) {
        FourMomentum mom_in_cone;
        foreach (const Particle& p, fs) {
          if (deltaR(p.momentum(), photon.momentum()) < 0.4) {
            mom_in_cone += p.momentum();
          }
        }
        if (mom_in_cone.Et()-photon.momentum().Et() < 2.0*GeV) {
          photons.push_back(photon);
        }
      }
      if (photons.size() != 1) {
        vetoEvent;
      }
      
      _h_Et_photon->fill(photons[0].momentum().Et(), weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_Et_photon, crossSection()/sumOfWeights()/2.0);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_h_Et_photon;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_2009_S8436959> plugin_CDF_2009_S8436959;


}
