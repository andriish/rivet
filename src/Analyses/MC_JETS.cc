// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class MC_JETS : public MC_JetAnalysis {
  public:

    MC_JETS() : MC_JetAnalysis("MC_JETS", 4, "Jets") 
    {
      setNeedsCrossSection(true);
    }


  public:

    void init() {
      FinalState fs;
      FastJets jetpro(fs, FastJets::KT, 0.7, 20.0*GeV);
      addProjection(jetpro, "Jets");
      
      MC_JetAnalysis::init();
    }


    void analyze(const Event& event) {
      MC_JetAnalysis::analyze(event);
    }


    void finalize() {
      MC_JetAnalysis::finalize();
    }

  };

  AnalysisBuilder<MC_JETS> plugin_MC_JETS;

}
