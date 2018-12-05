// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet {
  
  // Analysis
  class ALICE_2010_I880049 : public Analysis {
    
  public:
    
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2010_I880049);

    
    /// Book histograms and initialise projections before the run
    void init() {
      
      // Declare centrality projection
      declareCentrality(ALICE::V0MMultiplicity(),
        "ALICE_2015_PBPBCentrality", "V0M", "V0M");
      
      // Charged final states with |eta| < 0.5 and pT > 50 MeV
      const Cut& cut = Cuts::abseta < 0.5 && Cuts::pT > 50*MeV;
      const ChargedFinalState cfs(cut);
      addProjection(cfs,"CFS");
      
      // Histograms and variables initialization
      _histNchVsCentr = bookProfile1D(1, 1, 1);      
      _histNpartVsCentr = bookProfile1D(1, 1, 2);      
      
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");
      
      const CentralityProjection& centrProj = 
        apply<CentralityProjection>(event, "V0M");
      double centr = centrProj();
      double nch = charged.particles().size();
      if (centr > 80.)
        vetoEvent;
      _histNchVsCentr->fill(centr, nch, event.weight());
      // Attempt to extract Npart form GenEvent. TODO: Unclear how to handle this
      // in HepMC3
      const HepMC::HeavyIon* hi = event.genEvent()->heavy_ion();
      if (hi && hi->is_valid()) 
        _histNpartVsCentr->fill(centr, hi->Npart_proj() + hi->Npart_targ(),
	  event.weight());
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
      
    }

  private:
    
    Profile1DPtr _histNchVsCentr;
    Profile1DPtr _histNpartVsCentr;
  
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2010_I880049);
}
