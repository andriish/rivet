// -*- C++ -*-
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet {

  /// Analysis
  class ALICE_2012_I1127497 : public Analysis {
  
  public:
    
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2012_I1127497);
    
    /// Book histograms and initialise projections before the run
    void init() {
      
      // Declare centrality projection
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");
      
      // Charged final states with |eta| < 0.5 and pT > 150 MeV
      const Cut& cut = Cuts::abseta < 0.5 && Cuts::pT > 150*MeV;
      const ChargedFinalState cfs(cut);
      addProjection(cfs,"CFS");
      
      // Loop over all histograms
      for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
	
	// Initialize PbPb objects
	_histNch[1][ihist] = bookHisto1D(ihist+1, 1, 1);
	_sumOfWeights[1][ihist] = 0;
	
	// Initialize pp objects. In principle, only one pp histogram would be needed since centrality does not 
	// make any difference here. However, in some cases in this analysis the binning differ from each other, 
	// so this is easy-to-implement way to account for that.
	std::string namePbPb = _histNch[1][ihist]->name();
	std::string namePP = namePbPb + "-pp";
	_histNch[0][ihist] = bookHisto1D( namePP, refData(ihist+1, 1, 1) ); // binning taken from ref data
	_sumOfWeights[0][ihist] = 0;
	
	// Initialize R_AA histograms
	_histRAA[ihist] = bookScatter2D(ihist+16, 1, 1);
	
      }
      
      // Centrality regions, 2 following numbers are boundaries for a certain region. Note, that some
      // regions overlap with other regions.
      _centrRegions.clear();
      _centrRegions += {{0., 5., 5., 10., 10., 20., 20., 30., 30., 40., 40., 50., 50., 60., 60., 70., 70., 80., 0., 10., 0., 20., 20., 40., 40., 60., 40., 80., 60., 80.}};
      

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const double weight = event.weight();
      
      // Final state particles with at least pT = 150 MeV in eta range of |eta| < 0.5
      const ChargedFinalState& charged = applyProjection<ChargedFinalState>(event, "CFS");
      Particles chargedParticles = charged.particlesByPt();
      
      // Flag to select which histograms to fill
      size_t pp_AA;
      
      // Vector of indices to fill the histograms with right centrality
      vector<size_t> indices;
      indices.clear();
      
      // Check type of event. This may be not the perfect way to check for the type of event as 
      // there might be some weird conditions hidden inside. For example some HepMC versions check 
      // if number of hard collisions is equal to 0 and assign 'false' in that case, which is usually wrong.
      // This might be changed in the future
      if (event.genEvent()->heavy_ion()) {
	
	// Prepare centrality projection and value
	const CentralityProjection& centrProj = apply<CentralityProjection>(event, "V0M");
	double centr = centrProj();
	// Veto event for too large centralities since those are not used in the analysis at all
	if ((centr < 0.) || (centr > 80.))
	  vetoEvent;
	
	// Flag as PbPb event
	pp_AA = 1;
	
	// Select indices
        for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
          // Check centrality bins and push corresponding indices of AA histograms.
          if (inRange(centr, _centrRegions[2*ihist], _centrRegions[2*ihist+1])) {
            indices.push_back(ihist);
          }
        }
	
      }
      else {
	
	// Flag as pp event
	pp_AA = 0;
	
	// Select indices (push all indices in case of pp histograms)
        for(size_t ihist = 0; ihist < NHISTOS; ++ihist)
          indices.push_back(ihist);

      }

      // Fill the right histograms and add weights based on the pp_AA flag and vector of indices
      for (size_t iindex = 0; iindex < indices.size(); ++iindex) {
        _sumOfWeights[pp_AA][indices.at(iindex)] += weight;
        foreach (const Particle& p, chargedParticles) {
          float pT = p.pT()/GeV;
          if (pT < 50.) {
	    _histNch[pp_AA][indices.at(iindex)]->fill(pT, weight*(1/pT));
          }
        }
      }
      
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      
      // Check for the reentrant finalize
      reentrant_flag = true;
      // For each event type
      for (int itype = 0; itype < 2; itype++) {
	// For each centrality range
	for (int ihist = 0; ihist < NHISTOS; ihist++) {
	  if (_histNch[itype][ihist]->numEntries() <= 0)
	    reentrant_flag = false;
	}
      }
      
      // Normal finalize
      if (reentrant_flag == false) {
	
	// Right scaling of the histograms with their individual weights.
	for (size_t pp_AA = 0; pp_AA < 2; ++pp_AA ) {
	  for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
	    if (_sumOfWeights[pp_AA][ihist] > 0.)
	      scale(_histNch[pp_AA][ihist], ( 1./_sumOfWeights[pp_AA][ihist] / 2. /  M_PI / 1.6 ) );
	  }
	  
	}
	
      }
      // Postprocessing of the histograms
      else if (reentrant_flag == true) {
	
	for( size_t i = 0; i < NHISTOS; ++i) {
	  divide( _histNch[1][i], _histNch[0][i], _histRAA[i] );
	}
	
      }
    }
    
  private:
    
    static const int NHISTOS = 15;
    
    Histo1DPtr _histNch[2][NHISTOS];
    double _sumOfWeights[2][NHISTOS];
    
    Scatter2DPtr _histRAA[NHISTOS];
    std::vector<float> _centrRegions;
    
    bool reentrant_flag;
    
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I1127497);


}
