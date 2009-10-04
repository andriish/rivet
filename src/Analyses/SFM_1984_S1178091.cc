// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  class SFM_1984_S1178091 : public Analysis {
  public:

    /// Constructor
    SFM_1984_S1178091() : Analysis("SFM_1984_S1178091") {
      setBeams(PROTON, PROTON);
    }


    /// @name Analysis methods
    //@{

    void init() {
      // Projections
      addProjection(Beam(), "Beam");
      addProjection(ChargedFinalState(), "FS");

      // Histograms
      _hist_multiplicity_inel_30 = bookHistogram1D(1, 1, 1); 
      _hist_multiplicity_inel_45 = bookHistogram1D(1, 1, 2);
      _hist_multiplicity_inel_53 = bookHistogram1D(1, 1, 3);
      _hist_multiplicity_inel_63 = bookHistogram1D(1, 1, 4);
      _hist_multiplicity_nsd_30 = bookHistogram1D(2, 1, 1);
      _hist_multiplicity_nsd_45 = bookHistogram1D(2, 1, 2);
      _hist_multiplicity_nsd_53 = bookHistogram1D(2, 1, 3);
      _hist_multiplicity_nsd_63 = bookHistogram1D(2, 1, 4);  
    }
    
    
    void analyze(const Event& event) {
      const double weight = event.weight();
      const double sqrtS = applyProjection<Beam>(event, "Beam").sqrtS();
      const ChargedFinalState& fs = applyProjection<ChargedFinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      
      // Decide whether event is of diffractive type or not 
      // FIXME: it is not so clear in the paper how this distinction is made.
      // They seem to require either exactly one particle with Feynman x larger
      // than 0.8 to call an event diffractive or that there are no tracks
      // reconstructed in either of the two hemispheres. For the latter
      // they require in addition also the number of cahrged particles
      // to be smaller than 8.
      
      int n_left(0), n_right(0), n_large_x(0);
      foreach (const Particle& p, fs.particles()) {
        // Calculate the particles' Feynman x  
        const double x_feyn = 2.0 * (p.momentum().pz()/GeV) / sqrtS;
        if (fabs(x_feyn) > 0.8 ) n_large_x += 1;
        
        // Pseudorapidity
        const double eta = p.momentum().pseudorapidity();
        if (eta > 0.0) n_right += 1;
        else if (eta < 0.0) n_left += 1;
      }
      
      // Not sure about the "=="
      /// @todo Numerical precision problem!
      bool isDiffractive = false;
      if (n_large_x == 1) isDiffractive = true;
      
      // FIXME: Not sure about the "== 1", the paper says no charged particle
      // that was reconstructed so the incoming protons must run down the beam
      // pipe. Since we look a the complete final state here no particle being
      // reconstructed should be equal to one particle (proton) in each
      // hemisphere.  The "< 8" is also not certain.
      if ((n_left == 1 || n_right == 1) && numParticles < 8 ) {
        isDiffractive = true;
      }
      
      getLog() << Log::DEBUG << "N_left: " << n_left << ", N_right: " 
               << n_right << ", N_large_x: " << n_large_x << endl;
      
      
      // Fill histos of charged multiplicity distributions
      // The inelastic samples are said to contain also diffractive events.
      //
      if (fuzzyEquals(sqrtS, 30.4/GeV, 1E-1)) {
        if (isDiffractive) {
          _hist_multiplicity_nsd_30 ->fill(numParticles, weight);
          _hist_multiplicity_inel_30->fill(numParticles, weight);
        } else {
          _hist_multiplicity_inel_30->fill(numParticles, weight);
        }  
      } 
      else if (fuzzyEquals(sqrtS, 44/GeV, 1E-1)) {
        if (isDiffractive) {
          _hist_multiplicity_nsd_45 ->fill(numParticles, weight);
          _hist_multiplicity_inel_45->fill(numParticles, weight);
        } else {
          _hist_multiplicity_inel_45->fill(numParticles, weight);
        }  
      }
      else if (fuzzyEquals(sqrtS, 53/GeV, 1E-1)) {
        if (isDiffractive) {
          _hist_multiplicity_nsd_53 ->fill(numParticles, weight);
          _hist_multiplicity_inel_53->fill(numParticles, weight);
        } else {
          _hist_multiplicity_inel_53->fill(numParticles, weight);
        }  
      }
      else if (fuzzyEquals(sqrtS, 63/GeV, 1E-1)) {
        if (isDiffractive) {
          _hist_multiplicity_nsd_63 ->fill(numParticles, weight);
          _hist_multiplicity_inel_63->fill(numParticles, weight);
        }
        else {
          _hist_multiplicity_inel_63->fill(numParticles, weight);
        }  
      }
      
    }
    
    
    void finalize() {
      normalize(_hist_multiplicity_inel_30);
      normalize(_hist_multiplicity_inel_45);
      normalize(_hist_multiplicity_inel_53);
      normalize(_hist_multiplicity_inel_63);
      normalize(_hist_multiplicity_nsd_30 );
      normalize(_hist_multiplicity_nsd_45 );
      normalize(_hist_multiplicity_nsd_53 );
      normalize(_hist_multiplicity_nsd_63 );
    }
    //@}
    

  private:
    
    /// @name Histograms
    //@{

    AIDA::IHistogram1D *_hist_multiplicity_inel_30;
    AIDA::IHistogram1D *_hist_multiplicity_inel_45;
    AIDA::IHistogram1D *_hist_multiplicity_inel_53;
    AIDA::IHistogram1D *_hist_multiplicity_inel_63;
    AIDA::IHistogram1D *_hist_multiplicity_nsd_30;
    AIDA::IHistogram1D *_hist_multiplicity_nsd_45;
    AIDA::IHistogram1D *_hist_multiplicity_nsd_53;
    AIDA::IHistogram1D *_hist_multiplicity_nsd_63;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<SFM_1984_S1178091> plugin_SFM_1984_S1178091;

}
