// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/PrimaryParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/MixedFinalState.hh"
namespace Rivet {


  /// @brief Add a short analysis description here
  class ALICE_2016_I1507157 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2016_I1507157);


    /// @name Analysis methods
    //@{
        
    /// @brief Calculate angular distance between particles.
    double phaseDif(double a1, double a2){
      double dif = a1 - a2;
      while (dif < -M_PI/2)
        dif += 2*M_PI;
      while (dif > 3*M_PI/2)
        dif -= 2*M_PI;
      return dif;
    }

    /// Book histograms and initialise projections before the run
    void init() {

      double etamax = 0.8;
      double pTmin = 0.1; // GeV 

      // Charged tracks used to manage the mixing observable.
      ChargedFinalState cfsMult(Cuts::abseta < etamax);
      addProjection(cfsMult, "CFSMult");
      
      // Primary particles.

      PrimaryParticles pp({Rivet::PID::PIPLUS, Rivet::PID::KPLUS, 
	Rivet::PID::K0S, Rivet::PID::K0L, Rivet::PID::PROTON, 
	Rivet::PID::NEUTRON, Rivet::PID::LAMBDA, Rivet::PID::SIGMAMINUS,
       	Rivet::PID::SIGMAPLUS, Rivet::PID::XIMINUS, Rivet::PID::XI0, 
	Rivet::PID::OMEGAMINUS},Cuts::abseta < etamax && Cuts::pT > pTmin*GeV);
      addProjection(pp,"APRIM");

      // The event mixing projection
      declare(MixedFinalState(cfsMult, pp, 5, 0, 100, 10),"EVM");

      // Book histograms
      _h_protonRatio = bookScatter2D(2, 1, 2, true);
      _h_protonSignal = bookHisto1D("protonSignal",*_h_protonRatio,"protonSignal");
      _h_protonBackground = bookHisto1D("protonBackground",*_h_protonRatio,"protonBackground");
      nmp = 0;
      nsp = 0;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // The projections
      const PrimaryParticles& pp = applyProjection<PrimaryParticles>(event,"APRIM");
      const MixedFinalState& evm = applyProjection<MixedFinalState>(event, "EVM");

      // Get all mixing events
      vector<Particles> mixEvents = evm.getMixingEvents();
 
      // If there are not enough events to mix, don't fill any histograms
      if(mixEvents.size() == 0)
	      return;

      // Make a vector of mixed event particles
      vector<Particle> mixParticles;
      size_t pSize = 0;
      for(size_t i = 0; i < mixEvents.size(); ++i)
	      pSize+=mixEvents[i].size();
      mixParticles.reserve(pSize);
      for(size_t i = 0; i < mixEvents.size(); ++i)
      	      mixParticles.insert(mixParticles.end(), mixEvents[i].begin(), 
			      mixEvents[i].end());
	      
      // Shuffle the particles in the mixing events
      random_shuffle(mixParticles.begin(), mixParticles.end());
      
      for(const Particle& p1 : pp.particles()){
      	      // Start by doing the signal distribution
	      for(const Particle& p2 : pp.particles() ){
		      if(p1 == p2)
			   continue;
		      nsp+=1.0;
		      double dEta = abs(p1.eta() - p2.eta());
		      double dPhi = phaseDif(p1.phi(), p2.phi());
		      if(dEta < 1.3 && p1.pid() == 2212 && p2.pid() == -2212){
			      _h_protonSignal->fill(dPhi,weight);
		      }
	      }
	      // Then do the background distribution
	      for(const Particle& pMix : mixParticles){
	      	      nmp+=1.0;
		      double dEta = abs(p1.eta() - pMix.eta());
		      double dPhi = phaseDif(p1.phi(), pMix.phi());
		      if(dEta < 1.3 && p1.pid() == 2212 && pMix.pid() == -2212){
			      _h_protonBackground->fill(dPhi,weight);
	      		}
	      }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sc = nmp / nsp;
      scale(_h_protonSignal,sc);
      divide(_h_protonSignal,_h_protonBackground,_h_protonRatio);
      scale(_h_protonSignal,1.0/sc);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_protonSignal;
    Histo1DPtr _h_protonBackground;
    Scatter2DPtr _h_protonRatio;
    double nmp, nsp;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2016_I1507157);


}
