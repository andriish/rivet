// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 pi+ pi-
  class CLEO_2003_I633196 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2003_I633196);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");

      // Histograms
      book(_h_Kpim,1,1,1);
      book(_h_pipi,1,1,2);
      book(_h_Kpip,1,1,3);
      book(_dalitz, "dalitz",50,0.,3.,50,0.,2.);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode   = { { 310,1}, { 211,1},{-211,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	if (!D0.modeMatches(ix,3,mode) ) continue;
	int sign = D0.decaying()[ix].pid()/421;
	const Particles & pip= D0.decayProducts()[ix].at( sign*211);
	const Particles & pim= D0.decayProducts()[ix].at(-sign*211);
	const Particles & K0 = D0.decayProducts()[ix].at( 310);
	double mminus = (pim[0].momentum()+K0[0].momentum() ).mass2();
	double mplus  = (pip[0].momentum()+K0[0].momentum() ).mass2();
	double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	_h_Kpip->fill(mplus);
	_h_Kpim->fill(mminus);
	_h_pipi->fill(mpipi);
	_dalitz->fill(mminus,mpipi); 
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpip);
      normalize(_h_pipi);
      normalize(_h_Kpim);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpim,_h_pipi,_h_Kpip;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2003_I633196);

}
