// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi and psi(2S) -> pi+pi-pi0
  class BESIII_2012_I1088606 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2012_I1088606);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid== 443 || Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::PI0);
      declare(PSI,"PSI");
      // histos
      book(_h_pipi[0],1,1,1);
      book(_h_pipi[1],1,1,2);
      book(_dalitz[0], "dalitz_Jpsi" ,50,0.,9.,50,0.0,9.);
      book(_dalitz[1], "dalitz_psi2S",50,0.,9.,50,0.0,9.);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1},{-211,1}, {111,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for(unsigned int ix=0;ix<PSI.decaying().size();++ix) {
	if (!PSI.modeMatches(ix,3,mode)) continue;
	const Particle & pip = PSI.decayProducts()[ix].at( 211)[0];
	const Particle & pim = PSI.decayProducts()[ix].at(-211)[0];
	const Particle & pi0 = PSI.decayProducts()[ix].at( 111)[0];
	double mminus = (pim.momentum()+pi0.momentum()).mass2();
	double mplus  = (pip.momentum()+pi0.momentum()).mass2();
	double mneut  = (pip.momentum()+pim.momentum()).mass2();
	unsigned int iloc = PSI.decaying()[ix].pid()/100000;
	_h_pipi[iloc]->fill(sqrt(mneut ));
	_h_pipi[iloc]->fill(sqrt(mplus ));
	_h_pipi[iloc]->fill(sqrt(mminus));
	_dalitz[iloc]->fill(mplus,mminus);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_pipi  [ix]);
	normalize(_dalitz[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pipi[2];
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2012_I1088606);

}
