// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief J/psi and psi(2S) at 1.8 TeV
  class CDF_1997_I440101 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_1997_I440101);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_total[ix],1,1,1+ix);
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h_psi[ix][iy],2+2*ix+iy,1,1);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==443 or Cuts::pid==100443)) {
        double abseta = p.abseta();
        double xp = p.perp();
	if (xp<5. || abseta>0.6) continue;
	bool prompt = !p.fromBottom();
	unsigned int ipsi=p.pid()/100000;
	_h_total[ipsi]->fill(1.);
	_h_psi[prompt][ipsi]->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalisation factor
      double factor = crossSection()/nanobarn/sumOfWeights();
      // br to muons PDG 2021 (psi2s is e+e- due large errors on mu+mu-)
      vector<double> br = {0.05961,0.00793};
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_total[ix], factor*br[ix] );
	for(unsigned int iy=0;iy<2;++iy) {
	  scale(_h_psi[ix][iy],factor*br[iy]);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_psi[2][2],_h_total[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CDF_1997_I440101);

}
