// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief B meson production at 7 TeV
  class LHCB_2013_I1238809 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2013_I1238809);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      vector<double> ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      for(unsigned int ib=0;ib<3;++ib) {
	book(_h_sigma[ib],1,1,1+ib);
	book(_h_pT[ib]   ,5+ib,1,1);
	book(_h_y [ib]   ,8+ib,1,1);
	for(unsigned int iy=0;iy<5;++iy) {
	  Histo1DPtr tmp;
	  _h_B[ib].add(ybins[iy],ybins[iy+1],book(tmp,2+ib,1,1+iy));
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==511 or Cuts::pid==521 or Cuts::pid==531 )) {
        double absrap = p.absrap();
	if(absrap<2. || absrap>4.5) continue;
	unsigned int ib=(p.pid()%100)/10-1;
        double pT = p.perp();
	if(pT>40.) continue;
	_h_sigma[ib]->fill(sqrtS());
	_h_B[ib].fill(absrap,pT);
	_h_pT[ib]->fill(pT);
	_h_y[ib]->fill(absrap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/microbarn/sumOfWeights();
      for(unsigned int ib=0;ib<3;++ib) {
	_h_B[ib].scale(factor,this);
	scale(_h_sigma[ib],factor);
	scale(_h_y    [ib],factor);
	scale(_h_pT   [ib],factor);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_sigma[3];
    BinnedHistogram _h_B[3];
    Histo1DPtr _h_pT[3],_h_y[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2013_I1238809);

}
