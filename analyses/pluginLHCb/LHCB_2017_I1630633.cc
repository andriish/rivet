// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief B+ cross section at 7 and 13 TeV
  class LHCB_2017_I1630633 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2017_I1630633);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      int iloc=-1;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 0;
      }
      else if  (isCompatibleWithSqrtS(13000)) {
	iloc = 1;
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 7 or 13 TeV.");
      vector<double> ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      for(unsigned int iy=0;iy<5;++iy) {
	Histo1DPtr tmp;
	_h_B.add(ybins[iy],ybins[iy+1],book(tmp,1+iloc,1,1+iy));
      }
      book(_h_pT,3,1,1+iloc);
      book(_h_y ,4,1,1+iloc);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==521)) {
        double absrap = p.absrap();
	if(absrap<2. || absrap>4.5) continue;
        double pT = p.perp();
	if(pT>40.) continue;
	_h_B  .fill(absrap,pT);
	_h_pT->fill(pT);
	_h_y ->fill(absrap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1/2 due rapidity folding +/-
      double factor = 0.5*crossSection()/nanobarn/sumOfWeights();
      _h_B.scale(factor,this);
      scale(_h_pT,factor);
      scale(_h_y ,factor/1000.);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_B;
    Histo1DPtr _h_pT,_h_y;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2017_I1630633);

}
