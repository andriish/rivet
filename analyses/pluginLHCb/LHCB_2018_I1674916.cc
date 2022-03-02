// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief D_s asymmetry
  class LHCB_2018_I1674916 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2018_I1674916);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      vector<double> y = {2.,3.,3.5,4.5};
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  Histo1DPtr tmp;
	  _h[ix].add(y[iy],y[iy+1],book(tmp,"TMP/h_"+toString(ix+1)+"_"+toString(iy+1),
					refData(2,1,1+iy)));
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==431) ) {
	if(p.fromBottom()) continue;
	double pT = p.perp();
	double y  = p.absrap();
	bool anti = p.pid()<0;
	_h[anti].fill(y,pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      unsigned int is=0;
      if (isCompatibleWithSqrtS(7000)) {
      	is=1;
      }
      else if (isCompatibleWithSqrtS(8000)) {
      	is=2;
      }
      else  {
      	throw Error("Invalid CMS energy for LHCB_2018_I1674916");
      }
      // asymmetry
      for(unsigned int iy=0;iy<3;++iy) {
	Scatter2DPtr tmp;
	book(tmp,1+is,1,1+iy);
	asymm(_h[0].histos()[iy],_h[1].histos()[iy],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2018_I1674916);

}
