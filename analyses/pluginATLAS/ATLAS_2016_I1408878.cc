// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief 
  class ATLAS_2016_I1408878 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1408878);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<3;++ix)
	book(_h_tot[ix],1,1,1+ix);
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_pT[ix],2,1,1+ix);
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h_y[iy][ix],3+ix,1,1+iy);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==411 or Cuts::abspid==413 or Cuts::abspid==431)) {
	if(p.children().size()==1) continue;
	double eta = p.abseta();
	if (eta>2.1) continue;
	double pT = p.perp();
	if(pT<3.5 || pT>100.) continue;
	unsigned int itype=0;
	if (p.abspid()==411)
	  itype=1;
	else if(p.abspid()==431)
	  itype=2;
	_h_tot[itype]->fill(pT);
	if(itype<2) {
	  _h_pT[itype]->fill(pT);
	  if(pT<20.)
	    _h_y[itype][0]->fill(eta);
	  else
	    _h_y[itype][1]->fill(eta);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/microbarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_tot[ix],factor);
	scale(_h_pT[ix],factor);
	scale(_h_y[ix][0],factor);
	scale(_h_y[ix][1],1000.*factor);
      }
      scale(_h_tot[2],factor);
      // really total cross section in bin not differential so scale by binwidth
      for(unsigned int ix=0;ix<3;++ix) {
	for(size_t b = 0; b < _h_tot[ix]->numBins(); ++b) {
	  _h_tot[ix]->bins()[b].scaleW(_h_tot[ix]->bins()[b].xWidth());
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[2],_h_y[2][2],_h_tot[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1408878);

}
