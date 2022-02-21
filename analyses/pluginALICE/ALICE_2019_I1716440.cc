// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief D meson production at 5.02 and 7 TeV
  class ALICE_2019_I1716440 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2019_I1716440);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // histograms
      if (isCompatibleWithSqrtS(5020.0)) {
	for(unsigned int ix=0;ix<4;++ix) {
	  book(_h_prompt[ix],1+ix,1,1);
	  book(_h_ratio_num[ix],"TMP/h_ratio_num_"+toString(ix+1),refData(6+ix,1,1));
	  book(_h_ratio_den[ix],"TMP/h_ratio_den_"+toString(ix+1),refData(6+ix,1,1));
	}
	book(_h_incl,5,1,1);
	book(_h_total,18,1,1);
      }
      else if (isCompatibleWithSqrtS(7000.0)){
	for(unsigned int ix=0;ix<4;++ix) {
	  book(_h_ratio_num[ix],"TMP/h_ratio_num_"+toString(ix+1),refData(10+ix,1,1));
	  book(_h_ratio_den[ix],"TMP/h_ratio_den_"+toString(ix+1),refData(10+ix,1,1));
	}
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 5020 nor 7000 GeV.");
      for(unsigned int ix=0;ix<4;++ix) {
	book(_h_energy[ix],"h_energy_"+toString(ix+1),refData(14+ix,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==411 or Cuts::abspid==421 or
					     Cuts::abspid==431 or Cuts::abspid==413)) {
	// no mixing and |y|<0.5
	if(p.children().size()==1 || p.absrap()>0.5) continue;
	unsigned int imeson=0;
	if     (p.abspid()==411) imeson=1;
	else if(p.abspid()==413) imeson=2;
	else if(p.abspid()==431) imeson=3;
	double pT=p.perp();
	if(p.fromBottom()) {
	  if(imeson==0 && _h_incl) _h_incl->fill(pT);
	  continue;
	}
	// prompt at 5.02 TeV
	if (_h_prompt[imeson]) {
	  _h_prompt[imeson]->fill(pT);
	  _h_total->fill(imeson+1);
	}
	if(imeson==0) {
	  _h_ratio_den[0]->fill(pT);
	  _h_ratio_den[1]->fill(pT);
	  _h_ratio_den[2]->fill(pT);
	  if(_h_incl) _h_incl->fill(pT);
	}
	else if (imeson==1) {
	  _h_ratio_num[0]->fill(pT);
	  _h_ratio_den[3]->fill(pT);
	}
	else if(imeson==2) {
	  _h_ratio_num[1]->fill(pT);
	}
	else if(imeson==3) {
	  _h_ratio_num[2]->fill(pT);
	  _h_ratio_num[3]->fill(pT);
	}
	_h_energy[imeson]->fill(pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/microbarn/sumOfWeights();
      if(_h_prompt[0]) {
	for(unsigned int ix=0;ix<4;++ix)
	  scale(_h_prompt[ix],factor);
	scale(_h_total,factor);
	scale(_h_incl,factor);
      }
      int ioff(0);
      if (isCompatibleWithSqrtS(5020.0)) {
	ioff = 0;
      }
      else if (isCompatibleWithSqrtS(7000.0)) {
	ioff = 1;
      }
      for(unsigned int ix=0;ix<4;++ix) {
	Scatter2DPtr tmp;
	scale(_h_energy   [ix],factor);
	scale(_h_ratio_num[ix],factor);
	scale(_h_ratio_den[ix],factor);
	book(tmp,6+4*ioff+ix,1,1);
	divide(_h_ratio_num[ix],_h_ratio_den[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_prompt[4],_h_incl;
    Histo1DPtr _h_total;
    Histo1DPtr _h_ratio_num[4],_h_ratio_den[4];
    Histo1DPtr _h_energy[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2019_I1716440);

}
