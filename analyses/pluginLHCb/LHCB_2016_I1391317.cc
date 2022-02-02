// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Lambda_b and B0 production at 7 and 8 TeV
  class LHCB_2016_I1391317 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2016_I1391317);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // CMS energy
      int iloc=-1;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 0;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 1;
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 7 or 8 TeV.");
      // histograms
      vector<double> ybins={2.0,2.5,3.0,3.5,4.0,4.5};
      for(unsigned int iy=0;iy<5;++iy) {
	Histo1DPtr tmp;
	_h_Lambda.add(ybins[iy],ybins[iy+1],book(tmp,1+iloc,1,1+iy));
	_h_B     .add(ybins[iy],ybins[iy+1],book(tmp,3+iloc,1,1+iy));
      }
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_pT[ix],"TMP/h_pT_"+toString(ix),refData(5,1,1+iloc));
	book(_h_y [ix],"TMP/h_y_ "+toString(ix),refData(6,1,1+iloc));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::abspid==5122 || Cuts::abspid==511)) {
	// skip copies due mixing
	if(p.children().size()==1 && p.children()[0].abspid()==p.abspid()) continue;
	double rap=p.rapidity();
	if(rap<2. || rap>4.5) continue;
	double pT = p.perp();
	if(p.abspid()==5122) {
	  _h_Lambda.fill(rap,pT);
	  if(p.pid()>0) {
	    _h_pT[0]->fill(pT);
	    if(pT<20.) _h_y [0]->fill(rap);
	  }
	  else {
	    _h_pT[1]->fill(pT);
	    if(pT<20.) _h_y [1]->fill(rap);
	  }
	}
	else {
	  _h_B.fill(rap,pT);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // CMS energy
      int iloc=-1;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 0;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 1;
      }
      // branching ratios
      vector<double> br={3.2e-4,1.27e-3};
      // 0.5 particle/antiparticle
      double factor = 0.5*crossSection()/picobarn/sumOfWeights();
      _h_Lambda.scale(br[0]*factor,this);
      _h_B     .scale(br[1]*factor,this);
      Scatter2DPtr tmp;
      book(tmp,5,1,1+iloc);
      asymm(_h_pT[0],_h_pT[1],tmp);
      book(tmp,6,1,1+iloc);
      asymm(_h_y [0],_h_y [1],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_Lambda,_h_B;
    Histo1DPtr _h_pT[2],_h_y[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2016_I1391317);

}
