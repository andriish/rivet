// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_b at 7 TeV
  class CMS_2012_I1113442 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2012_I1113442);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_lam   ,1,1,1);
      book(_h_bar   ,"TMP/bar_total",refData(1,1,2));
      book(_h_pT_lam,2,1,1);
      book(_h_pT_bar,"TMP/bar_pT"   ,refData(2,1,2));
      book(_h_y_lam ,3,1,1);
      book(_h_y_bar ,"TMP/bar_y"    ,refData(3,1,2));
      for(unsigned int ix=0;ix<2;++ix)
	book(_c_lam[ix],"TMP/c_Bc_"+toString(ix+1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==5122) ) {
	// stuff for the branching ratio
	double sign = p.pid()/p.abspid();
	if(p.children().size()==2 && 
	   ( (p.children()[0].pid()==sign*3122 && p.children()[1].pid()==443) ||
	     (p.children()[1].pid()==sign*3122 && p.children()[0].pid()==443) ))
	  _c_lam[0]->fill(1.);
	_c_lam[1]->fill(1.);
	double pT = p.perp();
	double y  = p.absrap();
	if(y>2. || pT<10.) continue;
	if(p.pid()>0) {
	  _h_lam->fill(sqrtS());
	  _h_pT_lam->fill(pT);
	  _h_y_lam ->fill(y);
	}
	else {
	  _h_bar->fill(sqrtS());
	  _h_pT_bar->fill(pT);
	  _h_y_bar ->fill(y);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      // total cross sectiohn and ratio
      scale(_h_lam,fact);
      scale(_h_bar,fact);
      Scatter2DPtr tmp;
      book(tmp,1,1,2);
      divide(_h_bar,_h_lam,tmp);
      // d sig / dpT
      fact *= 1000.;
      scale(_h_pT_lam,fact);
      scale(_h_pT_bar,fact);
      book(tmp,2,1,2);
      divide(_h_pT_bar,_h_pT_lam,tmp);
      // d sig / dy
      scale(_h_y_lam,fact);
      scale(_h_y_bar,fact);
      book(tmp,3,1,2);
      divide(_h_y_bar,_h_y_lam,tmp);
      // now branching ratio for cross sections
      if (_c_lam[1]->effNumEntries()>0.) {
	scale(_h_lam   , *_c_lam[0]/ *_c_lam[1]);
	scale(_h_pT_lam, *_c_lam[0]/ *_c_lam[1]);
	scale(_h_y_lam , *_c_lam[0]/ *_c_lam[1]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_lam,_h_bar;
    Histo1DPtr _h_pT_lam,_h_pT_bar;
    Histo1DPtr _h_y_lam,_h_y_bar;
    CounterPtr _c_lam[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2012_I1113442);

}
