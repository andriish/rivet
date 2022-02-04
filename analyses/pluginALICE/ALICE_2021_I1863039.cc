// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief Xi_c0 at 5 TeV
  class ALICE_2021_I1863039 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I1863039);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_Xi_pT[ix],1+ix,1,1);
	book(_h_D0_pT[ix],"TMP/pT_D0_"+toString(ix),refData(3+ix,1,1));
      }
      book(_h_Xi_pT[2],"TMP/pT_Xi",refData(4,1,1));
      book(_h_sig,5,1,1);
      book(_c_D0,"TMP/c_D0");
      book(_c_Xi,"TMP/c_Xi");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for( const Particle & p : ufs.particles(Cuts::pid==4132 || Cuts::pid==421)) {
	// no mixing
	if(p.children().size()==1 || p.children()[0].abspid()==p.abspid()) continue;
	// rapidity cut
	if(p.absrap()>0.5) continue;
	double pT = p.perp();
	if(p.pid()==421) {
	  // only prompt
	  if(!p.fromBottom()) {
	    _h_D0_pT[0]->fill(pT);
	    _h_D0_pT[1]->fill(pT);
	  }
	  _c_D0->fill();
	}
	else {
	  // prompt
	  if(!p.fromBottom()) {
	    _c_Xi->fill();
	    _h_sig->fill(1);
	    _h_Xi_pT[0]->fill(pT);
	  }
	  // inclusive
	  _h_Xi_pT[1]->fill(pT);
	  _h_Xi_pT[2]->fill(pT);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/microbarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_Xi_pT[ix],factor);
	scale(_h_D0_pT [ix] ,factor);
      }
      scale(_h_Xi_pT[2],factor);
      scale(_h_sig,factor);
      Scatter2DPtr tmp;
      // ratio prompt Xi0/D0
      book(tmp,3,1,1);
      divide(_h_Xi_pT[0],_h_D0_pT[0],tmp);
      // ratio  inclusive Xi0/D0
      book(tmp,4,1,1);
      divide(_h_Xi_pT[2],_h_D0_pT[1],tmp);
      // ratio Xi0_D0 integrated
      Scatter1D s1("/ALICE_2021_I1863039/TMP/r_1");
      Scatter1DPtr s1d = registerAO(s1);
      divide(_c_Xi,_c_D0,s1d);
      Scatter2D temphisto(refData(6, 1, 1));
      Scatter2DPtr ratio;
      book(ratio, 6, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(1, x-ex2.first, x+ex2.second)) {
	  ratio->addPoint(x,s1d->points()[0].x(),ex,s1d->points()[0].xErrs());
	}
	else {
	  ratio->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Xi_pT[3],_h_D0_pT[2];
    Histo1DPtr _h_sig;
    CounterPtr _c_D0,_c_Xi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I1863039);

}
