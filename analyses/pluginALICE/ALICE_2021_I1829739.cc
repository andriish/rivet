// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_c+ at 5.02 TeV
  class ALICE_2021_I1829739 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2021_I1829739);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_Lambda,1,1,1);
      book(_h_D     ,"TMP/h_D",refData(4,1,1));
      book(_h_sig[0],7,1,1);
      book(_h_sig[1],7,1,2);
      book(_c_D,"TMP/c_D");
      book(_c_Lambda,"TMP/c_Lambda");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::abspid==4122 || Cuts::abspid==421)) {
	// prompt
	if(p.fromBottom()) continue;
	// skip copies due mixing
	if(p.children().size()==1 && p.children()[0].abspid()==p.abspid()) continue;
	if(p.absrap()>.5) continue;
	double pT = p.perp();
	if(p.abspid()==4122) {
	  _h_Lambda->fill(pT);
	  if(pT>1. && pT<12.) _h_sig[0]->fill(1);
	  _h_sig[1]->fill(1);
	  _c_Lambda->fill();
	}
	else {
	  _h_D->fill(pT);
	  _c_D->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/microbarn/sumOfWeights();
      scale(_h_Lambda,factor);
      scale(_h_D     ,factor);
      scale(_h_sig[0],factor);
      scale(_h_sig[1],factor);
      Scatter2DPtr tmp;
      // ratio prompt Lambda/D0
      book(tmp,4,1,1);
      divide(_h_Lambda,_h_D,tmp);
      // ratio lambda/D0 integrated
      Scatter1D s1("/ALICE_2021_I1829739/TMP/r_1");
      Scatter1DPtr s1d = registerAO(s1);
      divide(_c_Lambda,_c_D,s1d);
      Scatter2D temphisto(refData(8, 1, 1));
      Scatter2DPtr ratio;
      book(ratio, 8, 1, 1);
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
    Histo1DPtr _h_Lambda,_h_D;
    Histo1DPtr _h_sig[2];
    CounterPtr _c_D,_c_Lambda;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2021_I1829739);

}
