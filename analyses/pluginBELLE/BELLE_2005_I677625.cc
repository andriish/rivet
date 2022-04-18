// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> p pbar
  class BELLE_2005_I677625 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2005_I677625);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // book histos
      if(inRange(sqrtS()/GeV,2.025,4.)) {
       	book(_cProton,"TMP/nProton");
	if(inRange(sqrtS()/GeV, 2.075, 2.9) ||
	   inRange(sqrtS()/GeV, 3.1, 4.)) {
	  double sMin=2.075, sMax=2.1, step=0.1;
	  unsigned int ihist=3,iy=1;
	  while(sMin<4.) {
	    if(inRange(sqrtS()/GeV, sMin, sMax)) break;
	    sMin=sMax;
	    sMax+=step;
	    iy+=1;
	    if(iy==4) {
	      ihist+=1;
	      iy=1;
	    }
	    if(fuzzyEquals(2.9, sMin)) {
	      sMin=3.1;
	      sMax=3.5;
	      step=0.5;
	    }
	  }
	  book(_h_cTheta[0],ihist,1,iy);
	}
	if(inRange(sqrtS()/GeV, 2.075, 2.5))
	  book(_h_cTheta[1],2,1,1);
	else if(inRange(sqrtS()/GeV, 2.5, 3.))
	  book(_h_cTheta[1],2,1,2);
	else if(inRange(sqrtS()/GeV, 3., 4.))
	  book(_h_cTheta[1],2,1,3);
      }
      else
	throw Error("Invalid CMS energy for BELLE_2005_I677625");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = applyProjection<FinalState>(event,"FS").particles();
      if(part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false),foundM(false);
      for(const Particle & p : part) {
	if(p.pid()==PID::PROTON) {
	  foundP=true;
	  cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
	}
	else if(p.pid()==PID::ANTIPROTON)
	  foundM=true;
      }
      if(!foundP || !foundM) vetoEvent;
      if(cTheta<=0.6)    _cProton->fill();
      if(_h_cTheta[0] ) _h_cTheta[0]->fill(cTheta);
      if(_h_cTheta[1] ) _h_cTheta[1]->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix)
	if(_h_cTheta[ix] ) scale(_h_cTheta[ix] ,fact);
      double sigma = _cProton->val()*fact;
      double error = _cProton->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr mult;
      book(mult, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(sqrtS(), x-ex2.first, x+ex2.second)) {
	  mult->addPoint(x, sigma, ex, make_pair(error,error));
	}
	else {
	  mult->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_cTheta[2];
    CounterPtr _cProton;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2005_I677625);

}
