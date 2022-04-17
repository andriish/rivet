// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> K0S K0S
  class BELLE_2013_I1245023 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2013_I1245023);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // histos
      if(inRange(sqrtS()/GeV,1.05,4.0)) {
	if(sqrtS()>1.1) book(_nKK[0],"TMP/nKK_1");
	book(_nKK[1],"TMP/nKK_2");
	double sMin=1.1, step=0.01;
	unsigned int ihist=2,iy=1;
	while(sMin<3.3) {
	  if(inRange(sqrtS()/GeV, sMin, sMin+step)) {
	    break;
	  }
	  sMin+=step;
	  iy+=1;
	  if(iy==4) {
	    ihist+=1;
	    iy=1;
	  }
	  if(fuzzyEquals(1.9, sMin)) step=0.02;
	  else if(fuzzyEquals(2.4, sMin)) step=0.04;
	  else if(fuzzyEquals(2.6, sMin)) step=0.1;
	}
	if(ihist<=40) book(_h_cTheta,ihist,1,iy);
      }
      else
       	throw Error("Invalid CMS energy for BELLE_2013_I1245023");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = applyProjection<FinalState>(event,"FS").particles();
      if(part.size()!=2) vetoEvent;
      for(const Particle & p : part) if (p.pid()!=PID::K0S) vetoEvent;
      double cTheta = abs(part[0].momentum().z()/part[0].momentum().p3().mod());
      if(cTheta<=0.6)          _nKK[1]->fill();
      if(cTheta<=0.8&&_nKK[0]) _nKK[0]->fill();
      if(_h_cTheta ) _h_cTheta ->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      if(_h_cTheta ) scale(_h_cTheta ,fact);
      for(unsigned int ix=0;ix<2;++ix) {
	if(!_nKK[ix]) continue;
	double sigma = _nKK[ix]->val()*fact;
	double error = _nKK[ix]->err()*fact;
	Scatter2D temphisto(refData(1, 1, 1+ix));
	Scatter2DPtr mult;
	book(mult, 1, 1, 1+ix);
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
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nKK[2];
    Histo1DPtr _h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2013_I1245023);

}
