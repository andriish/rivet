// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi0 pi0
  class BELLE_2009_I815978 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2009_I815978);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // histos
      if(inRange(sqrtS()/GeV,0.6,3.3) ||
	 inRange(sqrtS()/GeV,3.6,4.1)) {
	if(sqrtS()>0.72) book(_npipi[0],"TMP/npipi_1");
	book(_npipi[1],"TMP/npipi_2");
	double sMin=0.6, step=0.02;
	unsigned int ihist=1,iy=1;
	while(sMin<4.1) {
	  if(inRange(sqrtS()/GeV, sMin, sMin+step)) {
	    break;
	  }
	  sMin+=step;
	  iy+=1;
	  if(iy==4) {
	    ihist+=1;
	    iy=1;
	  }
	  if(fuzzyEquals(1.8, sMin)) step=0.04;
	  else if(fuzzyEquals(2.4, sMin)) step=0.1;
	  else if(fuzzyEquals(3.2, sMin)) sMin=3.6;
	}
	if(!inRange(sqrtS()/GeV,3.2,3.6)) book(_h_cTheta,ihist,1,iy);
      }
      else
       	throw Error("Invalid CMS energy for BELLE_2009_I815978");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = applyProjection<FinalState>(event,"FS").particles();
      if(part.size()!=2) vetoEvent;
      for(const Particle & p : part) if (p.pid()!=PID::PI0) vetoEvent;
      double cTheta = abs(part[0].momentum().z()/part[0].momentum().p3().mod());
      if(cTheta<=0.6&&_npipi[0]) _npipi[0]->fill();
      if(cTheta<=0.8&&_npipi[1]) _npipi[1]->fill();
      if(_h_cTheta ) _h_cTheta ->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      if(_h_cTheta ) scale(_h_cTheta ,fact);
      for(unsigned int ix=0;ix<2;++ix) {
	if(!_npipi[ix]) continue;
	double sigma = _npipi[ix]->val()*fact;
	double error = _npipi[ix]->err()*fact;
	Scatter2D temphisto(refData(31, 1, 1+ix));
	Scatter2DPtr mult;
	book(mult, 31, 1, 1+ix);
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
    CounterPtr _npipi[2];
    Histo1DPtr _h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2009_I815978);

}
