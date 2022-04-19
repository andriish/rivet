// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief 
  class CRYSTAL_BALL_1982_I168793 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1982_I168793);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // histos
      if(inRange(sqrtS()/GeV,0.72,1.8)) {
	book(_npipi,"TMP/npipi");
	if(inRange(sqrtS()/GeV,1.04,1.48)) book(_h_cTheta,2,1,1);
      }
      else
       	throw Error("Invalid CMS energy for CRYSTAL_BALL_1990_I294492");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = applyProjection<FinalState>(event,"FS").particles();
      if(part.size()!=2) vetoEvent;
      for(const Particle & p : part) if (p.pid()!=PID::PI0) vetoEvent;
      double cTheta = abs(part[0].momentum().z()/part[0].momentum().p3().mod());
      if(cTheta<=0.7) _npipi->fill();
      if(_h_cTheta ) _h_cTheta->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      if(_h_cTheta) scale(_h_cTheta,fact);
      double sigma = _npipi->val()*fact;
      double error = _npipi->err()*fact;
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
    CounterPtr _npipi;
    Histo1DPtr _h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1982_I168793);

}
