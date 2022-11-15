// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> → D* D* pi0
  class BESIII_2015_I1382146 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1382146);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9030443);
      // projections
      declare(FinalState(), "FS");
      // histograms
      // histograms
      unsigned int iloc=0;
      if      (isCompatibleWithSqrtS(4.23)) iloc=0;
      else if (isCompatibleWithSqrtS(4.26)) iloc=1;
      book(_h[0],2,1,iloc+2);
      book(_h[1],2,1,1);
      book(_c,"TMP/counter");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles DD,other;
      for(const Particle & p : fs) {
	Particle parent=p;
	while(!parent.parents().empty()) {
	  parent=parent.parents()[0];
	  if(parent.abspid()==413 || parent.abspid()==423) break;
	}
	if(parent.abspid()==413 || parent.abspid()==423) {
	  bool found=false;
	  for (auto & D : DD) {
	    // D already in list 
	    if (fuzzyEquals(D.momentum(),parent.momentum())) {
	      found=true;
	      break;
	    }
	  }
	  if(!found) DD.push_back(parent);
	}
	else {
	  other.push_back(p);
	}
      }
      // D Dbar + neutral pion
      if(DD.size()!=2 || other.size()!=1) vetoEvent;
      if(DD[0].pid()!=-DD[1].pid()) vetoEvent;
      if(other[0].pid()!=111) vetoEvent;
      double mass = (DD[0].momentum()+DD[1].momentum()).mass();
      _h[0]->fill(mass);
      _h[1]->fill(mass);
      // parent Z0
      if(DD[0].parents()[0].pid()==_pid && DD[1].parents()[0].pid()==_pid &&
	 fuzzyEquals(DD[0].parents()[0].momentum(),DD[1].parents()[0].momentum()) ) _c->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // distributions
      for(unsigned int ix=0;ix<2;++ix)
      	normalize(_h[ix],1.,false);
      // cross section
      double fact = crossSection()/ sumOfWeights() /picobarn;
      double sigma = fact*_c->val();
      double error = fact*_c->err();
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr  mult;
      book(mult, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
      	const double x  = temphisto.point(b).x();
      	pair<double,double> ex = temphisto.point(b).xErrs();
      	pair<double,double> ex2 = ex;
      	if(ex2.first ==0.) ex2. first=0.0001;
      	if(ex2.second==0.) ex2.second=0.0001;
      	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
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
    CounterPtr _c;
    Histo1DPtr _h[2];
    int _pid;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1382146);

}
