// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+ e- -> pi0 pi0 J/psi
  class BESIII_2020_I1793431 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1793431);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9030443);
      // projections
      declare(FinalState(), "FS");
      // histos
      unsigned int iloc=0;
      if      (isCompatibleWithSqrtS(4.226)) iloc=1;
      else if (isCompatibleWithSqrtS(4.236)) iloc=2;
      else if (isCompatibleWithSqrtS(4.244)) iloc=3;
      else if (isCompatibleWithSqrtS(4.258)) iloc=4;
      if(iloc>0) {
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[iy],2,iloc,1+iy);
      }
      for(unsigned int iy=0;iy<3;++iy)
	book(_c[iy],"TMP/c_"+toString(iy+1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles Jpsi,other;
      for(const Particle & p : fs) {
	Particle parent=p;
	while(!parent.parents().empty()) {
	  parent=parent.parents()[0];
	  if(parent.abspid()==PID::JPSI) break;
	}
	if(parent.abspid()!=PID::JPSI) {
	  other.push_back(p);
	  continue;
	}
	bool found=false;
	for (auto & psi : Jpsi) {
	  // J/psi already in list 
	  if (fuzzyEquals(psi.momentum(),parent.momentum())) {
	    found=true;
	    break;
	  }
	}
	if(!found) Jpsi.push_back(parent);
      }
      if(Jpsi.size()!=1 || other.size()!=2) vetoEvent;
      if(other[0].pid()==PID::PI0 && other[1].pid()==PID::PI0) {
	_c[0]->fill();
	if(Jpsi[0].parents()[0].pid()==_pid) _c[2]->fill();
	if(_h[0]) {
	  for(unsigned int iy=0;iy<2;++iy)
	    _h[0]->fill((Jpsi[0].momentum()+other[iy].momentum()).mass());
	  _h[1]->fill((other[0].momentum()+other[1].momentum()).mass());
	}
      }
      else if(other[0].pid()==-other[1].pid() && other[0].abspid()==PID::PIPLUS) {
	_c[1]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=0;iy<2;++iy)
	if(_h[iy]) normalize(_h[iy],1.,false);
      double fact = crossSection()/ sumOfWeights() /picobarn;;
      for(unsigned int ii=0;ii<3;++ii) {
	double sigma;
	pair<double,double> error;
	unsigned int ix,iy;
	if(ii!=1) {
	  sigma = fact*_c[ii]->val();
	  error = make_pair(fact*_c[ii]->err(),
			    fact*_c[ii]->err());
	  iy=1;
	  ix=ii+1;
	}
	else {
	  ix=1;
	  iy=2;
	  Scatter1D R = *_c[0]/ *_c[1];
	  sigma = R.point(0).x();
	  error = R.point(0).xErrs();
	}
	Scatter2D temphisto(refData(ix, 1, iy));
	Scatter2DPtr  mult;
	book(mult, ix, 1, iy);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	    mult->addPoint(x, sigma, ex, error);
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
    int _pid;
    Histo1DPtr _h[2];
    CounterPtr _c[3];
    
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1793431);

}
