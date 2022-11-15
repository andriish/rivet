// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief e+e-> Z+- pi-+
  class BESIII_2015_I1391798 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1391798);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9044213);
      // projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      // histograms
      _iEnergy=0;
      if      (isCompatibleWithSqrtS(4.23)) _iEnergy=1;
      else if (isCompatibleWithSqrtS(4.26)) _iEnergy=2;
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h[ix],2,ix+1,_iEnergy);
	book(_h[2+ix],3,1,1+ix);
      }
      book(_c,"TMP/counter");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles DD,other;
      for(const Particle & p : fs) {
	Particle parent=p;
	while(!parent.parents().empty()) {
	  parent=parent.parents()[0];
	  if(parent.abspid()==411 || parent.abspid()==413 ||
	     parent.abspid()==421 || parent.abspid()==423) break;
	}
	if((parent.abspid()==411 || parent.abspid()==421)
	   && !parent.parents().empty()) {
	  Particle Dstar = parent.parents()[0]; 
	  if(Dstar.abspid()==413 || Dstar.abspid()==423)
	    parent=Dstar;
	}
	if(parent.abspid()==411 || parent.abspid()==413 ||
	   parent.abspid()==421 || parent.abspid()==423) {
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
      // D Dbar + charged pion
      if(DD.size()!=2 || other.size()!=1) vetoEvent;
      if(DD[0].pid()*DD[1].pid()>0) vetoEvent;
      if(other[0].abspid()!=211) vetoEvent;
      if(DD[0].abspid()%10!=3) swap(DD[0],DD[1]);
      if(DD[0].abspid()%10!=3 || DD[1].abspid()%10!=1) vetoEvent;
      double mass = (DD[0].momentum()+DD[1].momentum()).mass();
      double cTheta = abs(axis.dot(other[0].momentum().p3().unit()));
      unsigned int iloc=0;
      // D0 D*- pi+ +cc
      if((DD[0].pid()==-413 && DD[1].pid()== 421 && other[0].pid()== 211) ||
      	 (DD[0].pid()== 413 && DD[1].pid()==-421 && other[0].pid()==-211)) {
       	iloc=0;
      }
      // D- D*0 pi+ +cc
      else  if((DD[0].pid()== 423 && DD[1].pid()==-411 && other[0].pid()== 211) ||
	       (DD[0].pid()==-423 && DD[1].pid()== 411 && other[0].pid()==-211)) {
      	iloc=1;
      }
      // otherwise veto event
      else
      	vetoEvent;
      _h[  iloc]->fill(mass);
      // parent Z+/-
      if(DD[0].parents()[0].abspid()==_pid && DD[1].parents()[0].abspid()==_pid &&
      	 fuzzyEquals(DD[0].parents()[0].momentum(),DD[1].parents()[0].momentum()) ) {
	_c->fill();
	_h[2+iloc]->fill(cTheta);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // distributions
      for(unsigned int ix=0;ix<4;++ix)
      	normalize(_h[ix],1.,false);
      // cross section
      double fact = crossSection()/ sumOfWeights() /picobarn;
      double sigma = fact*_c->val();
      double error = fact*_c->err();
      for(unsigned int iy=1;iy<_iEnergy+1;++iy) {
	Scatter2D temphisto(refData(1, 1, iy));
	Scatter2DPtr  mult;
	book(mult, 1, 1, iy);
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
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c;
    Histo1DPtr _h[4];
    int _pid;
    unsigned int _iEnergy;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1391798);

}
