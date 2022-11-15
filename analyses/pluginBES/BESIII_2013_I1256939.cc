// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief e+e-> Z+- pi-+
  class BESIII_2013_I1256939 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1256939);


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
      book(_h[0],1,1,1);
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix+1],2,1,1+ix);
      book(_h[3],3,1,1);
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
      _h[1+iloc]->fill(mass);
      // parent Z+/-
      if(DD[0].parents()[0].abspid()==_pid && DD[1].parents()[0].abspid()==_pid &&
      	 fuzzyEquals(DD[0].parents()[0].momentum(),DD[1].parents()[0].momentum()) ) {
      	_h[0]->fill(4.26);
      	_h[3]->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h[0],crossSection()/ sumOfWeights() /picobarn);
      for(unsigned int ix=1;ix<4;++ix) {
      	normalize(_h[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    int _pid;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1256939);

}
