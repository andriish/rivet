// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e-> Z(4025)+- pi-+
  class BESIII_2013_I1247675 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1247675);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 9044213);
      // projections
      declare(FinalState(), "FS");
      // histograms
      book(_h[0],1,1,1);
      book(_h[1],2,1,1);
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
      // D Dbar + charged pion
      if(DD.size()!=2 || other.size()!=1) vetoEvent;
      if(DD[0].pid()*DD[1].pid()>0) vetoEvent;
      if(other[0].abspid()!=211) vetoEvent;
      if((DD[0].pid()==-413 && DD[1].pid()== 423 && other[0].pid()== 211) ||
      	 (DD[0].pid()== 413 && DD[1].pid()==-423 && other[0].pid()==-211) ||
	 (DD[0].pid()== 423 && DD[1].pid()==-413 && other[0].pid()== 211) ||
	 (DD[0].pid()==-423 && DD[1].pid()== 413 && other[0].pid()==-211)) {
	double mass = (DD[0].momentum()+DD[1].momentum()).mass();
	_h[1]->fill(mass);
	// parent Z+/-
	if(DD[0].parents()[0].abspid()==_pid && DD[1].parents()[0].abspid()==_pid &&
	   fuzzyEquals(DD[0].parents()[0].momentum(),DD[1].parents()[0].momentum()) ) {
	  _h[0]->fill(4.26);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h[0],crossSection()/ sumOfWeights() /picobarn);
      normalize(_h[1],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    int _pid;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1247675);

}
