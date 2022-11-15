// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > B(*) B* pi
  class BELLE_2015_I1411223 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2015_I1411223);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles BB,other;
      for(const Particle & p : fs) {
      	Particle parent=p;
       	while(!parent.parents().empty()) {
      	  parent=parent.parents()[0];
      	  if(parent.abspid()==511 || parent.abspid()==521 ||
	     parent.abspid()==513 || parent.abspid()==523) break;
      	}
	if((parent.abspid()==511 || parent.abspid()==521)
	   && !parent.parents().empty()) {
	  Particle Bstar = parent.parents()[0]; 
	  if(Bstar.abspid()==513 || Bstar.abspid()==523)
	    parent=Bstar;
	}
	if(parent.abspid()==511 || parent.abspid()==521 ||
	   parent.abspid()==513 || parent.abspid()==523) {
	  bool found=false;
	  for (auto & B : BB) {
	    // B already in list 
	    if (fuzzyEquals(B.momentum(),parent.momentum())) {
	      found=true;
	      break;
	    }
	  }
	  if(!found) BB.push_back(parent);
	}
       	else {
       	  other.push_back(p);
       	}
      }
      // B Bbar + charged pion
      if(BB.size()!=2 || other.size()!=1) vetoEvent;
      if(BB[0].pid()*BB[1].pid()>0) vetoEvent;
      if(other[0].abspid()!=211) vetoEvent;
      if(BB[0].abspid()%10!=3) swap(BB[0],BB[1]);
      // B0 B*- pi+ +cc
      if((BB[0].pid()==-513 && BB[1].pid()== 521 && other[0].pid()==-211) ||
      	 (BB[0].pid()== 513 && BB[1].pid()==-521 && other[0].pid()== 211) ||
	 (BB[0].pid()== 523 && BB[1].pid()==-511 && other[0].pid()==-211) ||
	 (BB[0].pid()==-523 && BB[1].pid()== 511 && other[0].pid()== 211)) {
	_h[0]->fill(10.866);
      }
      else if((BB[0].pid()==-513 && BB[1].pid()== 523 && other[0].pid()==-211) ||
	      (BB[0].pid()== 513 && BB[1].pid()==-523 && other[0].pid()== 211) ||
	      (BB[0].pid()== 523 && BB[1].pid()==-513 && other[0].pid()==-211) ||
	      (BB[0].pid()==-523 && BB[1].pid()== 513 && other[0].pid()== 211)) {
      	_h[1]->fill(10.866);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	scale(_h[ix],crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2015_I1411223);

}
