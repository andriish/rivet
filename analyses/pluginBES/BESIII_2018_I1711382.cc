// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+ e- -> Ds+ Dbar0(*) K-
  class BESIII_2018_I1711382 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1711382);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      // histograms
      for(unsigned int ix=0;ix<4;++ix)
	book(_h_sigma[ix],1,1,1+ix);
      book(_h_cTheta,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles DS,D0,other;
      // first find the kaon and Ds
      for(const Particle & p : fs) {
	Particle parent=p;
	while(!parent.parents().empty()) {
	  parent=parent.parents()[0];
	  if(parent.abspid()==PID::DSPLUS ||
	     parent.abspid()==PID::D0) break;
	}
	if(parent.abspid()!=PID::DSPLUS &&
	   parent.abspid()!=PID::D0 ) {
	  other.push_back(p);
	  continue;
	}
	bool found=false;
	for (auto & D : parent.abspid()==PID::DSPLUS ? DS : D0) {
	  // D already in list 
	  if (fuzzyEquals(D.momentum(),parent.momentum())) {
	    found=true;
	    break;
	  }
	}
	if(!found) {
	  (parent.abspid()==PID::DSPLUS ? DS : D0).push_back(parent);
	}
      }
      // Ds and D0, one particle and one anti particle
      if(DS.size()!=1 || D0.size()!=1 || DS[0].pid()*D0[0].pid()>0) vetoEvent;
      int iK = DS[0].pid()>0 ? -321 : 321;
      if(other.size()==1 && other[0].pid()==iK) {
	_h_sigma[1]->fill(4.6);
	if(other[0].parents()[0].abspid()==435 &&
	   D0   [0].parents()[0].abspid()==435) {
	  _h_sigma[3]->fill(4.6);
	  FourMomentum pDs2 = D0[0].momentum()+other[0].momentum();
	  LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(pDs2.betaVec());
	  Vector3 axis1 = pDs2.p3().unit();
	  Vector3 axis2 = boost.transform(other[0].momentum()).p3().unit();
	  _h_cTheta->fill(abs(axis1.dot(axis2)));
	}
      }
      else if(D0[0].parents()[0].abspid()==423) {
	bool Dstar = true;
	for( const Particle & p : other) {
	  if(p.parents()[0].abspid()!=423 && p.pid()!=iK) {
	    Dstar = false;
	    break;
	  }
	}
	if(Dstar) {
	  _h_sigma[0]->fill(4.6);
	  if(D0[0].parents()[0].parents()[0].abspid()==10433)
	    _h_sigma[2]->fill(4.2);
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<4;++ix)
	scale(_h_sigma[ix],fact);
      normalize(_h_cTheta);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_sigma[4],_h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1711382);

}
