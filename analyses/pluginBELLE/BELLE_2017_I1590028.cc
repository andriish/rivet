// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+ e- > J/psi D Dbar
  class BELLE_2017_I1590028 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2017_I1590028);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      // histograms
      for(unsigned int ix=0;ix<6;++ix) {
	book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      FourMomentum pep;
      if(beams.first.pid()>0) {
	axis = beams.first .momentum().p3().unit();
	pep  = beams.second.momentum();
      }
      else {
	axis = beams.second.momentum().p3().unit();
	pep  = beams.first .momentum();
      }
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles DD,JPSI,other;
      for(const Particle & p : fs) {
	Particle parent=p;
	while(!parent.parents().empty()) {
	  parent=parent.parents()[0];
	  if(parent.abspid()==PID::DPLUS ||
	     parent.abspid()==PID::D0||
	     parent.pid()   ==PID::JPSI) break;
	}
	if(parent.abspid()!=PID::DPLUS &&
	   parent.abspid()!=PID::D0    &&
	   parent.pid()   !=PID::JPSI) {
	  other.push_back(p);
	  continue;
	}
	bool found=false;
	for (auto & D : parent.pid()==PID::JPSI ? JPSI : DD) {
	  // D or J/psi already in list 
	  if (fuzzyEquals(D.momentum(),parent.momentum())) {
	    found=true;
	    break;
	  }
	}
	if(!found) {
	  (parent.pid()==PID::JPSI ? JPSI : DD).push_back(parent);
	}
      }
      // require J/psi D Dbar
      if(JPSI.size()!=1 || DD.size()!=2 || other.size()!=0) vetoEvent;
      if(DD[0].pid()!=-DD[1].pid()) vetoEvent;
      FourMomentum pDD = DD[0].momentum()+DD[1].momentum(); 
      _h[0]->fill(pDD.mass());
      double cProd = axis.dot(JPSI[0].momentum().p3().unit());
      _h[1]->fill(cProd);
      // finally the leptons from J/psi decay
      if(JPSI[0].children().size()!=2) vetoEvent;
      if(JPSI[0].children()[0].pid()!=-JPSI[0].children()[1].pid()) vetoEvent;
      if(JPSI[0].children()[0].abspid()!=PID::EMINUS &&
	 JPSI[0].children()[0].abspid()!=PID::MUON) vetoEvent;
      Particle lm = JPSI[0].children()[0];
      Particle lp = JPSI[0].children()[1];
      // variables in J/psi frame
      Vector3 axis2 = -pDD.p3().unit();
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(JPSI[0].momentum().betaVec());
      FourMomentum plm = boost1.transform(lm.momentum());
      Vector3 axis4 = plm.p3().unit();
      double cPsi = axis2.dot(axis4);
      _h[2]->fill(cPsi);
      Vector3 axis3 = boost1.transform(pep).p3().unit();
      Vector3 aep = axis3-axis3.dot(axis2)*axis2;
      Vector3 alm = axis4-axis4.dot(axis2)*axis2;
      double phil = atan2(-axis2.cross(aep).dot(alm),aep.dot(alm));
      _h[4]->fill(phil);
      // variables in DD frame
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDD.betaVec());
      FourMomentum pD = boost2.transform((DD[0].pid()>0 ? DD[0] : DD[1]).momentum());
      Vector3 axisD = pD.p3().unit();
      axis2 *= -1;
      double cX = axis2.dot(axisD);
      _h[3]->fill(cX);
      axis3 = boost2.transform(pep).p3().unit();
      aep = axis3-axis3.dot(axis2)*axis2;
      alm = axisD-axisD.dot(axis2)*axis2;
      phil = atan2(-axis2.cross(aep).dot(alm),aep.dot(alm));
      _h[5]->fill(phil);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<6;++ix) {
	normalize(_h[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[6];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BELLE_2017_I1590028);

}
