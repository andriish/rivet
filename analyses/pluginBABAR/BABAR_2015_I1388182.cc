// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief e+ e- > mu+ mu- gamma or pi+ pi- gamma
  class BABAR_2015_I1388182 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2015_I1388182);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(),"FS");
      // histograms
      for(unsigned int ix=0;ix<15;++ix) {
	Profile1DPtr tmp;
	book(tmp,3,1,1+ix);
	_h_pipi.push_back(tmp);
	if(ix>13) continue;
	book(tmp,1,1,1+ix);
	_h_mumu.push_back(tmp);
      }
      book(_A_mumu,2,1,1);
      book(_A_pipi,4,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming positron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis1 = beams.first .momentum().p3().unit();
      Vector3 axis2 = beams.second.momentum().p3().unit();
      if(beams.first.pid()<0) swap(axis1,axis2);
      // find the final state final state particles
      Particle mup,mum,pip,pim,gamma;
      Particles fs = apply<FinalState>(event,"FS").particles();
      // boost to CMF frame
      LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta((beams.first .momentum()+
									   beams.second.momentum()).betaVec());
      FourMomentum pGamma;
      for(const Particle & p : fs) {
	double theta = acos(p.p3().unit().dot(axis1));
	if(p.isCharged()) {
	  if(theta<.4 || theta>2.45) continue;
	  if(p.p3().mod()<1.) continue;
	  if(p.pid()==PID::MUON && mum.pid()!=PID::MUON)
	    mum = p;
	  else if(p.pid()==PID::ANTIMUON && mup.pid()!=PID::ANTIMUON)
	    mup = p;
	  else if(p.pid()==PID::PIPLUS  && pip.pid()!=PID::PIPLUS)
	    pip = p;
	  else if(p.pid()==PID::PIMINUS && pim.pid()!=PID::PIMINUS)
	    pim = p;
	}
	else if(p.pid()==PID::GAMMA) {
	  // angle cut on the photon
	  if(theta<.35 || theta>2.4) continue;
	  if(gamma.pid()!=PID::GAMMA) {
	    gamma = p;
	    pGamma = boost.transform(gamma.momentum());
	  }
	  else {
	    FourMomentum pNew = boost.transform(p.momentum());
	    if(pNew.E()>pGamma.E()) {
	      gamma = p;
	      pGamma = pNew;
	    }
	  }
	}
	else
	  vetoEvent;
      }
      if(gamma.pid()!=PID::GAMMA) vetoEvent;
      if(!( ((pip.pid()==PID::PIPLUS && pim.pid()==PID::PIMINUS ) ||
	     (mum.pid()==PID::MUON   && mup.pid()==PID::ANTIMUON))))
	vetoEvent;
      if (pip.pid()==PID::PIPLUS && mum.pid()==PID::MUON) vetoEvent;
      if(pGamma.E()<3.) vetoEvent;
      Vector3 axisZ = pGamma.p3().unit();
      Vector3 axisX = (axis2-axisZ.dot(axis2)*axisZ).unit();
      Vector3 axisY = axisZ.cross(axisX);
      FourMomentum pMinus,pPlus;
      if(mum.pid()==PID::MUON) {
	pMinus = boost.transform(mum.momentum());
	pPlus  = boost.transform(mup.momentum());
      }
      else {
	pMinus = boost.transform(pim.momentum());
	pPlus  = boost.transform(pip.momentum());
      }
      double phiM = atan2(pMinus.p3().dot(axisY),pMinus.p3().dot(axisX));
      if (phiM<0.) phiM+=2.*M_PI;
      double phiP = atan2(pPlus .p3().dot(axisY),pPlus .p3().dot(axisX));
      if (phiP<0.) phiP+=2.*M_PI;
      double mass = (pMinus+pPlus).mass();
      if(mum.pid()==PID::MUON) {
	if(mass>0.2 && mass<7.) {
	  unsigned int imass = int(mass/.5);
	  if(phiM<M_PI) {
	    _h_mumu[imass]->fill(cos(phiM), 1.);
	    _A_mumu->fill(mass, 1.5*cos(phiM));
	  }
	  else {
	    _h_mumu[imass]->fill(cos(phiP),-1.);
	    _A_mumu->fill(mass,-1.5*cos(phiP));
	  }
	}
      }
      else if(pip.pid()==PID::PIPLUS) {
	if(mass>0.3 && mass<1.8) {
	  unsigned int imass = int((mass-0.3)/.1);
	  if(phiM<M_PI) {
	    _h_pipi[imass]->fill(cos(phiM), 1.);
	    _A_pipi->fill(mass, 1.5*cos(phiM));
	  }
	  else {
	    _h_pipi[imass]->fill(cos(phiP),-1.);
	    _A_pipi->fill(mass,-1.5*cos(phiP));
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    }

    /// @}


    /// @name Histograms
    /// @{
    vector<Profile1DPtr> _h_mumu,_h_pipi;
    Profile1DPtr _A_mumu,_A_pipi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2015_I1388182);

}
