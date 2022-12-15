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
      Profile1DPtr tmp;
      book(tmp,"TMP/pi0",refData(3,1,1));
      _h_pipi.push_back(tmp);
      for(unsigned int ix=0;ix<15;++ix) {
	book(tmp,3,1,1+ix);
	_h_pipi.push_back(tmp);
	if(ix>13) continue;
	book(tmp,1,1,1+ix);
	_h_mumu.push_back(tmp);
      }
      book(tmp,"TMP/pi16",refData(3,1,15));
      _h_pipi.push_back(tmp);
      book(tmp,"TMP/pi17",refData(3,1,15));
      _h_pipi.push_back(tmp);
      book(_h_nopsi,"TMP/nopsi",refData(1,1,7));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming positron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      bool CMF = fuzzyEquals(-beams.first .momentum().z()/beams.second.momentum().z(),1);
      Vector3 axis1 = beams.first .momentum().p3().unit();
      Vector3 axis2 = beams.second.momentum().p3().unit();
      if(beams.first.pid()<0) swap(axis1,axis2);
      // find the final state final state particles
      Particle mup,mum,pip,pim,gamma;
      Particles fs = apply<FinalState>(event,"FS").particles();
      // boost to CMF frame
      LorentzTransform boost;
      if(!CMF) {
	boost = LorentzTransform::mkFrameTransformFromBeta((beams.first .momentum()+
							    beams.second.momentum()).betaVec());
      }
      else {
	double E1=3.5,E2 = 0.25*sqr(sqrtS())/E1;
	FourMomentum pnew(E1+E2,0,0,E1-E2);
	boost = LorentzTransform::mkFrameTransformFromBeta(-pnew.betaVec());
      }
      FourMomentum pGamma;
      for(const Particle & p : fs) {
	FourMomentum pLab,pCMF;
	if(CMF) {
	  pCMF = p.momentum();
	  pLab = boost.transform(p.momentum());
	}
	else {
	  pCMF = boost.transform(p.momentum());
	  pLab = p.momentum();
	}
	double theta  = acos(pLab.p3().unit().dot(axis1));
	if(p.isCharged()) {
	  if(theta<.4 || theta>2.45) continue;
	  if(pLab.p3().mod()<1.) continue;
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
	    pGamma = pCMF;
	  }
	  else {
	    if(pCMF.E()>pGamma.E()) {
	      gamma = p;
	      pGamma = pCMF;
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
      if(CMF) {
	if(mum.pid()==PID::MUON) {
	  pMinus = mum.momentum();
	  pPlus  = mup.momentum();
	}
	else {
	  pMinus = pim.momentum();
	  pPlus  = pip.momentum();
	}
      }
      else {
	if(mum.pid()==PID::MUON) {
	  pMinus = boost.transform(mum.momentum());
	  pPlus  = boost.transform(mup.momentum());
	}
	else {
	  pMinus = boost.transform(pim.momentum());
	  pPlus  = boost.transform(pip.momentum());
	}
      }
      double phiM = atan2(pMinus.p3().dot(axisY),pMinus.p3().dot(axisX));
      if (phiM<0.) phiM+=2.*M_PI;
      double phiP = atan2(pPlus .p3().dot(axisY),pPlus .p3().dot(axisX));
      if (phiP<0.) phiP+=2.*M_PI;
      double mass = (pMinus+pPlus).mass();
      if(mum.pid()==PID::MUON) {
	if(mass>0.2 && mass<7.) {
	  unsigned int imass = int(mass/.5);
	  if(phiM<M_PI) _h_mumu[imass]->fill(cos(phiM),-1.);
	  else          _h_mumu[imass]->fill(cos(phiP), 1.);
	  if(imass==6 && mass>3.2) {
	    if(phiM<M_PI) _h_nopsi->fill(cos(phiM),-1.);
	    else          _h_nopsi->fill(cos(phiP), 1.);
	  }
	}
      }
      else if(pip.pid()==PID::PIPLUS) {
	if(mass>0.2 && mass<2.0) {
	  unsigned int imass = int((mass-0.2)/.1);
	  if(phiM<M_PI) _h_pipi[imass]->fill(cos(phiM),-1.);
	  else          _h_pipi[imass]->fill(cos(phiP), 1.);
	}
      }
    }

    pair<double,double> calcA0(Profile1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (auto bin : hist->bins() ) {
	if(bin.numEntries()<2) continue;
	double Oi = bin.mean();
	double Ei = bin.stdErr();
	if(Ei==0.) continue;
	double xi = 0.5*(bin.xMin()+bin.xMax());
	sum1 += sqr(xi/Ei);
	sum2 += Oi*xi/sqr(Ei);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
      // muon assymetries
      Scatter2DPtr _A_mumu;
      book(_A_mumu,2,1,1);
      for(unsigned int ix=0;ix<14;++ix) {
	double xmin = 0.5*double(ix);
	double xmax = xmin+0.5;
	if(ix==6) xmin=3.2;
	double  x = 0.5*(xmax+xmin);
	double dx = 0.5*(xmax-xmin);
	pair<double,double> A0 = calcA0(ix!=6 ? _h_mumu[ix] : _h_nopsi);
	_A_mumu->addPoint(x,A0.first,make_pair(dx,dx),make_pair(A0.second,A0.second) );
      }
      // pion assymetries
      Scatter2DPtr _A_pipi;
      book(_A_pipi,4,1,1);
      for(unsigned int ix=0;ix<_h_pipi.size();++ix) {
	double  x = 0.25+0.1*double(ix);
	pair<double,double> A0 = calcA0(_h_pipi[ix]);
	_A_pipi->addPoint(x,A0.first,make_pair(0.05,0.05),make_pair(A0.second,A0.second) );
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    vector<Profile1DPtr> _h_mumu,_h_pipi;
    Profile1DPtr _h_nopsi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2015_I1388182);

}
