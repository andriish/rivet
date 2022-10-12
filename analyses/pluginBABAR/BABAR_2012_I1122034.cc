// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief gamma gamma -> X(3915) -> J/psi omega
  class BABAR_2012_I1122034 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1122034);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(),"FS");
      declare(UnstableParticles(Cuts::pid==223 or Cuts::pid==443), "UFS");
      // histograms
      book(_h_mass,1,1,1);
      for(unsigned int ix=0;ix<4;++ix) {
	book(_h_angle1[ix],2,1,1+ix);
	if(ix<3) book(_h_angle2[ix],3,1,1+ix);
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        } else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    bool findScattered(Particle beam, double& q2) {
      bool found = false;
      Particle scat = beam;
      while (!scat.children().empty()) {
        found = false;
        for (const Particle & p : scat.children()) {
          if (p.pid()==scat.pid()) {
            scat=p;
            found=true;
            break;
          }
        }
        if (!found) break;
      }
      if (!found) return false;
      q2 = -(beam.momentum() - scat.momentum()).mass2();
      return true;
    }
    
    void findChildren(const Particle & p, Particles & pim, Particles & pip,
		      Particles & pi0, unsigned int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.pid()==PID::PIPLUS) {
	  pip.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PIMINUS) {
	  pim.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PI0) {
	  pi0.push_back(child);
	  ncount+=1;
	}
	else if(child.children().empty()) {
	  ncount+=1;
	}
    	else
    	  findChildren(child,pim,pip,pi0,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find scattered leptons and calc Q2
      const Beam& beams = apply<Beam>(event, "Beams");
      double q12 = -1, q22 = -1;
      if (!findScattered(beams.beams().first,  q12)) vetoEvent;
      if (!findScattered(beams.beams().second, q22)) vetoEvent;
      // check the final state
      const FinalState & fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // find the J/psi
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      Particle omega,psi;
      bool found=false;
      for (const Particle &p1 : ufs.particles(Cuts::pid==443)) {
        if(p1.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ntotal);
	for(const Particle & p2 : ufs.particles(Cuts::pid==223)) {
	  if(p2.children().empty()) continue;
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  found = true;
	  for(auto const & val : nRes2) {
	    if(abs(val.first)==11) {
	      if(val.second!=1) {
		found = false;
		break;
	      }
	    }
	    else if(val.second!=0) {
	      found = false;
	      break;
	    }
	  }
	  if (found) {
	    psi   = p1;
	    omega = p2;
	    break;
	  }
	}
      }
      if(!found) vetoEvent;
      FourMomentum psum = omega.momentum()+psi.momentum();
      if(psum.pT()>0.2) vetoEvent;
      // mass distribution
      _h_mass->fill(psum.mass());
      // from now on we need specific decay modes of J/psi and omega
      // first J/psi -> l+l-
      if(psi.children().size()!=2) vetoEvent;
      if(psi.children()[0].pid()!=-psi.children()[1].pid()) vetoEvent;
      if(psi.children()[0].abspid()!=11 && psi.children()[0].abspid()!=13) vetoEvent;
      Particle ep = psi.children()[0];
      Particle em = psi.children()[1];
      if(ep.pid()>0) swap(ep,em);
      // omega decay
      unsigned int ncount=0;
      Particles pip,pim,pi0;
      findChildren(omega,pim,pip,pi0,ncount);
      if( ncount!=3 || !(pim.size()==1 && pip.size()==1 && pi0.size()==1)) vetoEvent;
      // boost to gamma gamma frame
      LorentzTransform boostCMS = LorentzTransform::mkFrameTransformFromBeta(psum.betaVec());
      FourMomentum pPsi   = boostCMS.transform(psi   .momentum());
      FourMomentum pOmega = boostCMS.transform(omega .momentum());
      FourMomentum pLp    = boostCMS.transform(ep    .momentum());
      FourMomentum pPip   = boostCMS.transform(pip[0].momentum());
      FourMomentum pPim   = boostCMS.transform(pim[0].momentum());
      Vector3 axis(0.,0.,1.);
      // lepton angle
      double cosL = pLp.p3().unit().dot(axis);
      _h_angle1[0]->fill(cosL);
      LorentzTransform boostOmega = LorentzTransform::mkFrameTransformFromBeta(pOmega.betaVec());
      pPip   = boostOmega.transform(pPip);
      pPim   = boostOmega.transform(pPim);
      // omega decay plane normal angle
      Vector3 axisOmega = pPip.p3().cross(pPim.p3()).unit();
      double cosN = axisOmega.dot(axis);
      _h_angle1[1]->fill(cosN);
      // angle lepton and omega
      _h_angle1[2]->fill(pLp.p3().unit().dot(axisOmega));
      // helicity angle
      _h_angle1[3]->fill(pPsi.p3().unit().dot(psum.p3().unit()));
      // now the new frame
      Vector3 axisZ = pOmega.p3().unit();
      Vector3 axisY = axisZ.cross(axisOmega);
      Vector3 axisX = axisY.cross(axisZ);
      // second set of angles
      LorentzTransform boostPsi = LorentzTransform::mkFrameTransformFromBeta(pPsi.betaVec());
      Vector3 axisL = boostPsi.transform(pLp).p3().unit();
      _h_angle2[0]->fill(axisZ.dot(axisOmega));
      _h_angle2[1]->fill(axisL.dot(pPsi.p3().unit()));
      axisZ *=-1.;
      axisX *=-1.;
      Vector3 axisnp = axisL.cross(axisZ).unit();
      double phiL = atan2(axisnp.dot(axisY),axisnp.dot(axisX));
      double phiN = atan2(axisOmega.dot(axisY),axisOmega.dot(axisX));
      _h_angle2[2]->fill(mapAngleMPiToPi(phiL-phiN)/M_PI*180.);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass,1.,false);
      for(unsigned int ix=0;ix<4;++ix) {
	normalize(_h_angle1[ix],1.,false);
	if(ix<3) normalize(_h_angle2[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass,_h_angle1[4],_h_angle2[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1122034);

}
