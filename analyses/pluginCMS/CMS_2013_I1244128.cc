// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi and psi(2s) polarization at 7 TeV
  class CMS_2013_I1244128 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2013_I1244128);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // loop over staes
      for(unsigned int ipsi=0;ipsi<2;++ipsi) {
	// rapidity intervals
	for(unsigned int iy=0;iy<2+ipsi;++iy) {
	  // frame defns
	  for(unsigned int iframe=0;iframe<3;++iframe) {
	    // 3 different moments
	    for(unsigned int iobs=0;iobs<3;++iobs) {
	      string name="TMP/POL_"+toString(ipsi)+"_"+toString(iy)
		+"_"+toString(iframe)+"_"+toString(iobs);
	      book(_p_psi[ipsi][iy][iframe][iobs],name,refData(1+24*ipsi,1,1));
	    }
	  }
	}
      }
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,  
                           Particles & mup, Particles & mum) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::MUON ) {
          ++nstable;
    	  mum.push_back(p);
    	}
        else if (id == PID::ANTIMUON) {
          ++nstable;
    	  mup.push_back(p);
        }
        else if (id == PID::PI0 || id == PID::K0S || id == PID::K0L ) {
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, mup, mum);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the beams
      const ParticlePair & beams = apply<Beam>(event, "Beams").beams();
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==443 or Cuts::pid==100443)) {
	// check mu+mu- decay and find muons
      	unsigned int nstable=0;
      	Particles mup,mum;
      	findDecayProducts(p,nstable,mup,mum);
      	if(mup.size()!=1 || mum.size()!=1 || nstable!=2) continue;
      	// pT and rapidity
      	double absrap = p.absrap();
      	double xp = p.perp();
	// state
	unsigned int ipsi = p.pid()/100000;
	// check in fiduical region
	if     (ipsi==0 && (xp<14 || xp>70. || absrap>1.2)) continue;
	else if(ipsi==1 && (xp<14 || xp>50. || absrap>1.5)) continue;
	// rapidity interval
	unsigned int iy = absrap>0.6;
	if(iy==1 && absrap>1.2) iy=2;
	// first the CS frame
	// first boost so upslion momentum =0 in z direction
	Vector3 beta = p.momentum().betaVec();
	beta.setX(0.);beta.setY(0.);
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(beta);
	FourMomentum pp = boost.transform(p.momentum());
	// and then transverse so pT=0
        beta = pp.betaVec();
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(beta);
	// get all the momenta in this frame
	Vector3 muDirn = boost2.transform(boost.transform(mup[0].momentum())).p3().unit();
	FourMomentum p1 = boost2.transform(boost.transform(beams. first.momentum()));
	FourMomentum p2 = boost2.transform(boost.transform(beams.second.momentum()));
	if(beams.first.momentum().z()<0.) swap(p1,p2);
	if(p.rapidity()<0.) swap(p1,p2);
	Vector3 axisy = (p1.p3().cross(p2.p3())).unit();
	Vector3 axisz(0.,0.,1.);
	Vector3 axisx = axisy.cross(axisz);
	double cTheta = axisz.dot(muDirn);
	double cPhi   = axisx.dot(muDirn);
	// fill the moments
	_p_psi[ipsi][iy][0][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_psi[ipsi][iy][0][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_psi[ipsi][iy][0][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
	// now for the HX frame
	beta = p.momentum().betaVec();
	boost = LorentzTransform::mkFrameTransformFromBeta(beta);
	axisz = pp.p3().unit();
	axisx = axisy.cross(axisz);
	cTheta = axisz.dot(muDirn);
	cPhi   = axisx.dot(muDirn);
	// fill the moments
	_p_psi[ipsi][iy][1][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_psi[ipsi][iy][1][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_psi[ipsi][iy][1][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
	// then PX
	axisz = (p1.p3().unit()+p2.p3().unit()).unit();
	axisx = axisy.cross(axisz);
	cTheta = axisz.dot(muDirn);
	cPhi   = axisx.dot(muDirn);
	// fill the moments
	_p_psi[ipsi][iy][2][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_psi[ipsi][iy][2][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_psi[ipsi][iy][2][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Loop over states
      for(unsigned int ipsi=0;ipsi<2;++ipsi) {
	// Loop over rapidity ranges
	for(unsigned int iy=0;iy<2+ipsi;++iy) {
	  // Loop over frame definition
	  for(unsigned int iframe=0;iframe<3;++iframe) {
	    // base no for the ihistos in rivet
	    unsigned int ibase = 24*ipsi+iy+(8+4*ipsi)*iframe;
	    // book scatters
	    Scatter2DPtr lTheta,lPhi,lThetaPhi,lTilde;
	    book(lTheta   ,ibase+1       ,1,1);
	    book(lPhi     ,ibase+3+  ipsi,1,1);
	    book(lThetaPhi,ibase+5+2*ipsi,1,1);
	    book(lTilde   ,ibase+7+3*ipsi,1,1);
	    // histos for the moments
	    Profile1DPtr moment[3];
	    for(unsigned int ix=0;ix<3;++ix)
	      moment[ix] = _p_psi[ipsi][iy][iframe][ix];
	    // loop over bins
	    Scatter2D temphisto(refData(1+24*ipsi, 1, 1));
	    for(unsigned int ibin=0;ibin<moment[0]->bins().size();++ibin) {
	      // extract moments and errors
	      double val[3],err[3];
	      // m1 = lTheta/(3+lTheta), m2 = lPhi/(3+lTheta), m3 = lThetaPhi/(3+lTheta)
	      for(unsigned int ix=0;ix<3;++ix) {
		val[ix] = moment[ix]->bins()[ibin].numEntries()>0 && moment[ix]->bins()[ibin].effNumEntries()>0 ? moment[ix]->bins()[ibin].mean()   : 0.;
		err[ix] = moment[ix]->bins()[ibin].numEntries()>1 && moment[ix]->bins()[ibin].effNumEntries()>1 ? moment[ix]->bins()[ibin].stdErr() : 0.;
	      }
	      // values of the lambdas and their errors
	      double l1 = 3.*val[0]/(1.-val[0]);
	      double e1 = 3./sqr(1.-val[0])*err[0];
	      double l2 = (3.+l1)*val[1];
	      double e2 = 3./sqr(1.-val[0])*sqrt(sqr(err[0]*val[1])+sqr(err[1]*(1.-val[0])));
	      double l3 = (3.+l1)*val[2];
	      double e3 = 3./sqr(1.-val[0])*sqrt(sqr(err[0]*val[1])+sqr(err[1]*(1.-val[0])));
	      double lt = (l1+3.*l2)/(1.-l2);
	      double et = 3./sqr(1.-val[0]-3*val[1])*sqrt(sqr(err[0])+9.*sqr(err[1]));
	      // fill the scatters
	      double              x  = temphisto.point(ibin).x();
	      pair<double,double> ex = temphisto.point(ibin).xErrs();
	      lTheta   ->addPoint(x,l1,ex,make_pair(e1,e1));
	      lPhi     ->addPoint(x,l2,ex,make_pair(e2,e2));
	      lThetaPhi->addPoint(x,l3,ex,make_pair(e3,e3));
	      lTilde   ->addPoint(x,lt,ex,make_pair(et,et));
	    }
	  }
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Profile1DPtr _p_psi[2][3][3][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2013_I1244128);

}
