// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Upsilon polarization at 7 and 8 TeV
  class LHCB_2017_I1621596 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2017_I1621596);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      int iloc=-1;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 0;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 1;
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 7 or 8 TeV.");
      // histograms
      _ybins={2.2,3.0,3.5,4.5};
      for(unsigned int iups=0;iups<3;++iups) {
	for(unsigned int iframe=0;iframe<3;++iframe) {
	  for(unsigned int imom=0;imom<3;++imom) {
	    for(unsigned int iy=0;iy<3;++iy) {
	      book(_p_Upsilon[iups][iframe][iy][imom],
		   "TMP/UPS_"+toString(iups)+"_"+toString(iframe)+"_"+toString(iy)+"_"+toString(imom),
		   refData(32*iups+4*iloc+8*iframe+1,1,iy+1));
	    }
	    book(_p_Upsilon[iups][iframe][3][imom],
		 "TMP/UPS_"+toString(iups)+"_"+toString(iframe)+"_3_"+toString(imom),
		 refData(32*iups+4*iloc+25,1,iframe+1));
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
      for (const Particle& p : ufs.particles(Cuts::pid==553 || Cuts::pid==100553 || Cuts::pid==200553)) {
      	// pT and rapidity
      	double rapidity = p.rapidity();
      	double xp = p.perp();
	if(rapidity<2.2 || rapidity >4.5) continue;
	// which upsilon
	unsigned int iups=p.pid()/100000;
	// polarization
      	unsigned int nstable=0;
      	Particles mup,mum;
      	findDecayProducts(p,nstable,mup,mum);
      	if(mup.size()!=1 || mum.size()!=1 || nstable!=2) continue;
	unsigned int iy=0;
	for(iy=0;iy<3;++iy) if(rapidity<_ybins[iy+1]) break;
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
	_p_Upsilon[iups][1][iy][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_Upsilon[iups][1][iy][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_Upsilon[iups][1][iy][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
	_p_Upsilon[iups][1][3 ][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_Upsilon[iups][1][3 ][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_Upsilon[iups][1][3 ][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
	// Gottfried-Jackson frame
	axisz = p1.p3().unit();
	axisx = axisy.cross(axisz);
	cTheta = axisz.dot(muDirn);
	cPhi   = axisx.dot(muDirn);
	// fill the moments
	_p_Upsilon[iups][2][iy][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_Upsilon[iups][2][iy][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_Upsilon[iups][2][iy][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
	_p_Upsilon[iups][2][3 ][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_Upsilon[iups][2][3 ][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_Upsilon[iups][2][3 ][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
	// now for the HX frame
	beta = p.momentum().betaVec();
	boost = LorentzTransform::mkFrameTransformFromBeta(beta);
	axisz = pp.p3().unit();
	axisx = axisy.cross(axisz);
	cTheta = axisz.dot(muDirn);
	cPhi   = axisx.dot(muDirn);
	// fill the moments
	_p_Upsilon[iups][0][iy][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_Upsilon[iups][0][iy][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_Upsilon[iups][0][iy][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
	_p_Upsilon[iups][0][3 ][0]->fill(xp, 1.25*(3.*sqr(cTheta)-1.));
	_p_Upsilon[iups][0][3 ][1]->fill(xp, 1.25*(1.-sqr(cTheta))*(2.*sqr(cPhi)-1.));
	_p_Upsilon[iups][0][3 ][2]->fill(xp, 2.5 *cTheta*sqrt(1.-sqr(cTheta))*cPhi);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      int iloc=-1;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 0;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 1;
      }
      // loop over upslion
      for(unsigned int iups=0;iups<3;++iups) {
	// loop over iframe
	for(unsigned int iframe=0;iframe<3;++iframe) {
	  unsigned int ibase = 32*iups+4*iloc+8*iframe;
	  unsigned int ibase2 = 32*iups+4*iloc+24;
	  // rapidity range
	  for(unsigned int iy=0;iy<4;++iy) {
	    // book scatters
	    Scatter2DPtr lTheta,lPhi,lThetaPhi,lTilde;
	    if(iy<3) {
	      book(lTheta   ,ibase+1,1,1+iy);
	      book(lPhi     ,ibase+3,1,1+iy);
	      book(lThetaPhi,ibase+2,1,1+iy);
	      book(lTilde   ,ibase+4,1,1+iy);
	    }
	    else {
	      book(lTheta   ,ibase2+1,1,1+iframe);
	      book(lPhi     ,ibase2+3,1,1+iframe);
	      book(lThetaPhi,ibase2+2,1,1+iframe);
	      book(lTilde   ,ibase2+4,1,1+iframe);
	    }
	    // histos for the moments
	    Profile1DPtr moment[3];
	    for(unsigned int ix=0;ix<3;++ix)
	      moment[ix] = _p_Upsilon[iups][iframe][iy][ix];
	    // loop over bins
	    Scatter2D temphisto;
	    if(iy<3) {
	      temphisto = refData(ibase+1,1,1+iy);
	    }
	    else {
	      temphisto = refData(ibase2+1,1,1+iframe);
	    }
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
    Profile1DPtr _p_Upsilon[3][3][4][3];
    vector<double>  _ybins;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2017_I1621596);

}
