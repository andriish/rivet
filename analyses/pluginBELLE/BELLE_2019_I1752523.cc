// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief azimuthal asymmetries in pi+/- (pi0, eta, pi+-) production
  class BELLE_2019_I1752523 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2019_I1752523);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs,"FS");
      declare(UnstableParticles(Cuts::pid==111 || Cuts::pid==221), "UFS");
      declare(Thrust(fs),"Thrust");
      declare(Beam(), "Beams");
      // histograms
      unsigned int nbin=20;
      // charged
      for(unsigned int ix=0;ix<3;++ix) {
      	// charged in pT bins
	for(unsigned int iy=0;iy<4;++iy) {
	  book(_h_charged_pt[ix][iy],"TMP/h_charged_pt_"+toString(ix+1)+"_"+toString(iy+1),nbin,0.,M_PI);
	  for(unsigned int iz=0;iz<4;++iz)
	    book(_b_charged_pt[ix][iy][iz],"TMP/b_charged_pt_"+toString(ix)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	}
	// charged in z bins
	for(unsigned int iy=0;iy<6;++iy) {
	  book(_h_charged_z[ix][iy],"TMP/h_charged_z_"+toString(ix+1)+"_"+toString(iy+1),nbin,0.,M_PI);
	  if(iy==5) continue;
	  for(unsigned int iz=0;iz<5;++iz) {
	    book(_b_charged_z[ix][iy][iz],"TMP/b_charged_z_"+toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	  }
	}
      }
      for(unsigned int ix=0;ix<2;++ix) {
	// pi0
	for(unsigned int iy=0;iy<4;++iy) {
	  book(_h_pi0_pt [ix][iy]     ,"TMP/h_pi0_pt_"      +toString(ix+1)+"_"+toString(iy+1),nbin,0.,M_PI);
	  book(_h_eta_pt [ix][iy]     ,"TMP/h_eta_pt_"      +toString(ix+1)+"_"+toString(iy+1),nbin,0.,M_PI);
	  book(_h_pi0B_pt[ix][iy]     ,"TMP/h_pi0B_pt_"     +toString(ix+1)+"_"+toString(iy+1),nbin,0.,M_PI);
	  for(unsigned int iz=0;iz<4;++iz) {
	    book(_b_pi0_pt[ix][iy][iz]      , "TMP/b_pi0_pt_"      +toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	    book(_b_eta_pt[ix][iy][iz]      , "TMP/b_eta_pt_"      +toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	    book(_b_pi0B_pt[ix][iy][iz]     , "TMP/b_pi0B_pt_"     +toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	    book(_b_charged_z_pt[ix][iy][iz], "TMP/b_charged_z_pt_"+toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	    book(_b_pi0_z_pt[ix][iy][iz]    , "TMP/b_pi0_z_pt_"    +toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	    if(iy>2) continue;
	    book(_b_eta_z_pt[ix][iy][iz]    , "TMP/b_eta_z_pt_"    +toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	    book(_b_pi0B_z_pt[ix][iy][iz]   , "TMP/b_pi0B_z_pt_"   +toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	  }
	}
	// pi0 in z bins
	for(unsigned int iy=0;iy<6;++iy) {
	  book(_h_pi0_z[ix][iy],"TMP/h_pi0_z_"+toString(ix+1)+"_"+toString(iy+1),nbin,0.,M_PI);
	  if(iy==5) continue;
	  for(unsigned int iz=0;iz<5;++iz) {
	    book(_b_pi0_z[ix][iy][iz],"TMP/b_pi0_z_"+toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	  }
	}
	// eta in z bins
	for(unsigned int iy=0;iy<5;++iy) {
	  book(_h_eta_z [ix][iy],"TMP/h_eta_z_" +toString(ix+1)+"_"+toString(iy+1),nbin,0.,M_PI);
	  book(_h_pi0B_z[ix][iy],"TMP/h_pi0B_z_"+toString(ix+1)+"_"+toString(iy+1),nbin,0.,M_PI);
	  if(iy>2) continue;
	  for(unsigned int iz=0;iz<3;++iz) {
	    book(_b_eta_z [ix][iy][iz], "TMP/b_eta_z_" +toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	    book(_b_pi0B_z[ix][iy][iz], "TMP/b_pi0B_z_"+toString(ix+1)+"_"+toString(iy+1)+"_"+toString(iz+1),nbin,0.,M_PI);
	  }
	}
      }
    }

    unsigned int iBin_z0(double z) {
      if     (z<.3) return 0;
      else if(z<.4) return 1;
      else if(z<.5) return 2;
      else if(z<.6) return 3;
      else if(z<.7) return 4;
      else          return 5;
    }
    
    unsigned int iBin_z1(double z) {
      if     (z<.2) return 0;
      else if(z<.3) return 1;
      else if(z<.5) return 2;
      else if(z<.7) return 3;
      else          return 4;
    }
    
    unsigned int iBin_pT1(double pT) {
      if     (pT<.15) return 0;
      else if(pT<.3 ) return 1;
      else if(pT<.5 ) return 2;
      else            return 3;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const double c03 = cos(0.3);
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis1;
      if(beams.first.pid()>0)
	axis1 = beams.second.momentum().p3().unit();
      else
	axis1 = beams.first .momentum().p3().unit();
      // apply thrust cuts  T > 0.8 
      Thrust thrust = apply<Thrust>(event,"Thrust");
      if(thrust.thrust()<=0.8) vetoEvent;
      ThreeVector t_z = thrust.thrustAxis();
      ThreeVector t_x = (axis1-t_z.dot(axis1)*t_z).unit();
      ThreeVector t_y = t_z.cross(t_x);
      double cThetaT = t_z.dot(axis1);
      static double cTheta[2]={cos(1.34),cos(2.03)};
      if(cThetaT<cTheta[1] || cThetaT>cTheta[0]) vetoEvent;
      Vector3 cross = axis1.cross(t_z);
      Particles charged = apply<FinalState>(event,"FS").particles(Cuts::abspid==PID::PIPLUS);
      Particles neutral = apply<UnstableParticles>(event,"UFS").particles();
      // first charged particle
      for(unsigned int ix=0;ix<charged.size();++ix) {
	// z and angle cut
	const double x1=2.*charged[ix].momentum().t()/sqrtS();
	double dot1 = t_z.dot(charged[ix].p3().unit());
	if(abs(dot1)<c03 || x1<0.1) continue;
	// second charged particle
	for(unsigned int iy=ix+1;iy<charged.size();++iy) {
	  const double x2=2.*charged[iy].momentum().t()/sqrtS();
	  double dot2 = t_z.dot(charged[iy].p3().unit());
	  if(abs(dot2)<c03 || x2<0.1 || dot1*dot2>0) continue;
	  // randomly order the particles
	  Particle p1=charged[ix], p2=charged[iy];
	  double z1(x1),z2(x2);
	  if(rand()/static_cast<double>(RAND_MAX) < 0.5 ) {
	    swap(p1,p2);
	    swap(z1,z2);
	  }
	  const double pT1 = sqrt(sqr(dot(p1.momentum().p3(), thrust.thrustMajorAxis())) +
				  sqr(dot(p1.momentum().p3(), thrust.thrustMinorAxis())));
	  const double pT2 = sqrt(sqr(dot(p2.momentum().p3(), thrust.thrustMajorAxis())) +
				  sqr(dot(p2.momentum().p3(), thrust.thrustMinorAxis())));
	  if(pT1>3. || pT2>3.) continue;
	  double phi12 = atan2(p1.p3().dot(t_y),p1.p3().dot(t_x))+atan2(p2.p3().dot(t_y),p2.p3().dot(t_x));
	  if(phi12>M_PI)  phi12 -= 2*M_PI;
	  if(phi12<-M_PI) phi12 += 2*M_PI;
	  if(phi12<0.) phi12 = -phi12;
	  if(z1>.2&&z2>.2) {
	    unsigned int ibin = iBin_z0(z1);
	    unsigned int ipT1 = iBin_pT1(pT1);
	    unsigned int ipT2 = iBin_pT1(pT2);
	    unsigned int ibin1 = iBin_z1(z1)-1;
	    if(p1.pid()*p2.pid()>0) {
	      _h_charged_z [0][ibin]->fill(phi12);
	      _h_charged_pt[0][ipT1]->fill(phi12);
	      _b_charged_pt[0][ipT1][ipT2]->fill(phi12);
	      _h_pi0_z[1][ibin]->fill(phi12);
	      _h_pi0_pt[1][ipT1]->fill(phi12);
	      _b_pi0_pt[1][ipT1][ipT2]->fill(phi12);
	      _b_pi0_z_pt[1][ibin1][ipT1]->fill(phi12);
	      if(z1>.3 && z2>.3) {
		_h_eta_z [1][ibin-1]->fill(phi12);
		_h_pi0B_z[1][ibin-1]->fill(phi12);
		_h_eta_pt [1][ipT1]->fill(phi12);
		_h_pi0B_pt[1][ipT1]->fill(phi12);
		_b_eta_pt [1][ipT1][ipT2]->fill(phi12);
		_b_pi0B_pt[1][ipT1][ipT2]->fill(phi12);
		_b_eta_z_pt [1][ibin1-1][ipT1]->fill(phi12);
		_b_pi0B_z_pt[1][ibin1-1][ipT1]->fill(phi12);
	      }
	    }
	    else {
	      _h_charged_z [1][ibin]->fill(phi12);
	      _h_charged_pt[1][ipT1]->fill(phi12);
	      _b_charged_pt[1][ipT1][ipT2]->fill(phi12);
	      _b_charged_z_pt[0][ibin1][ipT1]->fill(phi12);
	    }
	    _h_charged_z [2][ibin]->fill(phi12);
	    _h_charged_pt[2][ipT1]->fill(phi12);
	    _b_charged_pt[2][ipT1][ipT2]->fill(phi12);
	    _b_charged_z_pt[1][ibin1][ipT1]->fill(phi12);
	  }
	  unsigned int ibin1 = iBin_z1(z1), ibin2 = iBin_z1(z2);
	  if(p1.pid()*p2.pid()>0) {
	    _b_charged_z[0][ibin1][ibin2]->fill(phi12); 
	    _b_pi0_z[1][ibin1][ibin2]->fill(phi12);
	    if(z1>.3 && z2>.3) {
	      _b_eta_z [1][ibin1-2][ibin2-2]->fill(phi12);
	      _b_pi0B_z[1][ibin1-2][ibin2-2]->fill(phi12);
	    }
	  }
	  else
	    _b_charged_z[1][ibin1][ibin2]->fill(phi12);
	  _b_charged_z[2][ibin1][ibin2]->fill(phi12);
	}
	// pi0 or eta
	for(unsigned int iy=0;iy<neutral.size();++iy) {
	  const double x2=2.*neutral[iy].momentum().t()/sqrtS();
	  double dot2 = t_z.dot(neutral[iy].p3().unit());
	  if(abs(dot2)<c03 || x2<0.2 || dot1*dot2>0) continue;
	  // neutral particle first
	  Particle p2=charged[ix], p1=neutral[iy];
	  double z1(x2),z2(x1);
	  const double pT1 = sqrt(sqr(dot(p1.momentum().p3(), thrust.thrustMajorAxis())) +
	  			  sqr(dot(p1.momentum().p3(), thrust.thrustMinorAxis())));
	  const double pT2 = sqrt(sqr(dot(p2.momentum().p3(), thrust.thrustMajorAxis())) +
	  			  sqr(dot(p2.momentum().p3(), thrust.thrustMinorAxis())));
	  if(pT1>3. || pT2>3.) continue;
	  double phi12 = atan2(p1.p3().dot(t_y),p1.p3().dot(t_x))+atan2(p2.p3().dot(t_y),p2.p3().dot(t_x));
	  if(phi12>M_PI)  phi12 -= 2*M_PI;
	  if(phi12<-M_PI) phi12 += 2*M_PI;
	  if(phi12<0.) phi12 = -phi12;
	  if(z1>.2) {
	    unsigned int ibin = iBin_z0(z1);
	    unsigned int ipT1 = iBin_pT1(pT1);
	    unsigned int ipT2 = iBin_pT1(pT2);
	    unsigned int ibin1 = iBin_z1(z1)-1;
	    if(p1.pid()==111) {
	      _h_pi0_z[0][ibin]->fill(phi12);
	      _h_pi0_pt[0][ipT1]->fill(phi12);
	      _b_pi0_pt[0][ipT1][ipT2]->fill(phi12);
	      _b_pi0_z_pt[0][ibin1][ipT1]->fill(phi12);
	      if(z1>.3 && z2>.3) {
		_h_pi0B_z[0][ibin-1]->fill(phi12);
		_h_pi0B_pt[0][ipT1]->fill(phi12);
		_b_pi0B_pt[0][ipT1][ipT2]->fill(phi12);
		_b_pi0B_z_pt[0][ibin1-1][ipT1]->fill(phi12);
	      }
	    }
	    else if(z1>.3 && z2>.3) {
	      _h_eta_z [0][ibin-1]->fill(phi12);
	      _h_eta_pt [0][ipT1]->fill(phi12);
	      _b_eta_pt [0][ipT1][ipT2]->fill(phi12);
	      _b_eta_z_pt [0][ibin1-1][ipT1]->fill(phi12);
	    }
	  }
	  unsigned int ibin1 = iBin_z1(z1), ibin2 = iBin_z1(z2);
	  if(p1.pid()==111) {
	    _b_pi0_z[0][ibin1][ibin2]->fill(phi12);
	    if(z1>.3 && z2>.3) _b_pi0B_z[0][ibin1-2][ibin2-2]->fill(phi12);
	  }
	  else if(z1>.3 && z2>.3) {
	    _b_eta_z [0][ibin1-2][ibin2-2]->fill(phi12);
	  }
	}
      }
    }

    pair<double,double> calcAsymmetry(Scatter2DPtr hist,double fact=1.) {
      double sum1(0.),sum2(0.);
      for (auto point : hist->points() ) {
	double Oi = point.y();
	if(Oi==0. || std::isnan(Oi) ) continue;
	double ai = 1.;
	double bi = (sin(fact*point.xMax())-sin(fact*point.xMin()))/(point.xMax()-point.xMin())/fact;
	double Ei = point.yErrAvg();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      if(sum1==0.) return make_pair(0.,0.);
      return make_pair(sum2/sum1*1e2,sqrt(1./sum1)*1e2);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // charged in pT bins
      Scatter2D href(refData(1,1,1));
      Scatter2DPtr h_charged_UL,h_charged_UC;
      book(h_charged_UL,1,1,1);
      book(h_charged_UC,1,1,2);
      Scatter2DPtr h_pi0_UL,h_eta_UL,h_pi0B_UL;
      book(h_pi0_UL,6,1,1);
      book(h_eta_UL,7,1,1);
      book(h_pi0B_UL,8,1,1);
      for(unsigned int ix=0;ix<4;++ix) {
	// first the 1d dists
	for(unsigned int ii=0;ii<3;++ii) {
	  normalize(_h_charged_pt[ii][ix]);
	  if(ii==2) continue;
      	  normalize(_h_pi0_pt [ii][ix]);
      	  normalize(_h_eta_pt [ii][ix]);
      	  normalize(_h_pi0B_pt[ii][ix]);
	}
	Scatter2DPtr htemp;
	// UL
	// charged
	book(htemp,"TMP/R_charged_pt_UL_"+toString(ix));
	divide(_h_charged_pt[1][ix],_h_charged_pt[0][ix],htemp);
	pair<double,double> asym = calcAsymmetry(htemp);
	h_charged_UL->addPoint(href.points()[ix].x(),asym.first,href.points()[ix].xErrs(),make_pair(asym.second,asym.second));
	// UC 
	book(htemp,"TMP/R_charged_pt_UC_"+toString(ix));
	divide(_h_charged_pt[1][ix],_h_charged_pt[2][ix],htemp);
	asym = calcAsymmetry(htemp);
	h_charged_UC->addPoint(href.points()[ix].x(),asym.first,href.points()[ix].xErrs(),make_pair(asym.second,asym.second));
	// pi0
	book(htemp,"TMP/R_pi0_pt_UL_"+toString(ix));
	divide(_h_pi0_pt[0][ix],_h_pi0_pt[1][ix],htemp);
	asym = calcAsymmetry(htemp);
	h_pi0_UL->addPoint(href.points()[ix].x(),asym.first,href.points()[ix].xErrs(),make_pair(asym.second,asym.second));
	// eta
	book(htemp,"TMP/R_eta_pt_UL_"+toString(ix));
	divide(_h_eta_pt[0][ix],_h_eta_pt[1][ix],htemp);
	asym = calcAsymmetry(htemp);
	h_eta_UL->addPoint(href.points()[ix].x(),asym.first,href.points()[ix].xErrs(),make_pair(asym.second,asym.second));
	// pi0 z>.3
	book(htemp,"TMP/R_pi0B_pt_UL_"+toString(ix));
	divide(_h_pi0B_pt[0][ix],_h_pi0B_pt[1][ix],htemp);
	asym = calcAsymmetry(htemp);
	h_pi0B_UL->addPoint(href.points()[ix].x(),asym.first,href.points()[ix].xErrs(),make_pair(asym.second,asym.second));
	// then 2nd dists
	Scatter2DPtr h_charged_UL2,h_charged_UC2;
	book(h_charged_UL2,1,2+ix,1);
	book(h_charged_UC2,1,2+ix,2);
	Scatter2DPtr h_pi0_UL2,h_eta_UL2,h_pi0B_UL2;
	book(h_pi0_UL2,6,2+ix,1);
	book(h_eta_UL2,7,2+ix,1);
	book(h_pi0B_UL2,8,2+ix,1);
	for(unsigned int iy=0;iy<4;++iy) {
	  for(unsigned int ii=0;ii<3;++ii) {
	    normalize(_b_charged_pt[ii][ix][iy]);
	    if(ii==2) continue;
       	    normalize(_b_pi0_pt [ii][ix][iy]);
       	    normalize(_b_eta_pt [ii][ix][iy]);
       	    normalize(_b_pi0B_pt[ii][ix][iy]);
	  }
	  // UL 
	  book(htemp,"TMP/R_charged_pt_UL_"+toString(ix)+"_"+toString(iy));
	  divide(_b_charged_pt[1][ix][iy],_b_charged_pt[0][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_charged_UL2->addPoint(href.points()[iy].x(),asym.first,href.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	  // UC 
	  book(htemp,"TMP/R_charged_pt_UC_"+toString(ix)+"_"+toString(iy));
	  divide(_b_charged_pt[1][ix][iy],_b_charged_pt[2][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_charged_UC2->addPoint(href.points()[iy].x(),asym.first,href.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	  // pi0
	  book(htemp,"TMP/R_pi0_pt_UL_"+toString(ix)+"_"+toString(iy));
	  divide(_b_pi0_pt[0][ix][iy],_b_pi0_pt[1][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_pi0_UL2->addPoint(href.points()[iy].x(),asym.first,href.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	  // eta
	  book(htemp,"TMP/R_eta_pt_UL_"+toString(ix)+"_"+toString(iy));
	  divide(_b_eta_pt[0][ix][iy],_b_eta_pt[1][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_eta_UL2->addPoint(href.points()[iy].x(),asym.first,href.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	  // pi0 z>0.3
	  book(htemp,"TMP/R_pi0B_pt_UL_"+toString(ix)+"_"+toString(iy));
	  divide(_b_pi0B_pt[0][ix][iy],_b_pi0B_pt[1][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_pi0B_UL2->addPoint(href.points()[iy].x(),asym.first,href.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	}
      }
      // charged and pi0 in z bins
      href = Scatter2D(refData(2,1,1));
      Scatter2D href2(refData(2,2,1));
      book(h_charged_UL,2,1,1);
      book(h_charged_UC,2,1,2);
      book(h_pi0_UL,3,1,1);
      book(h_eta_UL,4,1,1);
      book(h_pi0B_UL,5,1,1);
      for(unsigned int ix=0;ix<6;++ix) {
      	// first the 1d dists
      	for(unsigned int ii=0;ii<3;++ii) {
	  normalize(_h_charged_z[ii][ix]);
       	  if(ii==2) continue;
       	  normalize(_h_pi0_z [ii][ix]);
	  if(ix==5) continue;
       	  normalize(_h_eta_z [ii][ix]);
       	  normalize(_h_pi0B_z[ii][ix]);
	}
	Scatter2DPtr htemp;
	// UL
	// charged
	book(htemp,"TMP/R_charged_z_UL_"+toString(ix));
	divide(_h_charged_z[1][ix],_h_charged_z[0][ix],htemp);
	pair<double,double> asym = calcAsymmetry(htemp);
	h_charged_UL->addPoint(href.points()[ix].x(),asym.first,href.points()[ix].xErrs(),make_pair(asym.second,asym.second));
	// UC 
	book(htemp,"TMP/R_charged_z_UC_"+toString(ix));
	divide(_h_charged_z[1][ix],_h_charged_z[2][ix],htemp);
	asym = calcAsymmetry(htemp);
	h_charged_UC->addPoint(href.points()[ix].x(),asym.first,href.points()[ix].xErrs(),make_pair(asym.second,asym.second));
	// pi0
	book(htemp,"TMP/R_pi0_z_UL_"+toString(ix));
	divide(_h_pi0_z[0][ix],_h_pi0_z[1][ix],htemp);
	asym = calcAsymmetry(htemp);
	h_pi0_UL->addPoint(href.points()[ix].x(),asym.first,href.points()[ix].xErrs(),make_pair(asym.second,asym.second));
	if(ix==5) continue;
	// eta
	book(htemp,"TMP/R_eta_z_UL_"+toString(ix));
	divide(_h_eta_z[0][ix],_h_eta_z[1][ix],htemp);
	asym = calcAsymmetry(htemp);
	h_eta_UL->addPoint(href.points()[ix+1].x(),asym.first,href.points()[ix+1].xErrs(),make_pair(asym.second,asym.second));
	// pi0 z>.3
	book(htemp,"TMP/R_pi0B_z_UL_"+toString(ix));
	divide(_h_pi0B_z[0][ix],_h_pi0B_z[1][ix],htemp);
	asym = calcAsymmetry(htemp);
	h_pi0B_UL->addPoint(href.points()[ix+1].x(),asym.first,href.points()[ix+1].xErrs(),make_pair(asym.second,asym.second));
	// then 2nd dists
	Scatter2DPtr h_charged_UL2,h_charged_UC2;
	book(h_charged_UL2,2,2+ix,1);
	book(h_charged_UC2,2,2+ix,2);
      	Scatter2DPtr h_pi0_UL2,h_eta_UL2,h_pi0B_UL2;
      	book(h_pi0_UL2,3,2+ix,1);
	if(ix<3) {
	  book(h_eta_UL2,4,2+ix,1);
	  book(h_pi0B_UL2,5,2+ix,1);
	}
      	for(unsigned int iy=0;iy<5;++iy) {
      	  for(unsigned int ii=0;ii<3;++ii) {
      	    normalize(_b_charged_z[ii][ix][iy]);
       	    if(ii==2) continue;
	    normalize(_b_pi0_z [ii][ix][iy]);
	    if(ix>2||iy>2) continue;
	    normalize(_b_eta_z [ii][ix][iy]);
	    normalize(_b_pi0B_z[ii][ix][iy]);
       	  }
	  // UL 
	  book(htemp,"TMP/R_charged_z_UL_"+toString(ix)+"_"+toString(iy));
	  divide(_b_charged_z[1][ix][iy],_b_charged_z[0][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_charged_UL2->addPoint(href2.points()[iy].x(),asym.first,href2.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	  // UC 
	  book(htemp,"TMP/R_charged_z_UC_"+toString(ix)+"_"+toString(iy));
	  divide(_b_charged_z[1][ix][iy],_b_charged_z[2][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_charged_UC2->addPoint(href2.points()[iy].x(),asym.first,href2.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	  // pi0
	  book(htemp,"TMP/R_pi0_z_UL_"+toString(ix)+"_"+toString(iy));
	  divide(_b_pi0_z[0][ix][iy],_b_pi0_z[1][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_pi0_UL2->addPoint(href2.points()[iy].x(),asym.first,href2.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	  if(ix>2||iy>2) continue;
	  // eta
	  book(htemp,"TMP/R_eta_z_UL_"+toString(ix)+"_"+toString(iy));
	  divide(_b_eta_z[0][ix][iy],_b_eta_z[1][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_eta_UL2->addPoint(href2.points()[iy+2].x(),asym.first,href2.points()[iy+2].xErrs(),make_pair(asym.second,asym.second));
	  // pi0 z>0.3
	  book(htemp,"TMP/R_pi0B_z_UL_"+toString(ix)+"_"+toString(iy));
	  divide(_b_pi0B_z[0][ix][iy],_b_pi0B_z[1][ix][iy],htemp);
	  asym = calcAsymmetry(htemp);
	  h_pi0B_UL2->addPoint(href2.points()[iy+2].x(),asym.first,href2.points()[iy+2].xErrs(),make_pair(asym.second,asym.second));
       	}
      }
      // finally z and pT
      href = Scatter2D(refData(9,2,1));
      // only 2d dists
      for(unsigned int ix=0;ix<4;++ix) {
	Scatter2DPtr h_charged_UC2;
	book(h_charged_UC2,9,1+ix,1);
	Scatter2DPtr h_pi0_UL2,h_eta_UL2,h_pi0B_UL2;
	book(h_pi0_UL2,10,1+ix,1);
	if(ix<3) {
	  book(h_eta_UL2,11,1+ix,1);
	  book(h_pi0B_UL2,12,1+ix,1);
	}
      	for(unsigned int iy=0;iy<4;++iy) {
      	  for(unsigned int ii=0;ii<2;++ii) {
	    normalize(_b_charged_z_pt[ii][ix][iy]);
	    normalize(_b_pi0_z_pt [ii][ix][iy]);
	    if(ix>2) continue;
      	    normalize(_b_eta_z_pt [ii][ix][iy]);
      	    normalize(_b_pi0B_z_pt[ii][ix][iy]);
       	  }
	  pair<double,double> asym;
	  Scatter2DPtr htemp;
      	  // UC (missing last point for first histo)
	  if(ix!=0 || iy!=3) {
	    book(htemp,"TMP/R_charged_z_pt_UC_"+toString(ix)+"_"+toString(iy));
	    divide(_b_charged_z_pt[1][ix][iy],_b_charged_z_pt[2][ix][iy],htemp);
	    asym = calcAsymmetry(htemp);
	    h_charged_UC2->addPoint(href.points()[iy].x(),asym.first,href.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	    // pi0
	    book(htemp,"TMP/R_pi0_z_pt_UL_"+toString(ix)+"_"+toString(iy));
	    divide(_b_pi0_z_pt[0][ix][iy],_b_pi0_z_pt[1][ix][iy],htemp);
	    asym = calcAsymmetry(htemp);
	    h_pi0_UL2->addPoint(href.points()[iy].x(),asym.first,href.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	  }
	  if(ix>2) continue;
      	  // eta
      	  book(htemp,"TMP/R_eta_z_pt_UL_"+toString(ix)+"_"+toString(iy));
      	  divide(_b_eta_z_pt[0][ix][iy],_b_eta_z_pt[1][ix][iy],htemp);
      	  asym = calcAsymmetry(htemp);
      	  h_eta_UL2->addPoint(href.points()[iy].x(),asym.first,href.points()[iy].xErrs(),make_pair(asym.second,asym.second));
      	  // pi0 z>0.3
      	  book(htemp,"TMP/R_pi0B_z_pt_UL_"+toString(ix)+"_"+toString(iy));
      	  divide(_b_pi0B_z_pt[0][ix][iy],_b_pi0B_z_pt[1][ix][iy],htemp);
      	  asym = calcAsymmetry(htemp);
      	  h_pi0B_UL2->addPoint(href.points()[iy].x(),asym.first,href.points()[iy].xErrs(),make_pair(asym.second,asym.second));
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_charged_pt[3][4],_h_charged_z[3][6];
    Histo1DPtr _b_charged_pt[3][4][4],_b_charged_z[3][5][5],_b_charged_z_pt[2][4][4];
    Histo1DPtr _h_pi0_z[2][6],_h_eta_z[2][5],_h_pi0B_z[2][5],_h_pi0_pt[2][4],_h_eta_pt[2][4],_h_pi0B_pt[2][4];
    Histo1DPtr _b_pi0_z[2][5][5],_b_eta_z[2][3][3],_b_pi0B_z[2][3][3],_b_pi0_pt[2][4][4];
    Histo1DPtr _b_eta_pt[2][4][4],_b_pi0B_pt[2][4][4],_b_pi0_z_pt[2][4][4],_b_pi0B_z_pt[2][3][4],_b_eta_z_pt[2][3][4];
    /// @}
  };


  RIVET_DECLARE_PLUGIN(BELLE_2019_I1752523);

}
