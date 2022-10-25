// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief azimuthal asymmetries in pipi
  class BABAR_2014_I1254862 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2014_I1254862);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs,"FS");
      declare(Thrust(fs),"Thrust");
      declare(Beam(), "Beams");
      // declare the histos for the distributions
      string charge[3] = {"Like","Opposite","All"};
      unsigned int nbin=20;
      for(unsigned int icharge=0;icharge<3;++icharge) {
	for(unsigned int ibin1=0;ibin1<6;++ibin1) {
	  for(unsigned int ibin2=0;ibin2<6;++ibin2) {
	    book(_h_thrust[icharge][ibin1][ibin2],"/TMP/h_thrust_"+charge[icharge]+"_" +toString(ibin1+1) + "_" + toString(ibin2+1),nbin,0.,M_PI);
	    book(_h_hadron[icharge][ibin1][ibin2],"/TMP/h_hadron_"+charge[icharge]+"_" +toString(ibin1+1) + "_" + toString(ibin2+1),nbin,0.,M_PI);
	  }
	}
      }
    }

    unsigned int iBin(double z) {
      if     (z<.2) return 0;
      else if(z<.3) return 1;
      else if(z<.4) return 2;
      else if(z<.5) return 3;
      else if(z<.7) return 4;
      else          return 5;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis1;
      if(beams.first.pid()>0)
	axis1 = beams.first .momentum().p3().unit();
      else
	axis1 = beams.second.momentum().p3().unit();
      // apply thrust cuts  T > 0.8 
      Thrust thrust = apply<Thrust>(event,"Thrust");
      if(thrust.thrust()<=0.8) vetoEvent;
      // construct x,y,z axes for thrust defn
      ThreeVector t_z = thrust.thrustAxis();
      ThreeVector t_x = (axis1-t_z.dot(axis1)*t_z).unit();
      ThreeVector t_y = t_z.cross(t_x);
      // loop over the particles
      Particles charged = apply<FinalState>(event,"FS").particles(Cuts::abspid==PID::PIPLUS);
      for(unsigned int ix=0;ix<charged.size();++ix) {
	// z and angle cut
	const double x1=2.*charged[ix].momentum().t()/sqrtS();
	if(x1<0.15||x1>.9) continue;
	if(abs(t_z.angle(charged[ix].momentum().p3()))>0.25*M_PI) continue;
	double dot1 = t_z.dot(charged[ix].p3());
	for(unsigned int iy=ix+1;iy<charged.size();++iy) {
	  const double x2=2.*charged[iy].momentum().t()/sqrtS();
	  // z and angle cut
	  if(x2<0.15||x2>.9) continue;
	  if(abs(t_z.angle(charged[ix].momentum().p3()))>0.25*M_PI) continue;
	  // different hemi
	  double dot2 = t_z.dot(charged[iy].p3());
	  if(dot1*dot2>0.) continue;
	  Particle p1=charged[ix], p2=charged[iy];
	  double z1(x1),z2(x2);
	  // randomly order the particles
	  if(rand()/static_cast<double>(RAND_MAX) < 0.5 ) {
	    swap(p1,p2);
	    swap(z1,z2);
	  }
	  // thrust def
	  double phi12 = atan2(p1.p3().dot(t_y),p1.p3().dot(t_x))+atan2(p2.p3().dot(t_y),p2.p3().dot(t_x));
	  if(phi12>M_PI)  phi12 -= 2*M_PI;
	  if(phi12<-M_PI) phi12 += 2*M_PI;
	  if(phi12<0.) phi12 = -phi12;
	  // hadron defn
	  ThreeVector h_z = p2.p3().unit();
	  ThreeVector h_x = (axis1-h_z.dot(axis1)*h_z).unit();
	  ThreeVector pt1 = p1.p3()-h_z.dot(p1.p3())*h_z;
	  double phi0 = pt1.angle(h_x);
	  if(phi0>M_PI)  phi0 -= 2*M_PI;
	  if(phi0<-M_PI) phi0 += 2*M_PI;
	  unsigned int ibin1=iBin(z1);
	  unsigned int ibin2=iBin(z2);
	  if(p1.pid()==p2.pid()) {
	    _h_thrust[0][ibin1][ibin2]->fill(phi12);
	    _h_hadron[0][ibin1][ibin2]->fill(phi0);
	  }
	  else {
	    _h_thrust[1][ibin1][ibin2]->fill(phi12);
	    _h_hadron[1][ibin1][ibin2]->fill(phi0);
	  }
	  _h_thrust[2][ibin1][ibin2]->fill(phi12);
	  _h_hadron[2][ibin1][ibin2]->fill(phi0);
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
      for(unsigned int ibin1=0;ibin1<6;++ibin1) {
	Scatter2D href(refData(1,1+ibin1,1));
	Scatter2DPtr hthrustUL,hhadronUL;
	book(hthrustUL,1,1+ibin1,1);
	book(hhadronUL,1,1+ibin1,2);
	Scatter2DPtr hthrustUC,hhadronUC;
	book(hthrustUC,2,1+ibin1,1);
	book(hhadronUC,2,1+ibin1,2);
	for(unsigned int ibin2=0;ibin2<6;++ibin2) {
	  const Point2D & pRef = href.points()[ibin2];
	  for(unsigned int icharge=0;icharge<3;++icharge) {
	    normalize(_h_thrust[icharge][ibin1][ibin2]);
	    normalize(_h_hadron[icharge][ibin1][ibin2]);
	  }
	  Scatter2DPtr htemp;
	  book(htemp,"TMP/R_thrust_UL_"+toString(ibin1)+"_"+toString(ibin2));
	  // UL thrust
	  divide(_h_thrust[1][ibin1][ibin2],_h_thrust[0][ibin1][ibin2],htemp);
	  pair<double,double> asym = calcAsymmetry(htemp);
	  hthrustUL->addPoint(pRef.x(),asym.first,pRef.xErrs(),make_pair(asym.second,asym.second) );
	  // UC thrust
	  book(htemp,"TMP/R_thrust_UC_"+toString(ibin1)+"_"+toString(ibin2));
	  divide(_h_thrust[1][ibin1][ibin2],_h_thrust[2][ibin1][ibin2],htemp);
	  asym = calcAsymmetry(htemp);
	  hthrustUC->addPoint(pRef.x(),asym.first,pRef.xErrs(),make_pair(asym.second,asym.second) );
	  // UL hadron
	  book(htemp,"TMP/R_hadron_UL_"+toString(ibin1)+"_"+toString(ibin2));
	  divide(_h_hadron[1][ibin1][ibin2],_h_hadron[0][ibin1][ibin2],htemp);
	  asym = calcAsymmetry(htemp);
	  hhadronUL->addPoint(pRef.x(),asym.first,pRef.xErrs(),make_pair(asym.second,asym.second) );
	  // UC hadron
	  book(htemp,"TMP/R_hadron_UC_"+toString(ibin1)+"_"+toString(ibin2));
	  divide(_h_hadron[1][ibin1][ibin2],_h_hadron[2][ibin1][ibin2],htemp);
	  asym = calcAsymmetry(htemp);
	  hhadronUC->addPoint(pRef.x(),asym.first,pRef.xErrs(),make_pair(asym.second,asym.second) );
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_thrust[3][6][6],_h_hadron[3][6][6];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2014_I1254862);

}
