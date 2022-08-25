// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_c -> Lambda0 or Sigma0 + (pi,K)+ decay asymmetries
  class BELLE_2022_I2138841 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2138841);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      for(unsigned int imode=0;imode<4;++imode) {
	if(imode<2) {
	  book(_h[imode][0],3,1,1+imode);
	  for(unsigned int iy=0;iy<2;++iy)
	    book(_h[imode][1+iy],4,1,1+iy+2*imode);
	}
	else {
	  for(unsigned int iy=0;iy<3;++iy) {
	    book(_h[imode][iy],"h_"+toString(imode+1)+"_"+toString(iy+1),10,-1,1);
	  }
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Lambda_c baryons
      for( const Particle& Lambdac : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==4122)) {
	int sign = Lambdac.pid()/4122;
	if(Lambdac.children().size()!=2) continue;
	Particle baryon1;
	int imeson=-1;
	if((Lambdac.children()[0].pid()==sign*3122 ||
	    Lambdac.children()[0].pid()==sign*3212) && 
	   Lambdac.children()[1].pid()==sign*321) {
	  baryon1 = Lambdac.children()[0];
	  imeson=0;
	}
	else if((Lambdac.children()[1].pid()==sign*3122 ||
		 Lambdac.children()[0].pid()==sign*3212) && 
		Lambdac.children()[0].pid()==sign*321) {
	  baryon1 = Lambdac.children()[1];
	  imeson=0;
	}
	else if((Lambdac.children()[0].pid()==sign*3122 ||
		 Lambdac.children()[0].pid()==sign*3212) && 
		Lambdac.children()[1].pid()==sign*211) {
	  baryon1 = Lambdac.children()[0];
	  imeson=1;
	}
	else if((Lambdac.children()[1].pid()==sign*3122 ||
		 Lambdac.children()[0].pid()==sign*3212) && 
		Lambdac.children()[0].pid()==sign*211) {
	  baryon1 = Lambdac.children()[1];
	  imeson=1;
	}
	else
	  continue;
	// Lambda0 case
	if(baryon1.abspid()==3122) {
	  Particle baryon2;
	  if(baryon1.children()[0].pid()== sign*2212 && 
	     baryon1.children()[1].pid()==-sign*211) {
	    baryon2 = baryon1.children()[0];
	  }
	  else if(baryon1.children()[1].pid()== sign*2212 && 
		  baryon1.children()[0].pid()==-sign*211) {
	    baryon2 = baryon1.children()[1];
	  }
	  else
	    continue;
	  // first boost to the Lambdac rest frame
	  LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lambdac.momentum().betaVec());
	  FourMomentum pbaryon1 = boost1.transform(baryon1.momentum());
	  FourMomentum pbaryon2 = boost1.transform(baryon2.momentum());
	  // to lambda rest frame
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	  Vector3 axis = pbaryon1.p3().unit();
	  FourMomentum pp = boost2.transform(pbaryon2);
	  // calculate angle
	  double cTheta = pp.p3().unit().dot(axis);
	  _h[imeson][0]->fill(cTheta);
	  if(baryon1.pid()>0)
	    _h[imeson][1]->fill(cTheta);
	  else
	    _h[imeson][2]->fill(cTheta);
	}
	// sigma0 case
	else {
	  Particle baryon2;
	  if(baryon1.children()[0].pid()== sign*3122 && 
	     baryon1.children()[1].pid()== 22) {
	    baryon2 = baryon1.children()[0];
	  }
	  else if(baryon1.children()[1].pid()== sign*3122 && 
		  baryon1.children()[0].pid()== 22) {
	    baryon2 = baryon1.children()[1];
	  }
	  else
	    continue;
	  Particle baryon3;
	  if(baryon2.children()[0].pid()== sign*2212 && 
	     baryon2.children()[1].pid()==-sign*211) {
	    baryon3 = baryon2.children()[0];
	  }
	  else if(baryon2.children()[1].pid()== sign*2212 && 
		  baryon2.children()[0].pid()==-sign*211) {
	    baryon3 = baryon2.children()[1];
	  }
	  else
	    continue;
	  // first boost to the Lambdac rest frame
	  LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lambdac.momentum().betaVec());
	  FourMomentum pbaryon1 = boost1.transform(baryon1.momentum());
	  FourMomentum pbaryon2 = boost1.transform(baryon2.momentum());
	  FourMomentum pbaryon3 = boost1.transform(baryon3.momentum());
	  // to  sigma rest frame
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	  Vector3 axis = pbaryon1.p3().unit();
	  FourMomentum pp  = boost2.transform(pbaryon2);
	  FourMomentum pp3 = boost2.transform(pbaryon3);
	  // calculate angle
	  double cTheta2 = pp.p3().unit().dot(axis);
	  // to lambda rest frame
	  LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pp.betaVec());
	  Vector3 axis2 = pp.p3().unit();
	  FourMomentum pp4 = boost3.transform(pp3);
	  // calculate angle
	  double cTheta3 = pp4.p3().unit().dot(axis);
	  double cTheta = cTheta2*cTheta3;
	  _h[imeson+2][0]->fill(cTheta);
	  if(baryon1.pid()>0)
	    _h[imeson+2][1]->fill(cTheta);
	  else
	    _h[imeson+2][2]->fill(cTheta);
	}
      }
    }

    pair<double,double> calcAlpha(Histo1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (auto bin : hist->bins() ) {
	double Oi = bin.area();
	if(Oi==0.) continue;
	double ai = 0.5*(bin.xMax()-bin.xMin());
	double bi = 0.5*ai*(bin.xMax()+bin.xMin());
	double Ei = bin.areaErr();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      pair<double,double> aLambda(0.7542,0.0022); 
      for(int imeson=0;imeson<4;++imeson) {
	for(int iy=0;iy<3;++iy) {
	  normalize(_h[imeson][iy]);
	  pair<double,double> alpha = calcAlpha(_h[imeson][iy]);
	  Scatter2DPtr _h_alpha1,_h_alpha2;
	  if(iy==0) {
	    book(_h_alpha1,1,1+imeson,1);
	    book(_h_alpha2,1,1+imeson,2);
	  }
	  else {
	    book(_h_alpha1,2,1+imeson,iy);
	    book(_h_alpha2,2,1+imeson,2+iy);
	  }
	  _h_alpha1->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
	  // divide out aLambda
	  alpha.second = alpha.first/aLambda.first*
	    sqrt(sqr(alpha.second/alpha.first) + sqr(aLambda.second/aLambda.first));
	  alpha.first /= aLambda.first;
	  _h_alpha2->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2138841);

}
