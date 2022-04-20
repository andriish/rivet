// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> K*K*
  class ARGUS_2000_I511512 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_2000_I511512);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if(inRange(sqrtS()/GeV,1.5,3.0)) {
	for(unsigned int ix=0;ix<9;++ix)
	  book(_nMeson[ix],"TMP/nMeson_"+toString(ix+1));
      }
      else
	throw Error("Invalid CMS energy for ARGUS_2000_I511512");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]-=1;
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal==4) {
	if( nCount[PID::KPLUS]==1 && nCount[PID::KMINUS]==1 &&
	    nCount[PID::PIPLUS]==1 && nCount[PID::PIMINUS]==1 )
	  _nMeson[6]->fill();
	else if( nCount[PID::K0S]==2 && 
		 nCount[PID::PIPLUS]==1 && nCount[PID::PIMINUS]==1 )
	  _nMeson[8]->fill();
	else if( nCount[PID::K0S]==1 && nCount[PID::PI0]==1 && 
		 ((nCount[PID::KPLUS ]==1 && nCount[PID::PIMINUS]==1) || 
		  (nCount[PID::KMINUS]==1 && nCount[PID::PIPLUS ]==1)))
	  _nMeson[7]->fill();

      }
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // find any K* mesons
      Particles Kstar=ufs.particles(Cuts::abspid==313 or Cuts::abspid==323);
      for (unsigned int ix=0;ix<Kstar.size();++ix) {
       	if(Kstar[ix].children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(Kstar[ix],nRes,ncount);
	int sign = Kstar[ix].pid()/Kstar[ix].abspid();
	// three body intermediate states
	if(ncount==2) {
	  // K*0 K- pi+ +ccd
	  if(Kstar[ix].abspid()==313) {
	    bool matched=true;
	    for(auto const & val : nRes) {
	      if (val.first==sign*211 || val.first==-sign*321) {
		if(val.second!=1) {
		  matched = false;
		  break;
		}
	      }
	      else {
		if(val.second!=0) {
		  matched = false;
		  break;
		}
	      }
	    }
	    if(matched) _nMeson[3]->fill();
	  }
	  else {
	    bool matched=false;
	    // K*+ K0S pi- + cc
	    for(auto const & val : nRes) {
	      if (val.first==-sign*211 || val.first==PID::K0S) {
		if(val.second!=1) {
		  matched = false;
		  break;
		}
	      }
	      else {
		if(val.second!=0) {
		  matched = false;
		  break;
		}
	      }
	    }
	    if(matched) _nMeson[4]->fill();
	    else {
	      // K*+ K- pi0 +cc
	      matched=false;
	      for(auto const & val : nRes) {
		if (val.first==-sign*321 || val.first==PID::PI0) {
		  if(val.second!=1) {
		    matched = false;
		    break;
		  }
		}
		else {
		  if(val.second!=0) {
		    matched = false;
		    break;
		  }
		}
	      }
	      if(matched) _nMeson[5]->fill();
	    }
	  }
	}
	       
	// K*K*
	for (unsigned int iy=ix+1;iy<Kstar.size();++iy) {
	  if(Kstar[iy].children().empty()) continue;
	  if(Kstar[ix].pid()!=-Kstar[iy].pid()) continue;
	  map<long,int> nRes2=nRes;
	  int ncount2 = ncount;
	  findChildren(Kstar[iy],nRes2,ncount2);
	  if(ncount2 !=0 ) continue;
	  bool matched2 = true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched2 = false;
	      break;
	    }
	  }
	  if(matched2) {
	    if(Kstar[ix].abspid()==313)
	      _nMeson[0]->fill();
	    else
	      _nMeson[1]->fill();
	    break;
	  }
	}
      }
      // finally the rho phi intermediate states 
      for (const Particle & p1 : ufs.particles(Cuts::pid==PID::RHO0)) {
       	if(p1.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p1,nRes,ncount);
	for (const Particle & p2 : ufs.particles(Cuts::pid==PID::PHI)) {
	  if(p2.children().empty()) continue;
	  map<long,int> nRes2=nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount2 !=0 ) continue;
	  bool matched = true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    _nMeson[2]->fill();
	    break;
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      // loop over tables in paper
      for(unsigned int ih=3;ih<11;++ih) {
	unsigned int imax=5;
	if(ih==6||ih==9||ih==10) imax=3;
	else if(ih==7||ih==8)  imax=4;
	for(unsigned int iy=1;iy<imax;++iy) {
	  unsigned int iloc=1000;
	  // K*0 K*bar0
	  if( (iy==1 && (ih==3||ih==4||ih==5||ih==9||ih==10)) ||
	      (iy==2&&ih==6))
	    iloc=0;
	  // K*+ K*-
	  else if ( (iy==1 && (ih==6 || ih==8)) ||
		    (iy==2 && (ih==9 || ih==10)))
	    iloc=1;
	  // rho phi
	  else if(iy==3 && (ih>=3&&ih<=5))
	    iloc=2;
	  // K*0 K- pi+ + cc
	  else if(iy==2 && (ih>=3&&ih<=5))
	    iloc=3;
	  // K+K-pi+pi-
	  else if(iy==4 && (ih>=3&&ih<=5))
	    iloc=6;
	  // K*+ K0S pi-+cc
	  else if( (ih==7&&iy==1) || (ih==8&&iy==2))
	    iloc=4;
	  // K*+ K- pi0 +cc
	  else if(ih==7&&iy==2)
	    iloc=5;
	  // KS0 K+ pi- pi0 +cc
	  else if(ih==7&&iy==3)
	    iloc=7;
	  // KS0 KS0 pi+ pi-
	  else if(ih==8&&iy==3)
	    iloc=8;
	  assert(iloc<=8);
	  double sigma = _nMeson[iloc]->val()*fact;
	  double error = _nMeson[iloc]->err()*fact;
	  Scatter2D temphisto(refData(ih, 1, iy));
	  Scatter2DPtr mult;
	  book(mult, ih, 1, iy);
	  for (size_t b = 0; b < temphisto.numPoints(); b++) {
	    const double x  = temphisto.point(b).x();
	    pair<double,double> ex = temphisto.point(b).xErrs();
	    pair<double,double> ex2 = ex;
	    if(ex2.first ==0.) ex2. first=0.0001;
	    if(ex2.second==0.) ex2.second=0.0001;
	    if (inRange(sqrtS(), x-ex2.first, x+ex2.second)) {
	      mult->addPoint(x, sigma, ex, make_pair(error,error));
	    }
	    else {
	      mult->addPoint(x, 0., ex, make_pair(0.,.0));
	    }
	  }
	}
      }
      // finally the ratio K*+/K*0
      if(_nMeson[0]->numEntries()>0) {
	Scatter1D R = *_nMeson[1]/ *_nMeson[0];
	double              rval = R.point(0).x();
	pair<double,double> rerr = R.point(0).xErrs();
	Scatter2D temphisto(refData(11, 1, 1));
	Scatter2DPtr mult;
	book(mult, 11, 1, 1);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS(), x-ex2.first, x+ex2.second)) {
	    mult->addPoint(x, rval, ex, rerr);
	  }
	  else {
	    mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nMeson[9];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_2000_I511512);

}
