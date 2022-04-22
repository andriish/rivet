// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-pi0
  class ARGUS_1997_I420421 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1997_I420421);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if(inRange(sqrtS()/GeV,0.8,2.1)) {
	for(unsigned int ix=0;ix<5;++ix)
	  book(_nMeson[ix],"TMP/nMeson_"+toString(ix+1));
      }
      else
	throw Error("Invalid CMS energy for ARGUS_1997_I420421");
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
      // check the 3 meson final state
      if(ntotal!=3) vetoEvent;
      if( nCount[PID::PI0]==1 && nCount[PID::PIPLUS]==1 && nCount[PID::PIMINUS]==1 )
	_nMeson[0]->fill();
      else
	vetoEvent;
      // now the intermediate states
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // f0 and f2 mesons and rho
      bool nonRes=true;
      for(const Particle & p : ufs.particles(Cuts::pid==225 ||
					     Cuts::pid==9010221 ||
					     Cuts::abspid==213)) {
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
	int sign = p.pid()/p.abspid();
	// id of the pion not from the resonance decay
	int idother = p.abspid()==213 ? -sign*211 : 111;
	bool matched=true;
	for(auto const & val : nRes) {
	  if (val.first==idother) {
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
	if(matched) {
	  nonRes=false;
	  if(p.abspid()==213)
	    _nMeson[2]->fill();
	  else if (p.pid()==225)
	    _nMeson[3]->fill();
	  else
	    _nMeson[4]->fill();
	  break;
	}
      }
      if(nonRes) _nMeson[1]->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      // loop over tables in paper
      for(unsigned int ix=0;ix<5;++ix) {
	unsigned int iy=1;
	if(ix==2) iy=5;
	else if(ix==3) iy=3;
	double sigma = _nMeson[ix]->val()*fact;
	double error = _nMeson[ix]->err()*fact;
	Scatter2D temphisto(refData(ix+1, 1, iy));
	Scatter2DPtr mult;
	book(mult, ix+1, 1, iy);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/MeV, x-ex2.first, x+ex2.second)) {
	    mult->addPoint(x, sigma, ex, make_pair(error,error));
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
    CounterPtr _nMeson[5];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1997_I420421);

}
