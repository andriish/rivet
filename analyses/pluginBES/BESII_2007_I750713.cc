// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > light hadrons
  class BESII_2007_I750713 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2007_I750713);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      for(unsigned int ix=1;ix<19;++ix) {
	stringstream ss;
	ss << "TMP/n" << ix;
	book(_nMeson[ix], ss.str());
      }

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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      
      for (const Particle& p : ufs.particles()) {
      	if(p.children().empty()) continue;
      	if(p.pid()!=221 && p.pid()!=333) continue;
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p,nRes,ncount);
	// eta
      	if(p.pid()==221) {
       	  if(ncount==4) {
      	    bool matched = true;
      	    for(auto const & val : nRes) {
      	      if(abs(val.first)==211) {
      		if(val.second!=2) {
      		  matched = false;
      		  break;
      		}
      	      }
      	      else if(val.second!=0) {
      		matched = false;
      		break;
      	      }
      	    }
      	    if(matched) _nMeson[12]->fill();
      	  }
	}
       	else if(p.pid()==333) {
       	  if(ncount!=1) continue;
       	  bool matched = true;
       	  for(auto const & val : nRes) {
       	    if(abs(val.first)==321) {
       	      if(val.second!=1) {
       		matched = false;
       		break;
       	      }
       	    }
       	    else if(val.second!=0) {
       	      matched = false;
       	      break;
       	    }
       	  }
       	  if(matched)
       	    _nMeson[7]->fill();
      	}
      }
      if(ntotal==3 &&  nCount[111]==1 &&
	 nCount[-2212] == 1 && nCount[ 2212]==1)
	_nMeson[16]->fill();
      else if(ntotal==4) {
	if(nCount[-211] == 2 && nCount[ 211]==2)
	  _nMeson[3]->fill();
	else if(nCount[-211] == 1 && nCount[ 211]==1 &&
		nCount[-321] == 1 && nCount[ 321]==1)
	  _nMeson[4]->fill();
	else if(nCount[-321] == 2 && nCount[ 321]==2)
	  _nMeson[6]->fill();
	else if(nCount[-211 ] == 1 && nCount[ 211 ]==1 &&
		nCount[-2212] == 1 && nCount[ 2212]==1)
	  _nMeson[8]->fill();
	else if(nCount[-321 ] == 1 && nCount[ 321 ]==1 &&
		nCount[-2212] == 1 && nCount[ 2212]==1)
	  _nMeson[9]->fill();
      }
      else if(ntotal==5 && nCount[111]==1) {
	if(nCount[-211] == 2 && nCount[ 211]==2)
	  _nMeson[13]->fill();
	else if(nCount[-211] == 1 && nCount[ 211]==1 &&
		nCount[-321] == 1 && nCount[ 321]==1)
	  _nMeson[14]->fill();
	else if(nCount[-321] == 2 && nCount[ 321]==2)
	  _nMeson[15]->fill();
	else if(nCount[-211 ] == 1 && nCount[ 211 ]==1 &&
		nCount[-2212] == 1 && nCount[ 2212]==1)
	  _nMeson[17]->fill();
      }
      else if(ntotal==6 && nCount[211]==3 && nCount[-211]==3)
	_nMeson[11]->fill();
      else if(ntotal==7 && nCount[111]==1 &&
	      nCount[211]==3 && nCount[-211]==3)
	_nMeson[18]->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=3;ix<19;++ix) {
        if(ix==5 || ix==10) continue;
	double sigma = _nMeson[ix]->val();
	double error = _nMeson[ix]->err();
    	sigma *= crossSection()/ sumOfWeights() /picobarn;
    	error *= crossSection()/ sumOfWeights() /picobarn; 
	Scatter2D temphisto(refData(1, 1, ix));
    	Scatter2DPtr  mult;
        book(mult, 1, 1, ix);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	    mult->addPoint(x, sigma, ex, make_pair(error,error));
	  }
	  else {
	    mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }
    //@}

    /// @name Histograms
    //@{
    CounterPtr _nMeson[19];
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESII_2007_I750713);


}
