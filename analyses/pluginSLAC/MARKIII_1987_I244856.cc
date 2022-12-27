// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e-> D_s D_s*
  class MARKIII_1987_I244856 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKIII_1987_I244856);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==431 or Cuts::abspid==433), "UFS");
      // histos
      book(_h,1,1,1);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
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
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) vetoEvent;
      Particles ds = apply<UnstableParticles>(event, "UFS").particles();
      for(const Particle & p1 : ds) {
	if(p1.abspid()!=431) continue;
      	// find the children
	bool matched = false;
      	map<long,int> nRes = nCount;
      	int ncount = ntotal;
      	findChildren(p1,nRes,ncount);
	for(const Particle & p2 : ds) {
	  if(p2.abspid()!=433 || p1.pid()*p2.pid()>0) continue;
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount2!=0) continue;
	  matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    break;
	  }
	}
	if(matched) {
	  _h->fill(sqrtS());
	  break;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double br = 0.045;
      double fact = br*crossSection()/ sumOfWeights()/picobarn;
      scale(_h,fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(MARKIII_1987_I244856);

}
