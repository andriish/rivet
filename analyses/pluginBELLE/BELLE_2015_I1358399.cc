// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2015_I1358399 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2015_I1358399);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histograms
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // extract botton hadrons
      Particles bHadrons=apply<FinalState>(event, "UFS").particles(Cuts::abspid==531 or Cuts::abspid==533);
      for(unsigned int ix=0;ix<bHadrons.size();++ix) {
       	if(bHadrons[ix].parents()[0].abspid()==531 ||
	   bHadrons[ix].parents()[0].abspid()==533 ) continue;
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(bHadrons[ix],nRes,ncount);
	bool matched=false;
	for(unsigned int iy=ix+1;iy<bHadrons.size();++iy) {
	  if(bHadrons[iy].parents()[0].abspid()==531 ||
	     bHadrons[iy].parents()[0].abspid()==533 ) continue;
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(bHadrons[iy],nRes2,ncount2);
	  if(ncount2!=0) continue;
	  matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
      	  if(matched) {
	    _h->fill(10.86);
	    break;
      	  }
      	}
      	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h,crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2015_I1358399);

}
