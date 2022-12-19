// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Measurement of R 3.65 and 3.872 GeV
  class BESII_2006_I735496 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2006_I735496);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==100443 ||
				Cuts::abspid==30443),"UFS");
      declare(FinalState(), "FS");
      // Book histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_c_hadrons[ix], "/TMP/sigma_hadrons_"+toString(ix+1));
      book(_c_muons,   "/TMP/sigma_muons");
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
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_c_muons->fill();
      // everything else
      else {
       	_c_hadrons[1]->fill();
	Particles psi = apply<UnstableParticles>(event,"UFS").particles();
	if(psi.empty())
	  _c_hadrons[0]->fill();
	else {
	  bool psi3770 = false;
	  for(const Particle & p : psi) psi3770 |= p.pid()==30443;
	  if(psi3770) _c_hadrons[0]->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	Scatter1D R = *_c_hadrons[ix]/ *_c_muons;
	double              rval = R.point(0).x();
	pair<double,double> rerr = R.point(0).xErrs();
	Scatter2D temphisto(refData(1, 1, 1+ix));
	Scatter2DPtr     mult;
	book(mult, 1, 1, 1+ix);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  pair<double,double> ex2 = ex;
	  if(ex2.first ==0.) ex2. first=0.0001;
	  if(ex2.second==0.) ex2.second=0.0001;
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second))
	    mult   ->addPoint(x, rval, ex, rerr);
	  else
	    mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_hadrons[2], _c_muons;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2006_I735496);

}
