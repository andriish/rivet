// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief 
  class BESIII_2020_I1808875 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1808875);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      declare(FinalState(), "FS");
      book(_c_muons, "/TMP/sigma_muons");
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
      if(nCount[-13]==1 && nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_c_muons->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sig_m = _c_muons  ->val()*fact;
      double err_m = _c_muons  ->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr muons;
      book(muons, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second))
	  muons  ->addPoint(x, sig_m, ex, make_pair(err_m,err_m));
	else
	  muons  ->addPoint(x, 0., ex, make_pair(0.,.0));
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_muons;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1808875);

}
