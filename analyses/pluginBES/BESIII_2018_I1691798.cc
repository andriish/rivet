// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Cross-section for $K^0_SK^\pm\pi^\mp$ between 3.8 and 4.6 GeV
  class BESIII_2018_I1691798 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1691798);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_nKKpi, "/TMP/nKKpi" );
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

      if(ntotal==3 && nCount[310]==1 &&
         ((nCount[ 321]=1 &&  nCount[-211] ==1) ||
          (nCount[-321]=1 &&  nCount[ 211] ==1)))
        _nKKpi->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nKKpi->val();
      double error = _nKKpi->err();
      sigma *= crossSection()/ sumOfWeights() /picobarn;
      error *= crossSection()/ sumOfWeights() /picobarn;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr  mult;
      book(mult, 1, 1, 1);
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

    //@}


    /// @name Histograms
    //@{
    CounterPtr _nKKpi;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2018_I1691798);


}
