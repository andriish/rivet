// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Measurement of the hadronic cross section for energies between 3.08 and 3.116 GeV
  class BES_1995_I39870 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_1995_I39870);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons");
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
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) vetoEvent;
      // everything else
      else _c_hadrons->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sig_h = _c_hadrons->val()*fact;
      double err_h = _c_hadrons->err()*fact;
      Scatter2D temphisto(refData(1, 1, 1));
      Scatter2DPtr hadrons;
      book(hadrons, 1, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
        const double x  = temphisto.point(b).x();
        pair<double,double> ex = temphisto.point(b).xErrs();
        pair<double,double> ex2 = ex;
        if(ex2.first ==0.) ex2. first=0.0001;
        if(ex2.second==0.) ex2.second=0.0001;
        if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
          hadrons->addPoint(x, sig_h, ex, make_pair(err_h,err_h));
        } else {
          hadrons->addPoint(x, 0., ex, make_pair(0.,.0));
        }
      }

    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c_hadrons;
    //@}

  };


  RIVET_DECLARE_PLUGIN(BES_1995_I39870);

}
