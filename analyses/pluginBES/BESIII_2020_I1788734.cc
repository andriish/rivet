// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > phi eta'
  class BESIII_2020_I1788734 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1788734);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nEtaPhi, "TMP/phietaprime");
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
      for (const Particle& p : ufs.particles(Cuts::pid==PID::PHI)) {
        if(p.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        for (const Particle& p2 : ufs.particles(Cuts::pid==PID::ETAPRIME)) {
          if(p2.parents()[0].isSame(p)) continue;
          map<long,int> nResB = nRes;
          int ncountB = ncount;
          findChildren(p2,nResB,ncountB);
          if(ncountB!=0) continue;
          bool matched2 = true;
          for(auto const & val : nResB) {
            if(val.second!=0) {
              matched2 = false;
              break;
            }
          }
          if(matched2) {
            _nEtaPhi->fill();
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double sigma = _nEtaPhi->val()*crossSection()/ sumOfWeights() /picobarn;
      double error = _nEtaPhi->err()*crossSection()/ sumOfWeights() /picobarn;
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

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _nEtaPhi;
    ///@}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1788734);

}
