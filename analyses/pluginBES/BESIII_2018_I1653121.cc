// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > KK J/Psi
  class BESIII_2018_I1653121 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1653121);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nKp, "TMP/Kp");
      book(_nK0, "TMP/K0");

    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
          --nRes[child.pid()];
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
        // find the J/psi
        if(p.pid()==443) {
          map<long,int> nRes = nCount;
          int ncount = ntotal;
          findChildren(p,nRes,ncount);
          // omega pi+pi-
          if(ncount!=2) continue;
          bool matched = true;
          for(auto const & val : nRes) {
            if(abs(val.first)==321 || abs(val.first)==310) {
              continue;
            }
            else if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            if(nRes[321]==1 && nRes[-321]==1) {
              _nKp->fill();
              break;
            }
            else if(nRes[310]==2) {
              _nK0->fill();
              break;
            }
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=1;iy<3;++iy) {
        double sigma,error;
        if(iy==1) {
          sigma = _nKp->val();
          error = _nKp->err();
        }
        else {
          sigma = _nK0->val();
          error = _nK0->err();
        }
        sigma *= crossSection()/ sumOfWeights() /picobarn;
        error *= crossSection()/ sumOfWeights() /picobarn;
        Scatter2D temphisto(refData(1, 1, iy));
        Scatter2DPtr  mult;
        book(mult, 1, 1, iy);
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
    CounterPtr _nKp, _nK0;
    //@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1653121);

}
