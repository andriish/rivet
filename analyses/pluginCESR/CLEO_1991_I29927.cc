// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CLEO_1991_I29927 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEO_1991_I29927);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      _c_B     = bookCounter("/TMP/sigma_B");
      _c_Bstar = bookCounter("/TMP/sigma_Bstar");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      unsigned int nBstar(0);
      // Get Bottom hadrons
      const Particles bhads = filter_select(ufs.particles(), isBottomHadron);
      // find the Bstars
      for (const Particle& p : bhads) {
        if(abs(p.pdgId())==513 || abs(p.pdgId())==523) {
          if(!p.hasDescendantWith(Cuts::pid == p.pdgId())) ++nBstar;
        }
      }
      if(!bhads.empty())
        _c_B->fill(event.weight());
      if(nBstar!=0)
        _c_Bstar->fill(nBstar*event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=1;ix<3;++ix) {
        double sig(0.),err(0.);
        if(ix==1) {
          sig = _c_B->val()*fact;
          err = _c_B->err()*fact;
        }
        else {
          sig = _c_Bstar->val()*fact;
          err = _c_Bstar->err()*fact;
        }
        Scatter2D    temphisto(refData(ix, 1, 1));
        Scatter2DPtr mult = bookScatter2D(ix, 1, 1);
        for (size_t b = 0; b < temphisto.numPoints(); b++) {
          const double x  = temphisto.point(b).x();
          pair<double,double> ex = temphisto.point(b).xErrs();
          pair<double,double> ex2 = ex;
          if(ex2.first ==0.) ex2. first=0.0001;
          if(ex2.second==0.) ex2.second=0.0001;
          if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
            mult->addPoint(x, sig, ex, make_pair(err,err));
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
    CounterPtr _c_B, _c_Bstar;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CLEO_1991_I29927);


}
