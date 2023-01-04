// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi- and Sigma*+/- spectra
  class HRS_1987_I246162 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1987_I246162);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projects
      declare(UnstableParticles(Cuts::abspid==3114 || Cuts::abspid==3224 ||
				Cuts::abspid==3312 || Cuts::abspid==3122), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      // histos
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) book(_h_x[ix][iy],1+ix,1,1+iy);
	book(_h_x[ix][2] ,4,1,1+ix);
	book(_h_ratio[ix],3,1,1+ix);
      }
      book(_c,"TMP/nLam");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const size_t numParticles = cfs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      for(const Particle & p : apply<UnstableParticles>(event,"UFS").particles()) {
	if(p.abspid()==3122) {
	  _c->fill();
	  continue;
	}
      	double xE = 2.*p.E()/sqrtS();
      	Vector3 mom3 = p.p3();
        const double energy = p.E();
      	double modp = mom3.mod();
      	double beta = modp/energy;
	unsigned int iloc = p.abspid()==3312 ? 1 : 0;
	_h_ratio[iloc]->fill(sqrtS());
	_h_x[iloc][0]->fill(xE);
	_h_x[iloc][1]->fill(xE);
	_h_x[iloc][2]->fill(xE,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_x[ix][0],0.8/sumOfWeights());
	scale(_h_x[ix][1],1.0/sumOfWeights());
	scale(_h_x[ix][2], sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
	scale(_h_ratio[ix], 1./ *_c);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_x[2][3],_h_ratio[2];
    CounterPtr _c;

    /// @}


  };


  RIVET_DECLARE_PLUGIN(HRS_1987_I246162);

}
