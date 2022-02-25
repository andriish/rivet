// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D meson production at 13 TeV
  class CMS_2021_I1876550 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2021_I1876550);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix) {
	book(_h_pT [ix],1,1,1+ix);
	book(_h_eta[ix],2,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==411 or Cuts::abspid==421 or
					     Cuts::abspid==413)) {
	// no mixing and |y|<2.1
	if(p.children().size()==1) continue;
	double eta = p.abseta();
	double pT=p.perp();
	if( eta>2.1 || pT<4. || pT>100.) continue;
	if(p.fromBottom()) continue;
	unsigned int imeson=0;
	if     (p.abspid()==421) imeson=1;
	else if(p.abspid()==411) imeson=2;
	_h_pT [imeson]->fill(pT);
	_h_eta[imeson]->fill(eta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection()/microbarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix) {
	scale(_h_pT[ix],factor);
	scale(_h_eta[ix],0.5*factor);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[3],_h_eta[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2021_I1876550);

}
