// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_c+ at 5.02 TeV
  class CMS_2020_I1738943 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2020_I1738943);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_Lambda,1,1,1);
      book(_h_D     ,"TMP/h_D",refData(4,1,1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::abspid==4122 || Cuts::abspid==421)) {
	// skip copies due mixing
	if(p.children().size()==1 && p.children()[0].abspid()==p.abspid()) continue;
	if(p.absrap()>1.) continue;
	double pT = p.perp();
	if(p.abspid()==4122)
	  _h_Lambda->fill(pT);
	else
	  _h_D->fill(pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 0.5 from particle/antiparticle
      double fact = 0.5*crossSection()/microbarn/sumOfWeights();
      scale(_h_Lambda,fact);
      scale(_h_D     ,fact);
      Scatter2DPtr tmp;
      book(tmp,4,1,1);
      divide(_h_Lambda,_h_D,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Lambda,_h_D;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2020_I1738943);

}
