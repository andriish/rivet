// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief J/psi at 7 TeV
  class CMS_2011_I878118 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_I878118);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      vector<double> raps={0.,1.2,1.6,2.4};
      for(unsigned int iy=1;iy<raps.size();++iy) {
	Histo1DPtr tmp;
	for(unsigned int iprompt=0;iprompt<2;++iprompt) {
	  _h_psi[iprompt].add(raps[iy-1],raps[iy],book(tmp,1+iy+iprompt*9,1,1));
	}
	_h_psi[2].add(raps[iy-1],raps[iy],book(tmp,"TMP/total"+toString(iy),refData(4+iy,1,1)));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::pid==443)) {
	// cuts on rapidity
	double y = p.absrap();
	if(y>2.4) continue;
	double pT = p.perp();
	// prompt
	unsigned int iprompt = p.fromBottom();
	_h_psi  [iprompt].fill(y,pT);
	_h_psi  [   2   ].fill(y,pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // branching ratio
      double br = 0.05961;
      // 0.5 due folded rapidity
      double factor = 0.5*br*crossSection() / nanobarn/ sumOfWeights();
      for(unsigned int iprompt=0;iprompt<3;++iprompt)
	_h_psi[iprompt].scale(factor,this);
      // non-prompt fraction
      for(unsigned int iy=0;iy<_h_psi[1].histos().size();++iy) {
	Scatter2DPtr tmp;
	book(tmp,5+iy,1,1);
	efficiency(_h_psi[1].histos()[iy],_h_psi[2].histos()[iy],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_psi[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2011_I878118);

}
