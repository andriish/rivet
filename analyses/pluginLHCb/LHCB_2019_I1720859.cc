// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B hadron fractions at 13 TeV
  class LHCB_2019_I1720859 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2019_I1720859);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<3;++ix) {
	book(_h_pT[ix],"TMP/h_pT_"+toString(ix),refData(1,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::abspid==5122 || Cuts::abspid==511 || Cuts::abspid==521 || Cuts::abspid==531)) {
	// skip copies due mixing
	if(p.children().size()==1 && p.children()[0].abspid()==p.abspid()) continue;
	double eta=p.abseta();
	if(eta<2. || eta>5.) continue;
	double pT = p.perp();
	if(pT<4. || pT>25.) continue;
	if(p.abspid()==5122)
	  _h_pT[0]->fill(pT);
	else if(p.abspid()==531)
	  _h_pT[1]->fill(pT);
	else
	  _h_pT[2]->fill(pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Scatter2DPtr tmp;
      book(tmp,1,1,2);
      divide(_h_pT[0],_h_pT[2],tmp);
      book(tmp,1,1,1);
      divide(_h_pT[1],_h_pT[2],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2019_I1720859);

}
