// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B_s fracion at 7 TeV
  class LHCB_2013_I1215607 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2013_I1215607);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
      	book(_h_pT[ix],"TMP/h_pT_"+toString(ix),refData(1,1,1));
	book(_h_y [ix],"TMP/h_y_" +toString(ix),refData(2,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::abspid==531 || Cuts::abspid==511)) {
      	// skip copies due mixing
      	if(p.children().size()==1 && p.children()[0].abspid()==p.abspid()) continue;
      	double eta=p.abseta();
      	if(eta<2. || eta>5.) continue;
      	double pT = p.perp();
      	if(pT<1.5 || pT>40.) continue;
      	if(p.abspid()==531) {
      	  _h_pT[0]->fill(pT);
      	  _h_y[0]->fill(eta);
      	}
      	else {
      	  _h_pT[1]->fill(pT);
      	  _h_y[1]->fill(eta);
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Scatter2DPtr tmp;
      book(tmp,1,1,1);
      divide(_h_pT[0],_h_pT[1],tmp);
      book(tmp,2,1,1);
      divide(_h_y [0],_h_y[1],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[2],_h_y[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2013_I1215607);

}
