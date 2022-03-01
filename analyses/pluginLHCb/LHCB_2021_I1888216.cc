// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_b asymmetry
  class LHCB_2021_I1888216 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2021_I1888216);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      //histograms
      for(unsigned int ix=0;ix<2;++ix) {
	book(_h_y [ix],"TMP/h_y_ "+toString(ix+1),refData(1,1,1));
	book(_h_pT[ix],"TMP/h_pT_"+toString(ix+1),refData(2,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==5122) ) {
	double pT = p.perp();
	double y  = p.absrap();
	if(y<2.15 ||  y>4.1) continue;
	if(pT<2   || pT>27 ) continue;
	bool anti = p.pid()<0;
	_h_pT[anti]->fill(pT);
	_h_y [anti]->fill(y);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      unsigned int is=0;
      if (isCompatibleWithSqrtS(7000)) {
      	is=1;
      }
      else if (isCompatibleWithSqrtS(8000)) {
      	is=2;
      }
      else  {
      	throw Error("Invalid CMS energy for LHCB_2021_I1888216");
      }
      // calculate asymmetry in %
      Scatter2DPtr tmp;
      book(tmp,1,1,is);
      asymm(_h_y [0], _h_y[1],tmp);
      tmp->scaleY(100.);
      book(tmp,2,1,is);
      asymm(_h_pT[0],_h_pT[1],tmp);
      tmp->scaleY(100.);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT[2],_h_y[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2021_I1888216);

}
