// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi from chi_c decays
  class CDF_1997_I440446 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_1997_I440446);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_total[0],"TMP/h_chi",refData(1,1,1));
      book(_h_total[1],"TMP/h_psi",refData(1,1,1));
      for(unsigned int ix=0;ix<3;++ix)
	book(_h_psi[ix],2+ix,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // first J/psi for denominator
      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
	if(p.fromBottom()) continue;
        double abseta = p.abseta();
        double xp = p.perp();
	if (xp<4. || abseta>0.6) continue;
	_h_total[1]->fill(sqrtS());
	// from those for higher charmonium
	Particle parent = p.parents()[0];
	if ( parent.pid()==100443 || parent.pid()==20443 || parent.pid()==445) continue;
	_h_psi[0]->fill(xp);
      }
      // chi_1 and chi_2 for numerator
      for (const Particle& p : ufs.particles(Cuts::pid==20443 || Cuts::pid==100443 || Cuts::pid==445)) {
	if(p.fromBottom()) continue;
	Particle Jpsi;
	bool found(false);
	for (const Particle & child : p.children()) {
	  if (child.pid()==443) {
	    found = true;
	    Jpsi=child;
	  }
	}
	if(!found) continue;
	double abseta=Jpsi.abseta();
	double xp = Jpsi.perp();
	if (xp<4. || abseta>0.6) continue;
	if (p.pid()==100443) {
	  _h_psi[2]->fill(xp);
	}
	else {
	  _h_psi[1]->fill(xp);
	  _h_total[0]->fill(sqrtS());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalisation factor
      double br = 0.05961;
      double factor = br*crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix)
	scale(_h_psi[ix],factor);   
      Scatter2DPtr tmp;
      book(tmp,1,1,1);
      efficiency(_h_total[0],_h_total[1],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_psi[3],_h_total[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CDF_1997_I440446);

}
