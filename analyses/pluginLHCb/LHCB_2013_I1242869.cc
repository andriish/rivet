// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_c at 7 TeV
  class LHCB_2013_I1242869 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2013_I1242869);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      for(unsigned int ichi=0;ichi<3;++ichi) {
	book(_h_chi[ichi],"TMP/h_CHI_"+toString(ichi),refData(1,1,1));
	book(_c_chi[ichi],"TMP/c_CHI_"+toString(ichi));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==10441 ||
					     Cuts::pid==20443 ||
					     Cuts::pid==445)) {
	// prompt
	if(p.fromBottom()) continue;
	// J/psi /gamma mode
	if(p.children().size()!=2) continue;
	Particle Jpsi;
	if(p.children()[0].pid()==22 && p.children()[1].pid()==443) {
	  Jpsi=p.children()[1];
	}
	else if(p.children()[1].pid()==22 && p.children()[0].pid()==443) {
	  Jpsi=p.children()[0];
	}
	else
	  continue;
	double absrap=Jpsi.absrap();
	if(absrap<2. || absrap>4.5) continue;
	unsigned int ichi = 0;
	if(p.pid()==20443) ichi=1;
	else if(p.pid()==445) ichi=2;
	double xp=Jpsi.perp();
	_h_chi[ichi]->fill(xp);
	if(xp>4. && xp<20.) _c_chi[ichi]->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // chi_c to gamma J/psi branching ratios from PDG 2021
      vector<double> br = {0.014,0.343,0.190};
      // divide out the branching ratio for mode used to get total rate
      for(unsigned int ichi=0;ichi<3;++ichi) {
	scale(_h_chi[ichi],1./br[ichi]);
	scale(_c_chi[ichi],1./br[ichi]);
      }
      // ratio chi_c2/chi_c2
      Scatter2DPtr tmp;
      book(tmp,1,1,1);
      divide(_h_chi[2],_h_chi[1],tmp);
      Scatter1D s1("/LHCB_2013_I1242869/TMP/r_1");
      Scatter1DPtr s1d = registerAO(s1);
      Scatter1D s2("/LHCB_2013_I1242869/TMP/r_2");
      Scatter1DPtr s2d = registerAO(s2);
      divide(_c_chi[0],_c_chi[2],s1d);
      divide(_c_chi[2],_c_chi[1],s2d);
      Scatter2D temphisto(refData(2, 1, 1));
      Scatter2DPtr ratio;
      book(ratio, 2, 1, 1);
      for (size_t b = 0; b < temphisto.numPoints(); b++) {
	const double x  = temphisto.point(b).x();
	pair<double,double> ex = temphisto.point(b).xErrs();
	pair<double,double> ex2 = ex;
	if(ex2.first ==0.) ex2. first=0.0001;
	if(ex2.second==0.) ex2.second=0.0001;
	if (inRange(1, x-ex2.first, x+ex2.second)) {
	  ratio->addPoint(x,s1d->points()[0].x(),ex,s1d->points()[0].xErrs());
	}
	else if (inRange(2, x-ex2.first, x+ex2.second)) {
	  ratio->addPoint(x,s2d->points()[0].x(),ex,s2d->points()[0].xErrs());
	}
	else {
	  ratio->addPoint(x, 0., ex, make_pair(0.,.0));
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_chi[3];
    CounterPtr _c_chi[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2013_I1242869);

}
