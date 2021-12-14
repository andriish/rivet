// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Quarkonium prouction at 13 TeV
  class CMS_2018_I1633431 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1633431);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // J/psi/psi(2s) histos
      vector<double> raps={0.,0.3,0.6,0.9,1.2};
      for(unsigned int ix=1;ix<3;++ix) {
	for(unsigned int iy=1;iy<5;++iy) {
	  Histo1DPtr tmp;
	  _h_onium[ix-1].add(raps[iy-1],raps[iy],book(tmp,ix,1,iy));
	}
	book(_h_total[ix-1],ix,1,5);
      }
      // Upsilon
      for(unsigned int ix=3;ix<6;++ix) {
	for(unsigned int iy=1;iy<3;++iy) {
	  Histo1DPtr tmp;
	  _h_onium[ix-1].add(raps[2*iy-2],raps[2*iy],book(tmp,ix,1,iy));
	}
	book(_h_total[ix-1],ix,1,3);
      }
      // histos for ratios
      for(unsigned int ix=0;ix<5;++ix) {
	if(ix<2)
	  book(_h_ratio[ix],"TMP/h_ratio_"+to_str(ix),refData(11,1,1));
	else
	  book(_h_ratio[ix],"TMP/h_ratio_"+to_str(ix),refData(12,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::pid==553 or Cuts::pid==100553 or
					      Cuts::pid==200553 or
					      Cuts::pid==443 or Cuts::pid==100443)) {
	// cuts on pT and rapidity
	double y = p.absrap();
	double pT = p.perp();
	if(y>1.2 || pT<20.) continue;
	unsigned int itype(0);
	if(p.pid()==443)          itype=0;
	else if(p.pid()==100443)  itype=1;
	else if (p.pid()==553)    itype=2;
	else if (p.pid()==100553) itype=3;
	else if (p.pid()==200553) itype=4;
	// for J/pis and psi(2s) check if prompt
	if(itype<=1 && p.fromBottom()) continue;
	_h_onium[itype]. fill(y,pT);
	_h_total[itype]->fill(pT);
	_h_ratio[itype]->fill(pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 0.5 due folded rapidity
      double factor = 0.5*crossSection() / picobarn/ sumOfWeights();
      // branching ratios to muons
      vector<double> brs={0.05961,0.00793,0.0248,0.0191,0.0218};
      for(unsigned int ix=0;ix<5;++ix) {
	_h_onium[ix].scale(factor * brs[ix] , this);
	scale(_h_total[ix],brs[ix]*factor/1.2);
	scale(_h_ratio[ix],brs[ix]*factor/1.2);
      }
      // ratio plots
      Scatter2DPtr tmp;
      book(tmp,11,1,1);
      divide(_h_ratio[1],_h_ratio[0],tmp);
      for(unsigned int ix=0;ix<2;++ix) {
	book(tmp,12,1,1+ix);
	divide(_h_ratio[3+ix],_h_ratio[2],tmp);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistogram _h_onium[5];
    Histo1DPtr _h_total[5];
    Histo1DPtr _h_ratio[5];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CMS_2018_I1633431);

}
