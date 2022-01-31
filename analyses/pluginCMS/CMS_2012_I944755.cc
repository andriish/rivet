// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief J/psi and psi(2S) at 7 TeV
  class CMS_2012_I944755 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2012_I944755);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // J/psi/psi(2s) histos
      for(unsigned int ipsi=0;ipsi<2;++ipsi) {
	for(unsigned int iprompt=0;iprompt<2;++iprompt) {
	  book(_h_pT[ipsi][iprompt],"TMP/pT_"+toString(ipsi)+"_"+toString(iprompt),refData(20,1,1+iprompt));
	}
	vector<double> raps;
	if(ipsi==0) raps={0.,0.9,1.2,1.6,2.1,2.4};
	else        raps={0.,1.2,1.6,2.4};
	for(unsigned int iy=1;iy<raps.size();++iy) {
	  Histo1DPtr tmp;
	  for(unsigned int iprompt=0;iprompt<2;++iprompt) {
	    _h_psi[ipsi][iprompt].add(raps[iy-1],raps[iy],book(tmp,iy+5*ipsi,1,1+iprompt));
	    if(ipsi==1) {
	      _h_ratio[0][iprompt].add(raps[iy-1],raps[iy],book(tmp,"TMP/ratio_0_"+toString(iprompt)+"_"+toString(iy),refData(16+iy,1,1+iprompt)));
	      _h_ratio[1][iprompt].add(raps[iy-1],raps[iy],book(tmp,"TMP/ratio_1_"+toString(iprompt)+"_"+toString(iy),refData(16+iy,1,1+iprompt)));
	    }
	  }
	  _h_psi[ipsi][2].add(raps[iy-1],raps[iy],book(tmp,"TMP/psi_"+toString(ipsi)+"_"+toString(iy),refData(iy+5*ipsi,1,1)));
	}
      }
    }
    
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over onium states
      for( const Particle & p : ufs.particles(Cuts::pid==443 or Cuts::pid==100443)) {
	// cuts on rapidity
	double y = p.absrap();
	if(y>2.4) continue;
	double pT = p.perp();
	// prompt
	unsigned int iprompt = p.fromBottom();
	unsigned int ipsi=p.pid()/100000;
	_h_psi  [ipsi][iprompt].fill(y,pT);
	_h_psi  [ipsi][   2   ].fill(y,pT);
	_h_ratio[ipsi][iprompt].fill(y,pT);
	_h_pT[ipsi][iprompt]->fill(pT);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 0.5 due folded rapidity
      double factor = 0.5*crossSection() / nanobarn/ sumOfWeights();
      // branching ratios to muons
      vector<double> brs={0.05961,0.00793};
      for(unsigned int ipsi=0;ipsi<2;++ipsi) {
	for(unsigned int iprompt=0;iprompt<2;++iprompt) {
	  _h_psi  [ipsi][iprompt].scale(factor * brs[ipsi] , this);
	  _h_ratio[ipsi][iprompt].scale(factor * brs[ipsi] , this);
	  scale(_h_pT[ipsi][iprompt],factor * brs[ipsi]);
	}
	_h_psi  [ipsi][3].scale(factor * brs[ipsi] , this);
      }
      // ratios psi(2S) over J/psi
      for(unsigned int iprompt=0;iprompt<2;++iprompt) {
	for(unsigned int iy=0;iy<_h_ratio[0][iprompt].histos().size();++iy) {
	  Scatter2DPtr tmp;
	  book(tmp,17+iy,1,1+iprompt);
	  divide(_h_ratio[1][iprompt].histos()[iy],_h_ratio[0][iprompt].histos()[iy],tmp);
	}
	Scatter2DPtr tmp;
	book(tmp,20,1,1+iprompt);
	divide(_h_pT[1][iprompt],_h_pT[0][iprompt],tmp);
      }
      // non-prompt fraction
      for(unsigned int ipsi=0;ipsi<2;++ipsi) {
	for(unsigned int iy=0;iy<_h_psi[ipsi][1].histos().size();++iy) {
	  Scatter2DPtr tmp;
	  book(tmp,21+iy+5*ipsi,1,1);
	  efficiency(_h_psi[ipsi][1].histos()[iy],_h_psi[ipsi][2].histos()[iy],tmp);
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistogram _h_psi[2][3];
    BinnedHistogram _h_ratio[2][2];
    Histo1DPtr _h_pT[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2012_I944755);

}
