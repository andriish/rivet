// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_b production at 7,8 TeV
  class LHCB_2014_I1308738 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2014_I1308738);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      unsigned int iloc=0;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 1;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 2;
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 7 or 8 TeV.");
      // decays to 1S
      for(unsigned int ix=0;ix<3;++ix)
	book(_h_pT_1S[ix],"TMP/Ups_1_"+toString(ix),refData(3*iloc-2,1,1+ix));
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_pT_1S[ix+3],"TMP/Ups_1_"+toString(ix+3),refData(3*iloc-2,1,2+ix));
      // decays to 2S
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h_pT_2S[ix+2*iy],"TMP/Ups_2_"+toString(ix+2*iy),refData(3*iloc-1,1,2+ix));
      // decays to 3S
      for(unsigned int ix=0;ix<2;++ix)
	book(_h_pT_3S[ix],"TMP/Ups_3_"+toString(ix),refData(3*iloc,1,3));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      // first the upsilon states for the denominator
      for (const Particle& p : ufs.particles(Cuts::pid==553 ||
					     Cuts::pid==100553 ||
					     Cuts::pid==200553)) {
	// rapidity cut
	double absrap=p.absrap();
	if(absrap<2. || absrap>4.5) continue;
	double xp = p.perp();
	if(p.pid()==553) {
	  for(unsigned int ix=0;ix<2;++ix) _h_pT_1S[3+ix]->fill(xp);
	}
	else if (p.pid()==100553) {
	  for(unsigned int ix=0;ix<2;++ix) _h_pT_2S[2+ix]->fill(xp);
	}
	else if (p.pid()==200553) {
	  _h_pT_3S[1]->fill(xp);
	}
      }
      // P states
      for (const Particle& p : ufs.particles(Cuts::pid==10551 || Cuts::pid==110551 || Cuts::pid==210551 ||
      					     Cuts::pid==20553 || Cuts::pid==120553 || Cuts::pid==220553 ||
      					     Cuts::pid==  555 || Cuts::pid==100555 || Cuts::pid==200555 )) {
	Particle Upsilon;
	if     (p.children()[0].pid()==22 && abs(p.children()[1].pid())%100000==553) {
	  Upsilon=p.children()[1];
	}
	else if(p.children()[1].pid()==22 && abs(p.children()[0].pid())%100000==553) {
	  Upsilon=p.children()[0];
	}
	else
	  continue;
	double absrap=Upsilon.absrap();
	if(absrap<2. || absrap>4.5) continue;
	unsigned int iups = p.children()[0].pid()/100000;
	unsigned int ichi = p.pid()/100000;
	double xp  = Upsilon.perp();
	if(iups==0)
	  _h_pT_1S[ichi]->fill(xp);
	else if(iups==1 && ichi>0)
	  _h_pT_2S[ichi-1]->fill(xp);
	else if(iups==2 && ichi==2)
	  _h_pT_3S[0]->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      unsigned int iloc=0;
      if (isCompatibleWithSqrtS(7000)) {
	iloc = 1;
      }
      else if  (isCompatibleWithSqrtS(8000)) {
	iloc = 2;
      }
      // calculate the ratios
      // for 1S
      for(unsigned int iy=0;iy<3;++iy) {
	Scatter2DPtr tmp;
	book(tmp,3*iloc-2,1,iy+1);
	if(iy<2)
	  efficiency(_h_pT_1S[iy],_h_pT_1S[3],tmp);
	else
	  efficiency(_h_pT_1S[iy],_h_pT_1S[4],tmp);
	tmp->scaleY(100.);
      }
      Scatter2DPtr tmp;
      book(tmp,3*iloc-1,1,2);
      efficiency(_h_pT_2S[0],_h_pT_2S[2],tmp);
      tmp->scaleY(100.);
      book(tmp,3*iloc-1,1,3);
      efficiency(_h_pT_2S[1],_h_pT_2S[3],tmp);
      tmp->scaleY(100.);
      book(tmp,3*iloc,1,3);
      efficiency(_h_pT_3S[0],_h_pT_3S[1],tmp);
      tmp->scaleY(100.);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pT_1S[5],_h_pT_2S[4],_h_pT_3S[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2014_I1308738);

}
