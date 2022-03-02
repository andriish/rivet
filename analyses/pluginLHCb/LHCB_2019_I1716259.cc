// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_b- 7,8 and 13 TeV
  class LHCB_2019_I1716259 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2019_I1716259);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // CMS energy
      if (isCompatibleWithSqrtS(7000) ||
	  isCompatibleWithSqrtS(8000)) {
	_rootS=7500;
      }
      else if (isCompatibleWithSqrtS(13000)) {
	_rootS=13000;
      }
      else
	throw UserError("Centre-of-mass energy of the given input is neither 7, 8, 13 TeV.");
      // histograms
      book(_h_xi ,"TMP/h_xi" ,refData(1,1,1));
      book(_h_lam,"TMP/h_lam",refData(1,1,1));
      for(unsigned int ix=0;ix<2;++ix) {
	book(_c_xi[ix] ,"TMP/c_xi_" +toString(ix+1));
	book(_c_lam[ix],"TMP/c_lam_"+toString(ix+1));
      }
    }
    
    void findDecayProducts(Particle mother, double sign, Particles & Lambda, Particles & Xi, Particles & Jpsi, unsigned int & nstable) {
      for(const Particle & p: mother.children()) {
	if(p.pid() == 3122*sign)
	  Lambda.push_back(p);
	else if(p.pid() == 3312*sign)
	  Xi.push_back(p);
	else if(p.pid()==443)
	  Jpsi.push_back(p);
	else if(p.pid()==111 || p.children().empty())
	  ++nstable;
	else 
	  findDecayProducts(p,sign,Lambda,Xi,Jpsi,nstable);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==5122 || Cuts::abspid==5132)) {
	// decay modes
	Particles Lambda,Xi,Jpsi;
	unsigned int nstable=0;
	double sign = p.pid()>0 ? 1. : -1.;
	findDecayProducts(p,sign,Lambda,Xi,Jpsi,nstable);
	if (p.abspid()==5112) {
	  _c_lam[1]->fill();
	  if(Lambda.size()==1 && Jpsi.size()==1 && nstable==0)
	    _c_lam[0]->fill();
	}
	else {
	  _c_xi[1]->fill();
	  if(Xi.size()==1 && Jpsi.size()==1 && nstable==0)
	    _c_xi[0]->fill();
	}
	// pT and rapidity cuts
	if(p.perp()>20.) continue;
	double eta=p.abseta();
	if(eta<2. || eta> 6.) continue;
	if (p.abspid()==5122) {
	  _h_lam->fill(_rootS);
	}
	else {
	  _h_xi->fill(_rootS);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // first the simple ratio
      Scatter2DPtr tmp;
      book(tmp,2,1,1);
      divide(_h_xi,_h_lam,tmp);
      // and the one with brs included
      book(tmp,1,1,1);
      if (_c_xi[1]->effNumEntries()>0.)
      	scale(_h_xi, *_c_xi[0]/ *_c_xi[1]);
      if (_c_lam[1]->effNumEntries()>0.)
      	scale(_h_lam, *_c_lam[0]/ *_c_lam[1]);
      divide(_h_xi,_h_lam,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_xi,_h_lam;
    CounterPtr _c_xi[2],_c_lam[2];
    double _rootS;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2019_I1716259);

}
