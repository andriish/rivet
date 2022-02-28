// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief xi_cc++ production
  class LHCB_2020_I1760788 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2020_I1760788);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_n_lambda_c,"TMP/nLambdac");
      for(unsigned int ix=0;ix<2;++ix)
	book(_n_Xicc[ix],"TMP/nXicc_"+toString(ix));
      book(_h_xi,1,1,1);
    }
    
    void findDecayProducts(Particle mother, double sign, Particles & lambdac, Particles & Km, Particles & pip, unsigned int & nstable) {
      for(const Particle & p: mother.children()) {
	if(p.pid()==4122*sign)
	  lambdac.push_back(p);
	else if(p.pid()==-321*sign)
	  Km.push_back(p);
	else if(p.pid()==211*sign)
	  pip.push_back(p);
	else if(p.pid()==111 || p.children().empty())
	  ++nstable;
	else 
	  findDecayProducts(p,sign,lambdac,Km,pip,nstable);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==4122 || Cuts::abspid==4422)) {
	if(p.abspid()==4422) {
	  Particles lambdac,Km,pip;
	  unsigned int nstable=0;
	  double sign = p.pid()>0 ? 1. : -1.;
	  findDecayProducts(p,sign,lambdac,Km,pip,nstable);
	  if(lambdac.size()==1 && Km.size()==1 && pip.size()==2 && nstable ==0)
	    _n_Xicc[0]->fill();
	  _n_Xicc[1]->fill();
	}
	// pT and rapidity cuts
	double pT=p.perp();
	if(pT<4. || pT>15.) continue;
	double y=p.absrap();
	if(y<2. || y> 4.5) continue;
	if (p.abspid()==4422) {
	  for(unsigned int ix=1;ix<4;++ix)
	    _h_xi->fill(ix);
	}
	else
	  _n_lambda_c->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // scale by Xi_cc br
      if (_n_Xicc[1]->effNumEntries()>0.)
      	scale(_h_xi, *_n_Xicc[0]/ *_n_Xicc[1]);
      // and by lambda_c rate for normalisation (10^6 as rate in units 10^-4)
      if (_n_lambda_c->effNumEntries()>0.)
	scale(_h_xi,1e4/ *_n_lambda_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_xi;
    CounterPtr _n_lambda_c, _n_Xicc[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2020_I1760788);

}
