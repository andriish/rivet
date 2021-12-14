// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Upsilon(1,2,3S) at 7 TeV
  class CMS_2013_I1225274 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2013_I1225274);

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_total[0],1,1,1);
      book(_h_total[1],1,1,2);
      vector<double> raps= {0.,0.4,0.8,1.2,1.6,2.0,2.4};
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  for(unsigned int iz=0;iz<6;++iz) {
	    Histo1DPtr tmp;
	    _h_pT_y[ix][iy].add(raps[iz],raps[iz+1], book(tmp,5+3*iz+ix,1,iy+1));
	  }
	  book(_h_pT[ix][iy],2 +ix,1,iy+1);
	  book(_h_y [ix][iy],23+ix,1,iy+1);
	  book(_h_r [ix][iy],"TMP/h_r_"+toString(ix)+"_"+toString(iy),refData(26,1,iy+1));
	}
	book(_h_pT_acc[ix],29+ix,1,1);
      }
      // branching ratios
      _br={0.0248,0.0193,0.0218};
    }

    void findChildren(const Particle & p,Particles & mum, Particles & mup, unsigned int & nstable) {
      for(const Particle & child : p.children()) {
	if(child.pid()==PID::MUON) {
	  mum.push_back(child);
	  ++nstable;
	}
	else if(child.pid()==PID::ANTIMUON) {
	  mup.push_back(child);
	  ++nstable;
	}
	else if(child.pid()==PID::PHOTON)
	  continue;
	else if(child.children().empty()) {
	  ++nstable;
	}
	else
	  findChildren(child,mum,mup,nstable);
      }
    }

    // from eqn 1 of paper
    bool acceptMuon(const Particle & p) {
      double abseta = p.abseta();
      double xp = p.perp();
      if(abseta<0.8)      return xp>3.75;
      else if(abseta<1.6) return xp>3.5;
      else if(abseta<2.4) return xp>3.0;
      else return false;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==553 || Cuts::pid==100553 || Cuts::pid==200553)) {
        double absrap = p.absrap();
	// rapidity cut
	if(absrap>2.4) continue;
        double xp = p.perp();
	unsigned int iloc=0;
	if     (p.pid()==   553) iloc=0;
	else if(p.pid()==100553) iloc=1;
	else if(p.pid()==200553) iloc=2;
	// acceptance corrected only hist
	if(absrap<1.2)
	  _h_pT_acc[iloc]->fill(xp);
	// check if children muons and within acceptance
	unsigned int imin=1;
	// find the children
	Particles mum,mup;
	unsigned int nstable(0);
	findChildren(p,mum,mup,nstable);
	if(mup.size()==1 && mup.size()==1 && nstable==2) {
	  if(acceptMuon(mup[0]) && acceptMuon(mum[0])) imin=0;
	}
	// fill the histos
	for(unsigned int ix=imin;ix<2;++ix) {
	  _h_pT_y[iloc][ix] .fill(absrap,xp);
	  _h_pT  [iloc][ix]->fill(xp);
	  _h_r   [iloc][ix]->fill(xp);
	  if(xp<50.) _h_y[iloc][ix]->fill(absrap);
	  if(ix==0)
	    _h_total[ix]->fill(iloc+1);
	  else
	    _h_total[ix]->fill(iloc+1,_br[iloc]);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double factor = crossSection() / nanobarn/ sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) {
	// total cross sections, just the factor
	scale(_h_total[ix],factor);
	for(unsigned int iy=0;iy<3;++iy) {
	  double factor2 = factor;
	  if(ix==1) factor2*=_br[iy];
	  // pT integrated over y, just the factor
	  scale(_h_pT    [iy][ix],factor2);
	  scale(_h_r     [iy][ix],factor2);
	  if(ix==1) scale(_h_pT_acc[iy],factor2);
	  // not integrated over y, alsso undo y +/- folding
	  scale(_h_y[iy][ix],0.5*factor2);
	  _h_pT_y[iy][ix].scale(0.5*factor2,this);
	}
	// ratios
	Scatter2DPtr tmp;
	// ups 3/ ups1
	book(tmp,26,1,ix+1);
	divide(_h_r[2][ix],_h_r[0][ix],tmp);
	// ups 2/ ups1
	book(tmp,27,1,ix+1);
	divide(_h_r[1][ix],_h_r[0][ix],tmp);
	// ups 3/ ups2
	book(tmp,28,1,ix+1);
	divide(_h_r[2][ix],_h_r[1][ix],tmp);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistogram _h_pT_y[3][2];
    Histo1DPtr _h_total[2],_h_pT[3][2],_h_pT_acc[3],_h_y[3][2],_h_r[3][2];
    vector<double> _br;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CMS_2013_I1225274);

}
