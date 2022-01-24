// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief Upsilon(1,2,3S) production at 7TeV
  class ATLAS_2013_I1204994 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1204994);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // branching ratios
      _br = {0.0248,0.0193,0.0218};
      // total cross section
      book(_h_total,1,1,1);
      // double differential
      vector<double> ybins = {0.,1.2,2.25};
      for(unsigned int ix=0;ix<3;++ix) {
      	book(_h_Upsilon_fid_y[ix],5,1,ix+1);
      	book(_h_Upsilon_y    [ix],9,1,ix+1);
      	book(_h_Upsilon_y_r  [ix],"TMP/Ups"+toString(ix),refData(12,1,1));
      	for(unsigned int iy=0;iy<2;++iy) {
      	  Histo1DPtr tmp;
      	  _h_Upsilon_fid[ix].add(ybins[iy],ybins[iy+1],book(tmp,ix+2,1,iy+1));
      	  _h_Upsilon    [ix].add(ybins[iy],ybins[iy+1],book(tmp,ix+6,1,iy+1));
      	  _h_Upsilon_r  [ix].add(ybins[iy],ybins[iy+1],book(tmp,"TMP/Ups"+toString(ix)+"_"+toString(iy), refData(11,1,iy+1)));
      	}
      }
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,  
                           Particles & mup, Particles & mum, unsigned int & ngamma) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::MUON ) {
          ++nstable;
    	  mum.push_back(p);
    	}
        else if (id == PID::ANTIMUON) {
          ++nstable;
    	  mup.push_back(p);
        }
        else if (id == PID::PI0 || id == PID::K0S || id == PID::K0L ) {
          ++nstable;
        }
    	else if (id == PID::GAMMA && p.children().empty() ) {
    	  ++ngamma;
    	  ++nstable;
    	}
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, mup, mum,ngamma);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::pid==553 or Cuts::pid==100553 or Cuts::pid==200553)) {
      	// check if in fiducal region
      	unsigned int nstable=0,ngamma=0;
      	Particles mup,mum;
      	findDecayProducts(p,nstable,mup,mum,ngamma);
      	bool fiduical = false;
      	if(mup.size()==1 && mum.size()==1 and nstable==ngamma+2) {
      	  fiduical = true;
      	  if(mup[0].perp()<4. || mup[0].abseta()>2.4) fiduical = false;
      	  if(mum[0].perp()<4. || mum[0].abseta()>2.4) fiduical = false;
      	}
      	// pT and rapidity
      	double absrap = p.absrap();
      	double xp = p.perp();
      	// type of upsilon
      	unsigned int iups=p.pid()/100000;
      	if(fiduical) {
      	  if(xp<70.) _h_Upsilon_fid_y[iups]->fill(absrap);
      	  _h_Upsilon_fid  [iups]. fill(absrap,xp);
      	}
      	if(xp<70.) {
	  _h_Upsilon_y[iups]->fill(absrap);
	  _h_Upsilon_y_r[iups]->fill(absrap);
	}
      	_h_Upsilon  [iups]. fill(absrap,xp);
      	_h_Upsilon_r[iups]. fill(absrap,xp);
      	if(absrap<2.25 && xp<70) _h_total->fill(iups+1,_br[iups]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_total,crossSection()/nanobarn/sumOfWeights());
      for(unsigned int ix=0;ix<3;++ix) {
      	scale(_h_Upsilon_fid_y[ix],0.5*crossSection()/picobarn/sumOfWeights());
      	scale(_h_Upsilon_y    [ix],0.5*_br[ix]*crossSection()/picobarn/sumOfWeights());
      	scale(_h_Upsilon_y_r  [ix],0.5*_br[ix]*crossSection()/picobarn/sumOfWeights());
      	_h_Upsilon_fid[ix].scale(0.5*crossSection()/femtobarn/sumOfWeights(),this);
      	_h_Upsilon    [ix].scale(0.5*_br[ix]*crossSection()/femtobarn/sumOfWeights(),this);
      	_h_Upsilon_r  [ix].scale(0.5*_br[ix]*crossSection()/femtobarn/sumOfWeights(),this);
      }
      // ratios
      for(unsigned int iy=0;iy<_h_Upsilon[0].histos().size();++iy) {
      	Scatter2DPtr tmp;
      	book(tmp,10,1,iy+1);
      	divide(_h_Upsilon_r[1].histos()[iy],_h_Upsilon_r[0].histos()[iy],tmp);
      	book(tmp,11,1,iy+1);
      	divide(_h_Upsilon_r[2].histos()[iy],_h_Upsilon_r[0].histos()[iy],tmp);
      }
      for(unsigned int iy=0;iy<2;++iy) {
      	Scatter2DPtr tmp;
      	book(tmp,12,1,iy+1);
      	divide(_h_Upsilon_y_r[iy+1],_h_Upsilon_y_r[0],tmp);
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_total,_h_Upsilon_fid_y[3],_h_Upsilon_y[3],_h_Upsilon_y_r[3];
    BinnedHistogram _h_Upsilon[3],_h_Upsilon_fid[3],_h_Upsilon_r[3];
    vector<double> _br;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1204994);

}
