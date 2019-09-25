// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/GammaGammaFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class OPAL_2003_I611415 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_2003_I611415);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // get the hadronic final state
      const GammaGammaKinematics& diskin = declare(GammaGammaKinematics(), "Kinematics");
      const FinalState & fs = declare(GammaGammaFinalState(diskin, GammaGammaFinalState::LAB), "FS");
      declare(FastJets(fs, FastJets::KT,1.),"Jets");
      _h_theta[0]    = bookHisto1D( 1,1,1);
      _h_theta[1]    = bookHisto1D( 2,1,1);
      _h_ET[0]       = bookHisto1D( 3,1,1);
      _h_ET[1]       = bookHisto1D( 4,1,1);
      _h_ET[2]       = bookHisto1D( 5,1,1);
      _h_xg[0][0]    = bookHisto1D( 6,1,1);
      _h_xg[0][1]    = bookHisto1D( 7,1,1);
      _h_xg[1][0]    = bookHisto1D( 9,1,1);
      _h_xg[1][1]    = bookHisto1D(10,1,1);
      _h_xg[2][0]    = bookHisto1D(11,1,1);
      _h_xg[2][1]    = bookHisto1D(12,1,1);
      _h_xg_high     = bookHisto1D( 8,1,1);
      _h_xlog[0]     = bookHisto1D(13,1,1);
      _h_xlog[1]     = bookHisto1D(14,1,1);
      _h_xlog[2]     = bookHisto1D(15,1,1);
      _h_eta_diff[0] = bookHisto1D(16,1,1);
      _h_eta_diff[1] = bookHisto1D(17,1,1);
      _h_eta_min[0]  = bookHisto1D(18,1,1);
      _h_eta_min[1]  = bookHisto1D(19,1,1);
      _h_eta_max[0]  = bookHisto1D(20,1,1);
      _h_eta_max[1]  = bookHisto1D(21,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double weight = event.weight();
      // need at least two jets with |eta|<2 and pT>3 
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::Et > 3.*GeV and Cuts::abseta < 2.);
      if(jets.size()<2) vetoEvent;
      if(jets[0].Et()<jets[1].Et()) swap(jets[0],jets[1]);
      // Ets of jets
      double Et1 = jets[0].Et(), Et2 = jets[1].Et();
      // average Et
      double Etbar = 0.5*(Et1+Et2);
      double etaBar = 0.5*(jets[0].eta()+jets[1].eta());
      if(Etbar<5.) vetoEvent;
      // assymetry cut
      if((Et1-Et2)/(Et1+Et2)>.25) vetoEvent;
      // calculate x_gamma
      FourMomentum psum;
      for(const Particle & part : apply<FinalState>(event,"FS").particles()) {
	psum += part.momentum();
      }
      FourMomentum pj = jets[0].momentum()+jets[1].momentum();
      double xp = (pj.E()+pj.pz())/(psum.E()+psum.pz());
      double xm = (pj.E()-pj.pz())/(psum.E()-psum.pz());
      double cost = tanh(0.5*(jets[0].eta()-jets[1].eta()));
      // cost distributions
      if(pj.mass()>15.*GeV && etaBar<=1.) {
	if(xp>0.75 && xm>0.75)
	  _h_theta[0]->fill(abs(cost),weight);
	else if(xp<0.75 && xm<0.75)
	  _h_theta[1]->fill(abs(cost),weight);
      }
      // ET distributions
      _h_ET[0]->fill(Etbar,weight);
      if((xp<0.75 && xm>0.75)|| (xm<0.75&&xp>0.75))
	_h_ET[1]->fill(Etbar,weight);
      else if(xp<0.75 && xm <0.75)
	_h_ET[2]->fill(Etbar,weight);
      if(Etbar>=5.&&Etbar<7.) {
	_h_xg[0][0]->fill(xp,weight);
	_h_xg[0][0]->fill(xm,weight);
	_h_xlog[0]->fill(log(xp),weight);
	_h_xlog[0]->fill(log(xm),weight);
	if((xp<0.75 && xm>0.75)|| (xm<0.75&&xp>0.75)) {
	  _h_xg[1][0]->fill(xp,weight);
	  _h_xg[1][0]->fill(xm,weight);
	  _h_xlog[1]->fill(log(xp),weight);
	  _h_xlog[1]->fill(log(xm),weight);
	}
	else if(xp<0.75 && xm <0.75) {
	  _h_xg[2][0]->fill(xp,weight);
	  _h_xg[2][0]->fill(xm,weight);
	  _h_xlog[2]->fill(log(xp),weight);
	  _h_xlog[2]->fill(log(xm),weight);
	}
      }
      else if(Etbar>=7.&& Etbar<11.) {
	_h_xg[0][1]->fill(xp,weight);
	_h_xg[0][1]->fill(xm,weight);
	if((xp<0.75 && xm>0.75)|| (xm<0.75&&xp>0.75)) {
	  _h_xg[1][1]->fill(xp,weight);
	  _h_xg[1][1]->fill(xm,weight);
	}
	else if(xp<0.75 && xm <0.75) {
	  _h_xg[2][1]->fill(xp,weight);
	  _h_xg[2][1]->fill(xm,weight);
	}
      }
      else if(Etbar>=11.&& Etbar<25.) {
	_h_xg_high->fill(xp,weight);
	_h_xg_high->fill(xm,weight);
      }
      // vs eta
      double etaMin = min(abs(jets[0].eta()),abs(jets[1].eta()));
      double etaMax = max(abs(jets[0].eta()),abs(jets[1].eta()));
      if((xp<0.75 && xm>0.75)|| (xm<0.75&&xp>0.75)) {
	_h_eta_diff[0]->fill(abs(jets[0].eta()-jets[1].eta()),weight);
	_h_eta_min[0]->fill(etaMin,weight);
	_h_eta_max[0]->fill(etaMax,weight);
      }
      else if(xp<0.75 && xm <0.75) {
	_h_eta_diff[1]->fill(abs(jets[0].eta()-jets[1].eta()),weight);
	_h_eta_min[1]->fill(etaMin,weight);
	_h_eta_max[1]->fill(etaMax,weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      for(unsigned int ix=0;ix<2;++ix) {
	scale(_h_theta[ix], fact);
	scale(_h_eta_diff[ix], fact);
	scale(_h_eta_min[ix], fact);
	scale(_h_eta_max[ix], fact);
	for(unsigned int iy=0;iy<3;++iy) {
	  scale(_h_xg[iy][ix],fact);
	}
      }
      for(unsigned int ix=0;ix<3;++ix) {
	scale(_h_ET[ix],fact);
	scale(_h_xlog[ix],fact);
      }
      scale(_h_xg_high,fact);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_theta[2],_h_ET[3],_h_xg[3][2],_h_xg_high;
    Histo1DPtr _h_xlog[3],_h_eta_diff[2],_h_eta_min[2],_h_eta_max[2];
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2003_I611415);


}
