// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief flavour seperate pi,K,p spectra
  class DELPHI_1998_I473409 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1998_I473409);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(InitialQuarks(), "IQF");

      // Book histograms
      _h_all_pi  = std::make_shared<YODA::Histo1D>(Histo1D(refData( 4,1,1)));
      _h_all_K   = std::make_shared<YODA::Histo1D>(Histo1D(refData( 5,1,1)));
      _h_all_p   = std::make_shared<YODA::Histo1D>(Histo1D(refData( 6,1,1)));
      _h_all_Kp  = std::make_shared<YODA::Histo1D>(Histo1D(refData( 7,1,1)));
      _d_all     = std::make_shared<YODA::Histo1D>(Histo1D(refData( 4,1,1)));
      
      _h_bot_pi  = std::make_shared<YODA::Histo1D>(Histo1D(refData( 8,1,1)));
      _h_bot_K   = std::make_shared<YODA::Histo1D>(Histo1D(refData( 9,1,1)));
      _h_bot_p   = std::make_shared<YODA::Histo1D>(Histo1D(refData(10,1,1)));
      _h_bot_Kp  = std::make_shared<YODA::Histo1D>(Histo1D(refData(11,1,1)));
      _d_bot     = std::make_shared<YODA::Histo1D>(Histo1D(refData( 8,1,1)));
      
      _h_lgt_pi  = std::make_shared<YODA::Histo1D>(Histo1D(refData(12,1,1)));
      _h_lgt_K   = std::make_shared<YODA::Histo1D>(Histo1D(refData(13,1,1)));
      _h_lgt_p   = std::make_shared<YODA::Histo1D>(Histo1D(refData(14,1,1)));
      _h_lgt_Kp  = std::make_shared<YODA::Histo1D>(Histo1D(refData(15,1,1)));
      _d_lgt     = std::make_shared<YODA::Histo1D>(Histo1D(refData(12,1,1)));

      _h_all_ch_p = bookHisto1D(16,1,1);
      _h_all_ch_x = bookHisto1D(17,1,1);
      _h_all_pi_p = bookHisto1D(18,1,1);
      _h_all_pi_x = bookHisto1D(19,1,1);
      _h_all_K_p  = bookHisto1D(20,1,1);
      _h_all_k_x  = bookHisto1D(21,1,1);
      _h_all_p_p  = bookHisto1D(22,1,1);
      _h_all_p_x  = bookHisto1D(23,1,1);

      _h_bot_ch_p = bookHisto1D(24,1,1);
      _h_bot_ch_x = bookHisto1D(25,1,1);
      _h_bot_pi_p = bookHisto1D(26,1,1);
      _h_bot_pi_x = bookHisto1D(27,1,1);
      _h_bot_K_p  = bookHisto1D(28,1,1);
      _h_bot_k_x  = bookHisto1D(29,1,1);
      _h_bot_p_p  = bookHisto1D(30,1,1);
      _h_bot_p_x  = bookHisto1D(31,1,1);

      _h_lgt_ch_p = bookHisto1D(32,1,1);
      _h_lgt_ch_x = bookHisto1D(33,1,1);
      _h_lgt_pi_p = bookHisto1D(34,1,1);
      _h_lgt_pi_x = bookHisto1D(35,1,1);
      _h_lgt_K_p  = bookHisto1D(36,1,1);
      _h_lgt_k_x  = bookHisto1D(37,1,1);
      _h_lgt_p_p  = bookHisto1D(38,1,1);
      _h_lgt_p_x  = bookHisto1D(39,1,1);

      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<5;++iy) {
	  ostringstream title;
	  title << "/TMP/MULT_" << ix << "_" << iy;
	  _mult[ix][iy] = bookCounter(title.str());
	}
      }
      _wLgt = _wBot = _wAll =0.; 
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<ChargedFinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");


      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      }
      else {
        map<int, double> quarkmap;
        foreach (const Particle& p, iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
           flavour = i;
          }
        }
      }

      // Get event weight for histo filling
      const double weight = event.weight();
      _wAll += weight;
      if(flavour<=3)      _wLgt += weight;
      else if(flavour==5) _wBot += weight;

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // loop over the charged particles
      foreach (const Particle& p, fs.particles()) {
	double modp = p.p3().mod();
        double xp = modp/meanBeamMom;
	int id = abs(p.pdgId());
	_d_all->fill(modp,weight);
	_mult[0][0]->fill(weight);
	_h_all_ch_p->fill(modp,weight);
	_h_all_ch_x->fill(xp  ,weight);
	if(flavour<=3) {
	  _d_lgt->fill(modp,weight);
	  _mult[2][0]->fill(weight);
	  _h_lgt_ch_p->fill(modp,weight);
	  _h_lgt_ch_x->fill(xp  ,weight);
	}
	else if(flavour==5) {
	  _d_bot  ->fill(modp,weight);
	  _mult[1][0]->fill(weight);
	  _h_bot_ch_p->fill(modp,weight);
	  _h_bot_ch_x->fill(xp  ,weight);
	}
	if(id==211) {
	  _h_all_pi ->fill(modp,weight);
	  _mult[0][1]->fill(weight);
	  _h_all_pi_p->fill(modp,weight);
	  _h_all_pi_x->fill(xp  ,weight);
	  if(flavour<=3) {
	    _h_lgt_pi ->fill(modp,weight); 
	    _mult[2][1]->fill(weight);
	    _h_lgt_pi_p->fill(modp,weight);
	    _h_lgt_pi_x->fill(xp  ,weight);
	  }
	  else if(flavour==5) {
	    _h_bot_pi ->fill(modp,weight);
	    _mult[1][1]->fill(weight);
	    _h_bot_pi_p->fill(modp,weight);
	    _h_bot_pi_x->fill(xp  ,weight);
	  }
	}
	else if(id==321) {
	  _h_all_K ->fill(modp,weight);
	  _h_all_Kp->fill(modp,weight);
	  _mult[0][2]->fill(weight);
	  _mult[0][4]->fill(weight);
	  _h_all_K_p ->fill(modp,weight);
	  _h_all_k_x ->fill(xp  ,weight);
	  if(flavour<=3) {
	    _h_lgt_K->fill(modp,weight);
	    _h_lgt_Kp->fill(modp,weight);
	    _mult[2][2]->fill(weight);
	    _mult[2][4]->fill(weight);
	    _h_lgt_K_p ->fill(modp,weight);
	    _h_lgt_k_x ->fill(xp  ,weight);
	  }
	  else if(flavour==5) {
	    _h_bot_K ->fill(modp,weight);
	    _h_bot_Kp->fill(modp,weight);
	    _mult[1][2]->fill(weight);
	    _mult[1][4]->fill(weight);
	    _h_bot_K_p ->fill(modp,weight);
	    _h_bot_k_x ->fill(xp  ,weight);
	  }
	}
	else if(id==2212) {
	  _h_all_p ->fill(modp,weight);
	  _h_all_Kp->fill(modp,weight);
	  _mult[0][3]->fill(weight);
	  _mult[0][4]->fill(weight);
	  _h_all_p_p ->fill(modp,weight);
	  _h_all_p_x ->fill(xp  ,weight);
	  if(flavour<=3) {
	    _h_lgt_p ->fill(modp,weight);
	    _h_lgt_Kp->fill(modp,weight);
	    _mult[2][3]->fill(weight);
	    _mult[2][4]->fill(weight);
	    _h_lgt_p_p ->fill(modp,weight);
	    _h_lgt_p_x ->fill(xp  ,weight);
	  }
	  else if(flavour==5) {
	    _h_bot_p ->fill(modp,weight);
	    _h_bot_Kp->fill(modp,weight); 
	    _mult[1][3]->fill(weight);
	    _mult[1][4]->fill(weight);
	    _h_bot_p_p ->fill(modp,weight);
	    _h_bot_p_x ->fill(xp  ,weight);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {


      // // Book histograms
      scale(_h_all_pi,100.);
      scale(_h_all_K ,100.);
      scale(_h_all_p ,100.);
      scale(_h_all_Kp,100.);
      divide(_h_all_pi, _d_all, bookScatter2D( 4,1,1));
      divide(_h_all_K , _d_all, bookScatter2D( 5,1,1));
      divide(_h_all_p , _d_all, bookScatter2D( 6,1,1));
      divide(_h_all_Kp, _d_all, bookScatter2D( 7,1,1));
      
      scale(_h_bot_pi,100.);
      scale(_h_bot_K ,100.);
      scale(_h_bot_p ,100.);
      scale(_h_bot_Kp,100.);
      divide(_h_bot_pi, _d_bot, bookScatter2D( 8,1,1));
      divide(_h_bot_K , _d_bot, bookScatter2D( 9,1,1));
      divide(_h_bot_p , _d_bot, bookScatter2D(10,1,1));
      divide(_h_bot_Kp, _d_bot, bookScatter2D(11,1,1));
      
      scale(_h_lgt_pi,100.);
      scale(_h_lgt_K ,100.);
      scale(_h_lgt_p ,100.);
      scale(_h_lgt_Kp,100.);
      divide(_h_lgt_pi, _d_lgt, bookScatter2D(12,1,1));
      divide(_h_lgt_K , _d_lgt, bookScatter2D(13,1,1));
      divide(_h_lgt_p , _d_lgt, bookScatter2D(14,1,1));
      divide(_h_lgt_Kp, _d_lgt, bookScatter2D(15,1,1));

      scale(_h_all_ch_p, 1./_wAll);
      scale(_h_all_ch_x, 1./_wAll);
      scale(_h_all_pi_p, 1./_wAll);
      scale(_h_all_pi_x, 1./_wAll);
      scale(_h_all_K_p , 1./_wAll);
      scale(_h_all_k_x , 1./_wAll);
      scale(_h_all_p_p , 1./_wAll);
      scale(_h_all_p_x , 1./_wAll);

      scale(_h_bot_ch_p, 1./_wBot);
      scale(_h_bot_ch_x, 1./_wBot);
      scale(_h_bot_pi_p, 1./_wBot);
      scale(_h_bot_pi_x, 1./_wBot);
      scale(_h_bot_K_p , 1./_wBot);
      scale(_h_bot_k_x , 1./_wBot);
      scale(_h_bot_p_p , 1./_wBot);
      scale(_h_bot_p_x , 1./_wBot);

      scale(_h_lgt_ch_p, 1./_wLgt);
      scale(_h_lgt_ch_x, 1./_wLgt);
      scale(_h_lgt_pi_p, 1./_wLgt);
      scale(_h_lgt_pi_x, 1./_wLgt);
      scale(_h_lgt_K_p , 1./_wLgt);
      scale(_h_lgt_k_x , 1./_wLgt);
      scale(_h_lgt_p_p , 1./_wLgt);
      scale(_h_lgt_p_x , 1./_wLgt);

      // multiplicities
      vector<double> scales = {_wAll,_wBot,_wLgt};
      for(unsigned int ix=0;ix<3;++ix) {
	if(scales[ix]<=0.) continue;
	for(unsigned int iy=0;iy<5;++iy) {
	  Scatter2DPtr scatter = bookScatter2D(ix+1, 1, iy+1, true);
	  scale(_mult[ix][iy],1./scales[ix]);
	  scatter->point(0).setY(_mult[ix][iy]->val(),_mult[ix][iy]->err());
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_all_pi , _h_all_K  , _h_all_p  , _h_all_Kp , _d_all;
    Histo1DPtr _h_bot_pi , _h_bot_K  , _h_bot_p  , _h_bot_Kp , _d_bot;
    Histo1DPtr _h_lgt_pi , _h_lgt_K  , _h_lgt_p  , _h_lgt_Kp , _d_lgt;
    Histo1DPtr _h_all_ch_p, _h_all_ch_x , _h_all_pi_p , _h_all_pi_x ;
    Histo1DPtr _h_all_K_p , _h_all_k_x  , _h_all_p_p  , _h_all_p_x  ;
    Histo1DPtr _h_bot_ch_p , _h_bot_ch_x , _h_bot_pi_p , _h_bot_pi_x;
    Histo1DPtr _h_bot_K_p  , _h_bot_k_x  , _h_bot_p_p  , _h_bot_p_x ;
    Histo1DPtr _h_lgt_ch_p , _h_lgt_ch_x , _h_lgt_pi_p , _h_lgt_pi_x;
    Histo1DPtr _h_lgt_K_p  , _h_lgt_k_x  , _h_lgt_p_p  , _h_lgt_p_x ;
    CounterPtr _mult[3][5];

    double _wLgt, _wBot, _wAll;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_1998_I473409);


}
