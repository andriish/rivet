// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {

  class ATLAS_2014_I1306294 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructors
    ATLAS_2014_I1306294(std::string name="ATLAS_2014_I1306294")
      : Analysis(name)
    {
      _mode = 1;
      setNeedsCrossSection(true);
    }

    //@}

  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;

      Cut cuts = Cuts::etaIn(-2.5,2.5) & (Cuts::pT > 20.0*GeV);

      ZFinder zfinder(fs, cuts, _mode==1? PID::ELECTRON : PID::MUON, 76.0*GeV, 106.0*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfinder, "ZFinder");

      //FastJets jetpro1( getProjection<ZFinder>("ZFinder").remainingFinalState(), FastJets::ANTIKT, 0.4);
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(getProjection<ZFinder>("ZFinder"));
      FastJets jetpro1(jet_fs, FastJets::ANTIKT, 0.4);
      jetpro1.useInvisibles();
      declare(jetpro1, "AntiKtJets04");
      declare(HeavyHadrons(), "BHadrons");

      //Histograms with data binning
      book(_h_bjet_Pt      , 3, 1, 1);
      book(_h_bjet_Y       , 5, 1, 1);
      book(_h_bjet_Yboost  , 7, 1, 1);
      book(_h_bjet_DY20    , 9, 1, 1);
      book(_h_bjet_ZdPhi20 ,11, 1, 1);
      book(_h_bjet_ZdR20   ,13, 1, 1);
      book(_h_bjet_ZPt     ,15, 1, 1);
      book(_h_bjet_ZY      ,17, 1, 1);
      book(_h_2bjet_dR     ,21, 1, 1);
      book(_h_2bjet_Mbb    ,23, 1, 1);
      book(_h_2bjet_ZPt    ,25, 1, 1);
      book(_h_2bjet_ZY     ,27, 1, 1);
    }

	  //==========================================================================================


    /// Perform the per-event analysis
    void analyze(const Event& e) {


      //---------------------------
      // -- check we have a Z:
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");

      if(zfinder.bosons().size() != 1)  vetoEvent;

      const ParticleVector boson_s =  zfinder.bosons();
      const Particle       boson_f =  boson_s[0]      ;
      const ParticleVector zleps   =  zfinder.constituents();
      //---------------------------


      //---------------------------
      //------------- stop processing the event if no true b-partons or hadrons are found
      const Particles& allBs = apply<HeavyHadrons>(e, "BHadrons").bHadrons(5.0*GeV);
      Particles stableBs;
      foreach(Particle p, allBs) {
        if(p.abseta() < 2.5)  stableBs += p;
      }
      if( stableBs.empty() )  vetoEvent;


      //---------------------------
      // -- get the b-jets:
      const Jets& jets = apply<JetAlg>(e, "AntiKtJets04").jetsByPt(Cuts::pT >20.0*GeV && Cuts::abseta <2.4);
      Jets b_jets;
      foreach(const Jet& jet, jets) {
        //veto overlaps with Z leptons:
        bool veto = false;
        foreach(const Particle& zlep, zleps) {
          if(deltaR(jet, zlep) < 0.5)  veto = true;
        }
        if(veto) continue;

        foreach(const Particle& bhadron, stableBs) {
          if( deltaR(jet, bhadron) <= 0.3 ) {
            b_jets.push_back(jet);
            break; // match
          }
	      } // end loop on b-hadrons
      }

      //and make sure we have at least 1:
      if(b_jets.empty())  vetoEvent;

      //---------------------------
      // fill the plots:
      const double ZpT = boson_f.pT()/GeV;
      const double ZY  = boson_f.absrap();

      _h_bjet_ZPt->fill(ZpT);
      _h_bjet_ZY ->fill(ZY);

      foreach(const Jet& jet, b_jets) {

        _h_bjet_Pt->fill(jet.pT()/GeV);
        _h_bjet_Y ->fill(jet.absrap());

        const double Yboost = 0.5 * fabs(boson_f.rapidity() + jet.rapidity());

        _h_bjet_Yboost->fill(Yboost);

        if(ZpT > 20.) {

          const double ZBDY   = fabs( boson_f.rapidity() - jet.rapidity() );
          const double ZBDPHI = fabs( deltaPhi(jet.phi(), boson_f.phi()) );
          const double ZBDR   = deltaR(jet, boson_f, RAPIDITY);
          _h_bjet_DY20->fill(   ZBDY);
          _h_bjet_ZdPhi20->fill(ZBDPHI);
          _h_bjet_ZdR20->fill(  ZBDR);
        }

      } //loop over b-jets

      if (b_jets.size() < 2) return;

      _h_2bjet_ZPt->fill(ZpT);
      _h_2bjet_ZY ->fill(ZY);

      const double BBDR = deltaR(b_jets[0], b_jets[1], RAPIDITY);
      const double Mbb  = (b_jets[0].momentum() + b_jets[1].momentum()).mass();

      _h_2bjet_dR ->fill(BBDR);
      _h_2bjet_Mbb->fill(Mbb);

    } // end of analysis loop


    /// Normalise histograms etc., after the run
    void finalize() {

      const double normfac = crossSection() / sumOfWeights();

      scale( _h_bjet_Pt,      normfac);
      scale( _h_bjet_Y,       normfac);
      scale( _h_bjet_Yboost,  normfac);
      scale( _h_bjet_DY20,    normfac);
      scale( _h_bjet_ZdPhi20, normfac);
      scale( _h_bjet_ZdR20,   normfac);
      scale( _h_bjet_ZPt,     normfac);
      scale( _h_bjet_ZY,      normfac);
      scale( _h_2bjet_dR,     normfac);
      scale( _h_2bjet_Mbb,    normfac);
      scale( _h_2bjet_ZPt,    normfac);
      scale( _h_2bjet_ZY,     normfac);
    }

    //@}


  protected:

    // Data members like post-cuts event weight counters go here
    size_t _mode;


  private:

    Histo1DPtr _h_bjet_Pt;
    Histo1DPtr _h_bjet_Y;
    Histo1DPtr _h_bjet_Yboost;
    Histo1DPtr _h_bjet_DY20;
    Histo1DPtr _h_bjet_ZdPhi20;
    Histo1DPtr _h_bjet_ZdR20;
    Histo1DPtr _h_bjet_ZPt;
    Histo1DPtr _h_bjet_ZY;
    Histo1DPtr _h_2bjet_dR;
    Histo1DPtr _h_2bjet_Mbb;
    Histo1DPtr _h_2bjet_ZPt;
    Histo1DPtr _h_2bjet_ZY;

  };


  class ATLAS_2014_I1306294_EL : public ATLAS_2014_I1306294 {
  public:
    ATLAS_2014_I1306294_EL()
      : ATLAS_2014_I1306294("ATLAS_2014_I1306294_EL")
    {
      _mode = 1;
    }
  };

  class ATLAS_2014_I1306294_MU : public ATLAS_2014_I1306294 {
  public:
    ATLAS_2014_I1306294_MU()
      : ATLAS_2014_I1306294("ATLAS_2014_I1306294_MU")
    {
      _mode = 2;
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1306294);
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1306294_MU);
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1306294_EL);

}
