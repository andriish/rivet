// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/LossyFinalState.hh"

namespace Rivet {


  /// @todo Replace with SmearedParticles
  class STARRandomFilter {
  public:

    STARRandomFilter() { }

    // Return true to throw away a particle
    bool operator()(const Particle& p) {
      /// @todo Use a better RNG?
      size_t idx = int(floor(p.pT()/MeV/50));
      if (idx > 11) idx = 11;
      return (rand()/static_cast<double>(RAND_MAX) > _trkeff[idx]);
    }

    CmpState compare(const STARRandomFilter& other) const {
      return CmpState::GT; // @todo really?
    }

  private:

    const static double _trkeff[12];

  };


  // Here we have the track reconstruction efficiencies for tracks with pT from 0 to 600 MeV
  // in steps of 50 MeV. The efficiency is assumed to be 0.88 for pT >= 600 MeV
  const double STARRandomFilter::_trkeff[12] = {0,0,0.38,0.72,0.78,0.81,0.82,0.84,0.85,0.86,0.87,0.88};



  class STAR_2008_S7869363 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    STAR_2008_S7869363()
      : Analysis("STAR_2008_S7869363")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs((Cuts::etaIn(-0.5, 0.5) && Cuts::pT >=  0.2*GeV));
      const LossyFinalState<STARRandomFilter> lfs(cfs, STARRandomFilter());
      declare(lfs, "FS");

      book(_h_dNch           ,1, 1, 1);
      book(_h_dpT_Pi         ,2, 1, 1);
      book(_h_dpT_Piplus     ,2, 1, 2);
      book(_h_dpT_Kaon       ,2, 1, 3);
      book(_h_dpT_Kaonplus   ,2, 1, 4);
      book(_h_dpT_AntiProton ,2, 1, 5);
      book(_h_dpT_Proton     ,2, 1, 6);
      // book(nCutsPassed, "nCutsPassed");
      // book(nPi, "nPi");
      // book(nPiPlus, "nPiPlus");
      // book(nKaon, "nKaon");
      // book(nKaonPlus, "nKaonPlus");
      // book(nProton, "nProton");
      // book(nAntiProton, "nAntiProton");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& charged = apply<FinalState>(event, "FS");

      // Vertex reconstruction efficiencies as a function of charged multiplicity.
      // For events with more than 23 reconstructed tracks the efficiency is 100%.
      double vtxeffs[24] = { 0.000000,0.512667,0.739365,0.847131,0.906946,0.940922,0.959328,0.96997,
                             0.975838,0.984432,0.988311,0.990327,0.990758,0.995767,0.99412,0.992271,
                             0.996631,0.994802,0.99635,0.997384,0.998986,0.996441,0.994513,1.000000 };

      double vtxeff = 1.0;
      if (charged.particles().size() < 24) {
        vtxeff = vtxeffs[charged.particles().size()];
      }

      const double weight = vtxeff;

      for (const Particle& p : charged.particles()) {
        double pT = p.pT()/GeV;
        double y = p.rapidity();
        if (fabs(y) < 0.1) {
          // nCutsPassed->fill(weight);
          const PdgId id = p.pid();
          switch (id) {
          case -211:
            _h_dpT_Pi->fill(pT, weight/(TWOPI*pT*0.2));
            // nPi->fill(weight);
            break;
          case 211:
            _h_dpT_Piplus->fill(pT, weight/(TWOPI*pT*0.2));
            // nPiPlus->fill(weight);
            break;
          case -321:
            _h_dpT_Kaon->fill(pT, weight/(TWOPI*pT*0.2));
            // nKaon->fill(weight);
            break;
          case 321:
            _h_dpT_Kaonplus->fill(pT, weight/(TWOPI*pT*0.2));
            // nKaonPlus->fill(weight);
            break;
          case -2212:
            _h_dpT_AntiProton->fill(pT, weight/(TWOPI*pT*0.2));
            // nAntiProton->fill(weight);
            break;
          case 2212:
            _h_dpT_Proton->fill(pT, weight/(TWOPI*pT*0.2));
            // nProton->fill(weight);
            break;
          }
        }
        else {
          continue;
        }
      }
      _h_dNch->fill(charged.particles().size(), weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //double nTot = nPi + nPiPlus + nKaon + nKaonPlus + nProton + nAntiProton;
      normalize(_h_dNch);

      /// @todo Norm to data!
      normalize(_h_dpT_Pi        , 0.389825 );
      normalize(_h_dpT_Piplus    , 0.396025 );
      normalize(_h_dpT_Kaon      , 0.03897  );
      normalize(_h_dpT_Kaonplus  , 0.04046  );
      normalize(_h_dpT_AntiProton, 0.0187255);
      normalize(_h_dpT_Proton    , 0.016511 );
    }


  private:


    Histo1DPtr _h_dNch;

    Histo1DPtr _h_dpT_Pi, _h_dpT_Piplus;
    Histo1DPtr _h_dpT_Kaon, _h_dpT_Kaonplus;
    Histo1DPtr _h_dpT_AntiProton, _h_dpT_Proton;

    Profile1DPtr _h_pT_vs_Nch;
    //CounterPtr nCutsPassed, nPi, nPiPlus, nKaon, nKaonPlus, nProton, nAntiProton;
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_2008_S7869363);

}
