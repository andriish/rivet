// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {

  /// @brief inclusive Z cross sections at 7 TeV
  class ATLAS_2016_I1502620_Z : public Analysis {
  public:

    /// Constructor
    ATLAS_2016_I1502620_Z(string name="ATLAS_2016_I1502620_Z")
      : Analysis(name) {
      // using electron channel by default
      _mode = 0;
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;

      Cut cuts = Cuts::pT >= 20.0*GeV;
      ZFinder zfinder(fs, cuts, (_mode ? PID::MUON : PID::ELECTRON), 46.0*GeV, 150*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfinder, "ZFinder");

      // Book histograms
      _h_Zcenlow_y_dressed   = bookHisto1D(11, 1, _mode + 1);
      _h_Zcenpeak_y_dressed  = bookHisto1D(12, 1, _mode + 1);
      _h_Zcenhigh_y_dressed  = bookHisto1D(13, 1, _mode + 1);
      _h_Zfwdpeak_y_dressed  = bookHisto1D(14, 1, _mode + 1);
      _h_Zfwdhigh_y_dressed  = bookHisto1D(15, 1, _mode + 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinder  = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1 ) vetoEvent;

      const Particle& Zboson = zfinder.boson();
      const double zrap  = Zboson.absrap();
      const double zmass = Zboson.mass();

      // Get/cut on Z boson leptons
      const ParticleVector& leptons = zfinder.constituents();
      const double eta1 = leptons[0].abseta();
      const double eta2 = leptons[1].abseta();

      const double weight = event.weight();

      // separation into central/forward and three mass bins
      if (eta1 < 2.5 && eta2 < 2.5) {
        if (zmass < 66.0*GeV)        _h_Zcenlow_y_dressed->fill(zrap, weight);
        else if (zmass < 116.0*GeV)  _h_Zcenpeak_y_dressed->fill(zrap, weight);
        else                         _h_Zcenhigh_y_dressed->fill(zrap, weight);
      } 
      else if ((eta1 < 2.5 && 2.5 < eta2 && eta2 < 4.9) || (eta2 < 2.5 && 2.5 < eta1 && eta1 < 4.9)) {
        if (zmass < 66.0*GeV)   vetoEvent;
        if (zmass < 116.0*GeV)  _h_Zfwdpeak_y_dressed->fill(zrap, weight);
        else                    _h_Zfwdhigh_y_dressed->fill(zrap, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Print summary info
      const double xs_pb(crossSection() / picobarn);
      const double sumw(sumOfWeights());
      MSG_DEBUG("Cross-Section/pb: " << xs_pb      );
      MSG_DEBUG("Sum of weights  : " << sumw       );
      MSG_DEBUG("nEvents         : " << numEvents());

      // Normalise, scale and otherwise manipulate histograms here
      const double sf(0.5 * xs_pb / sumw); // 0.5 accounts for rapidity bin width
      scale(_h_Zcenlow_y_dressed, sf);
      scale(_h_Zcenpeak_y_dressed, sf);
      scale(_h_Zcenhigh_y_dressed, sf);
      scale(_h_Zfwdpeak_y_dressed, sf);
      scale(_h_Zfwdhigh_y_dressed, sf);

    }

    //@}

  protected:
    size_t _mode;

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Zcenlow_y_dressed;
    Histo1DPtr _h_Zcenpeak_y_dressed;
    Histo1DPtr _h_Zcenhigh_y_dressed;
    Histo1DPtr _h_Zfwdpeak_y_dressed;
    Histo1DPtr _h_Zfwdhigh_y_dressed;
    //@}

  };


  class ATLAS_2016_I1502620_Z_EL : public ATLAS_2016_I1502620_Z {
  public:
    ATLAS_2016_I1502620_Z_EL()
      : ATLAS_2016_I1502620_Z("ATLAS_2016_I1502620_Z_EL")
    {
      _mode = 0;
    }
  };


  class ATLAS_2016_I1502620_Z_MU : public ATLAS_2016_I1502620_Z {
  public:
    ATLAS_2016_I1502620_Z_MU()
      : ATLAS_2016_I1502620_Z("ATLAS_2016_I1502620_Z_MU")
    {
      _mode = 1;
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1502620_Z);
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1502620_Z_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1502620_Z_MU);

}
