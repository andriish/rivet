// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {

  /// @brief inclusive W cross sections at 7 TeV
  class ATLAS_2016_I1502620_W : public Analysis {
  public:

    /// Constructor
    ATLAS_2016_I1502620_W(string name="ATLAS_2016_I1502620_W")
      : Analysis(name) {
      // using electron channel by default
      _mode = 0;
    }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      ///Initialise and register projections here
      const FinalState fs;

      Cut cut = Cuts::pT >= 25*GeV; // minimum lepton pT

      WFinder wfinder_dressed(fs, cut, _mode? PID::MUON : PID::ELECTRON, 40*GeV, 13*TeV, 25*GeV, 0.1, 
                              WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);

      declare(wfinder_dressed, "WFinder_dressed");

      /// Book histograms here
      _h_Wp_eta = bookHisto1D(   9, 1, _mode + 1);
      _h_Wm_eta = bookHisto1D(  10, 1, _mode + 1);
      _h_W_asym = bookScatter2D(35, 1, _mode + 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const WFinder& wfinder = apply<WFinder>(event, "WFinder_dressed");
      if (wfinder.bosons().size() != 1)  vetoEvent;

      const double weight = event.weight();
      const Particle lep = wfinder.constituentLeptons()[0];
      if (lep.charge3() > 0)  _h_Wp_eta->fill(lep.abseta(), weight);
      else                    _h_Wm_eta->fill(lep.abseta(), weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta)
      //divide(*_h_Wp_eta - *_h_Wm_eta, *_h_Wp_eta + *_h_Wm_eta, _h_W_asym);
      for (size_t i = 0; i < _h_Wp_eta->numBins(); ++i) {
        YODA::HistoBin1D& bp = _h_Wp_eta->bin(i);
        YODA::HistoBin1D& bm = _h_Wm_eta->bin(i);
        const double sum  = bp.height() + bm.height();
        //const double xerr = 0.5 * bp.xWidth();
        double val = 0., yerr = 0.;

        if (sum) {
          const double pos2  = bp.height() * bp.height();
          const double min2  = bm.height() * bm.height();
          const double errp2 = bp.heightErr() * bp.heightErr();
          const double errm2 = bm.heightErr() * bm.heightErr();
          val = (bp.height() - bm.height()) / sum;
          yerr = 2. * sqrt(errm2 * pos2 + errp2 * min2) / (sum * sum);
        }

        _h_W_asym->addPoint(bp.midpoint(), val, 0.5*bp.xWidth(), yerr);
      }

      // Print summary info
      const double xs_pb(crossSection() / picobarn);
      const double sumw(sumOfWeights());
      MSG_DEBUG( "Cross-section/pb     : " << xs_pb       );
      MSG_DEBUG( "Sum of weights       : " << sumw        );
      MSG_DEBUG( "nEvents              : " << numEvents() );

      ///  Normalise, scale and otherwise manipulate histograms here
      const double sf = 0.5 * xs_pb / sumw; // 0.5 accounts for rapidity bin width
      scale(_h_Wp_eta, sf);
      scale(_h_Wm_eta, sf);

    }

    //@}


  protected:
    size_t _mode;

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Wp_eta, _h_Wm_eta;
    Scatter2DPtr _h_W_asym;

    //@}

  };


  class ATLAS_2016_I1502620_W_EL : public ATLAS_2016_I1502620_W {
  public:
    ATLAS_2016_I1502620_W_EL()
      : ATLAS_2016_I1502620_W("ATLAS_2016_I1502620_W_EL")
    {
      _mode = 0;
    }
  };


  class ATLAS_2016_I1502620_W_MU : public ATLAS_2016_I1502620_W {
  public:
    ATLAS_2016_I1502620_W_MU()
      : ATLAS_2016_I1502620_W("ATLAS_2016_I1502620_W_MU")
    {
      _mode = 1;
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1502620_W);
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1502620_W_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1502620_W_MU);

}
