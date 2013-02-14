// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetYODA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  class ALICE_2012_I1181770 : public Analysis {
  public:

    ALICE_2012_I1181770()
      : Analysis("ALICE_2012_I1181770")
    {    }

  public:

    void init() {
      addProjection(ChargedFinalState(),"CFS");

      if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3)) {
        _h_frac_sd_inel      = bookScatter2D(1, 1, 1);
        _h_frac_dd_inel      = bookScatter2D(2, 1, 1);
        _h_xsec_sd           = bookHisto1D  (3, 1, 1);
        _h_xsec_dd           = bookHisto1D  (4, 1, 1);
        _h_xsec_inel         = bookHisto1D  (5, 1, 1);
      } else if (fuzzyEquals(sqrtS()/GeV, 2760, 1E-3)) {
        _h_frac_sd_inel      = bookScatter2D(1, 1, 2);
        _h_frac_dd_inel      = bookScatter2D(2, 1, 2);
        _h_xsec_sd           = bookHisto1D  (3, 1, 2);
        _h_xsec_dd           = bookHisto1D  (4, 1, 2);
        _h_xsec_inel         = bookHisto1D  (5, 1, 2);
      } else if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) {
        _h_frac_sd_inel      = bookScatter2D(1, 1, 3);
        _h_frac_dd_inel      = bookScatter2D(2, 1, 3);
        _h_xsec_sd           = bookHisto1D  (3, 1, 3);
        _h_xsec_dd           = bookHisto1D  (4, 1, 3);
        _h_xsec_inel         = bookHisto1D  (5, 1, 3);
      }
    }

    void analyze(const Event& event) {
      const double weight = event.weight();

      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");

      // fill INEL plots for each event
      if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3)) {
        _h_xsec_inel->fill( 900/GeV, weight );
      } else if (fuzzyEquals(sqrtS()/GeV, 2760, 1E-3)) {
        _h_xsec_inel->fill( 2760/GeV, weight );
      } else if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) {
        _h_xsec_inel->fill( 7000/GeV, weight );
      }

      double Eslowest(0.0), pslowest(0.0), yslowest(999.);
      double Efastest(0.0), pfastest(0.0), yfastest(-999.);
      int pidslowest(0), pidfastest(0);

      double LRG(0.0), etapre(0.0), gapbwd(0.0), gapfwd(0.0);
      unsigned int num_p(0);

      FourMomentum leadP(0.,0.,0.,0.);
      double Elead(0.0), plead(0.0);

      foreach(const Particle& p, cfs.particlesByEta()) { //sorted from minus to plus
        const PdgId pid = p.pdgId();
        double y = p.momentum().rapidity();
        double eta = p.momentum().eta();
        //SD case
        if (y < yslowest) {
          Eslowest = p.momentum().E();
          pslowest = p.momentum().vector3().mod();
          yslowest = p.momentum().rapidity();
          pidslowest = pid;
        }
        if (y > yfastest) {
          Efastest = p.momentum().E();
          pfastest = p.momentum().vector3().mod();
          yfastest = p.momentum().rapidity();
          pidfastest = pid;
        }

        num_p += 1;
        // DD case
        if (num_p==1) {
          etapre = p.momentum().eta();
        } else if (num_p > 1) {
          if (num_p==2) gapbwd = fabs(eta-etapre);

          double gap = fabs(eta-etapre);
          LRG = (gap > LRG ? gap : LRG); // largest gap

          if (num_p==cfs.size()) gapfwd = fabs(eta-etapre);
          etapre = eta;
        }
      }

      // Mx calculation
      if (pidslowest==2212 && pidfastest==2212) {
        if (fabs(yslowest) > fabs(yfastest)) {
          Elead = Eslowest;
          plead = pslowest;
        } else if (fabs(yslowest) < fabs(yfastest)) {
          Elead = Efastest;
          plead = pfastest;
        } else {
          Elead = Eslowest;
          plead = pslowest;
          if ( (double)rand() / (double)RAND_MAX > 0.5) { // generate random number in [0.,1.] range and make decision randomly
            Elead = Efastest;
            plead = pfastest;
          }
        }
      } else if (pidslowest==2212) {
        Elead = Eslowest;
        plead = pslowest;
      } else if (pidfastest==2212) {
        Elead = Efastest;
        plead = pfastest;
      }

      double Mx = sqrt((sqrtS()/GeV-Elead-plead)*(sqrtS()/GeV-Elead+plead));
      bool singleDiff = false;

      // Fill SD
      if (Mx < 200.) {
        singleDiff = true;
        if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3)) {
          _h_xsec_sd->fill( 900/GeV, weight);
        } else if (fuzzyEquals(sqrtS()/GeV, 2760, 1E-3)) {
          _h_xsec_sd->fill( 2760/GeV, weight);
        } else if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) {
          _h_xsec_sd->fill( 7000/GeV, weight);
        }
      }

      if ( singleDiff ) vetoEvent; // DD events are defined as NSD with large gap.

      // also remove SD-like events in NSD events
      if ( std::abs(gapbwd-LRG) < std::numeric_limits<double>::epsilon() || std::abs(gapfwd-LRG) < std::numeric_limits<double>::epsilon() ) vetoEvent;

      // Fill DD plots
      if (LRG > 3.) {
        if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3)) {
          _h_xsec_dd->fill( 900/GeV, weight);
        } else if (fuzzyEquals(sqrtS()/GeV, 2760, 1E-3)) {
          _h_xsec_dd->fill( 2760/GeV, weight);
        } else if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) {
          _h_xsec_dd->fill( 7000/GeV, weight);
        }
      }
    }

    void finalize() {

      // get the ratio plots: SD/inel, DD/inel
      divide(_h_xsec_sd , _h_xsec_inel, _h_frac_sd_inel);
      divide(_h_xsec_sd , _h_xsec_inel, _h_frac_dd_inel);

      scale(_h_xsec_sd,   crossSection()/millibarn/sumOfWeights());
      scale(_h_xsec_dd,   crossSection()/millibarn/sumOfWeights());
      scale(_h_xsec_inel, crossSection()/millibarn/sumOfWeights());

    }

  private:

    Scatter2DPtr _h_frac_sd_inel;
    Scatter2DPtr _h_frac_dd_inel;
    Histo1DPtr   _h_xsec_sd;
    Histo1DPtr   _h_xsec_dd;
    Histo1DPtr   _h_xsec_inel;

  };

  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I1181770);

}
