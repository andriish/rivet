// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S7662670.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2008_S7662670::D0_2008_S7662670()
  {
    setBeams(PROTON, ANTIPROTON);
    setNeedsCrossSection(true);

    //full final state
    FinalState fs(-5.0, 5.0);
    addProjection(fs, "FS");

    FastJets jetpro(fs, FastJets::D0ILCONE, 0.7, 6*GeV);
    addProjection(jetpro, "Jets");
  }



  // Book histograms
  void D0_2008_S7662670::init() 
  {
    const string basetitle = "Inclusive jet $p_\\perp$, ";
    const string xlabel = "Jet $p_\\perp$ / GeV";
    const string ylabel = "$\\mathrm{d}{\\sigma}/\\mathrm{d}{p_\\perp}$";

    _h_dsigdptdy_y00_04 = 
      bookHistogram1D(1, 1, 1, basetitle + "$0.0 < |y| < 0.4$", xlabel, ylabel);
    _h_dsigdptdy_y04_08 =
      bookHistogram1D(2, 1, 1, basetitle + "$0.4 < |y| < 0.8$", xlabel, ylabel);
    _h_dsigdptdy_y08_12 =
      bookHistogram1D(3, 1, 1, basetitle + "$0.8 < |y| < 1.2$", xlabel, ylabel);
    _h_dsigdptdy_y12_16 =
      bookHistogram1D(4, 1, 1, basetitle + "$1.2 < |y| < 1.6$", xlabel, ylabel);
    _h_dsigdptdy_y16_20 =
      bookHistogram1D(5, 1, 1, basetitle + "$1.6 < |y| < 2.0$", xlabel, ylabel);
    _h_dsigdptdy_y20_24 =
      bookHistogram1D(6, 1, 1, basetitle + "$2.0 < |y| < 2.4$", xlabel, ylabel);
  }



  // Do the analysis 
  void D0_2008_S7662670::analyze(const Event & event) {
    const double weight = event.weight();

    // Skip if the event is empty
    const FinalState& fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty()) {
      getLog() << Log::DEBUG << "Empty event!" << endl;
      vetoEvent;
    }

    // Find the jets
    const JetAlg& jetpro = applyProjection<JetAlg>(event, "Jets");
    // If there are no jets, skip the event
    if (jetpro.jets().size() == 0) {
      getLog() << Log::DEBUG << "No jets found" << endl;
      vetoEvent;
    }

    // Fill histo for each jet
    foreach (const Jet& j, jetpro.jets()) {
      const double pt = j.momentum().pT();
      const double y = fabs(j.momentum().rapidity());
      if (pt/GeV > 50) {
        getLog() << Log::TRACE << "Filling histos: pT = " << pt/GeV 
                 << ", |y| = " << y << endl;
        if (y < 0.4) {
          _h_dsigdptdy_y00_04->fill(pt/GeV, weight);
        } else if (y < 0.8) {
          _h_dsigdptdy_y04_08->fill(pt/GeV, weight);
        } else if (y < 1.2) {
          _h_dsigdptdy_y08_12->fill(pt/GeV, weight);
        } else if (y < 1.6) {
          _h_dsigdptdy_y12_16->fill(pt/GeV, weight);
        } else if (y < 2.0) {
          _h_dsigdptdy_y16_20->fill(pt/GeV, weight);
        } else if (y < 2.4) {
          _h_dsigdptdy_y20_24->fill(pt/GeV, weight);
        }
      }
    }

  }


  // Finalize
  void D0_2008_S7662670::finalize() {
    /// Scale by L_eff = sig_MC * L_exp / num_MC
    const double lumi_mc = sumOfWeights() / crossSection();
    const double scalefactor =  1 / lumi_mc;
    scale(_h_dsigdptdy_y00_04, scalefactor);
    scale(_h_dsigdptdy_y04_08, scalefactor);
    scale(_h_dsigdptdy_y08_12, scalefactor);
    scale(_h_dsigdptdy_y12_16, scalefactor);
    scale(_h_dsigdptdy_y16_20, scalefactor);
    scale(_h_dsigdptdy_y20_24, scalefactor);
  }

}
