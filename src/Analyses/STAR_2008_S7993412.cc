// -*- C++ -*-
#include "Rivet/Analyses/STAR_2008_S7993412.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  STAR_2008_S7993412::STAR_2008_S7993412()
    : Analysis("STAR_2008_S7993412")
  {
    setBeams(PROTON, PROTON);

    //full final state
    ChargedFinalState fs(-1.0, 1.0, 1.0*GeV);
    addProjection(fs, "FS");
  } 


  // Book histograms
  void STAR_2008_S7993412::init() {
    _h_Y_jet_trigger = bookProfile1D(1, 1, 1);
    _h_Y_jet_associated = bookProfile1D(2, 1, 1);
  }



  // Do the analysis 
  void STAR_2008_S7993412::analyze(const Event & event) {
    double weight = event.weight();

    // Skip if the event is empty
    const FinalState& fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty()) {
      getLog() << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
               << " because no final state found " << endl;
      vetoEvent;
    }

    foreach (const Particle& tp, fs.particles()) {
      const double triggerpT = tp.momentum().pT();
      if (triggerpT >= 2.0 && triggerpT < 5.0) {
        int N_associated = 0;
        foreach (const Particle& ap, fs.particles()) {
          if (ap.momentum().pT() > 1.5 &&
              ap.momentum().pT() < triggerpT &&
              deltaPhi(tp.momentum().phi(), ap.momentum().phi()) < 1 &&
              fabs(tp.momentum().pseudorapidity() - ap.momentum().pseudorapidity()) < 1.75) {
            N_associated += 1;
          }
        }
        //const double dPhidEta = 2 * 2*1.75;
        //_h_Y_jet_trigger->fill(triggerpT, N_associated/dPhidEta, weight);
        _h_Y_jet_trigger->fill(triggerpT, N_associated, weight);
      }
    }
  }



  // Finalize
  void STAR_2008_S7993412::finalize() {
    /// @todo Use the generator cross-section
    //_h_total_cross_section->fill(crossSection());
    //normalize(_h_jet_pT_MB, 16603100);
    //normalize(_h_jet_pT_HT, 1808234);
  }

}
