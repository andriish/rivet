// -*- C++ -*-
#include "Rivet/Analyses/D0_2009_S8202443.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2009_S8202443::D0_2009_S8202443()
  {
    setBeams(PROTON, ANTIPROTON);
    
    /// @todo Use cross-section from generator
    //setNeedsCrossSection(true);

    //leptons in tracking acceptance
    std::vector<std::pair<double, double> > etaRanges;
    etaRanges.push_back(make_pair(-2.5, -1.5));
    etaRanges.push_back(make_pair(-1.1, 1.1));
    etaRanges.push_back(make_pair(1.5, 2.5));
    ZFinder zfinder(etaRanges, 25.0*GeV, ELECTRON,
                    65.0*GeV, 115.0*GeV, 0.2);
    addProjection(zfinder, "ZFinder");
    
    D0ILConeJets conefinder(zfinder.remainingFinalState(), 0.5, 20.0*GeV);
    addProjection(conefinder, "ConeFinder");
  } 



  // Book histograms
  void D0_2009_S8202443::init() {
    /// @todo
  }



  // Do the analysis 
  void D0_2009_S8202443::analyze(const Event & e) {
    double weight = e.weight();

    const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
    if (zfinder.particles().size()!=1) {
      getLog() << Log::DEBUG << "Skipping event " << e.genEvent().event_number()
               << " because no unique lepton pair found." << endl;
      vetoEvent(e);
    }
    const D0ILConeJets& jetpro = applyProjection<D0ILConeJets>(e, "ConeFinder");
    const Jets& jets = jetpro.jets();
    Jets jets_cut;
    foreach (const Jet& j, jets) {
      if (fabs(j.momentum().pseudorapidity()) < 2.5) {
        jets_cut.push_back(j);
      }
    }

    // Sort by pT:
    sort(jets_cut.begin(), jets_cut.end(), cmpJetsByPt);

//     if (jets_cut.size()>0) {
//       _h_jet1_pT->fill(jets_cut[0].momentum().pT()/GeV, weight);
//     }
//     if (jets_cut.size()>1) {
//       _h_jet2_pT->fill(jets_cut[1].momentum().pT()/GeV, weight);
//     }
//     if (jets_cut.size()>2) {
//       _h_jet3_pT->fill(jets_cut[2].momentum().pT()/GeV, weight);
//     }
  }



  // Finalize
  void D0_2009_S8202443::finalize() {
    /// @todo Use the generator cross-section
  }

}
