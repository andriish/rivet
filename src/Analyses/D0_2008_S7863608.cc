// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S7863608.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2008_S7863608::D0_2008_S7863608()
  {
    setBeams(PROTON, ANTIPROTON);
    
    /// @todo Use cross-section from generator
    //setNeedsCrossSection(true);

    ZFinder zfinder(-1.7, 1.7, 15.0*GeV, MUON, 65.0*GeV, 115.0*GeV, 0.2);
    addProjection(zfinder, "ZFinder");
    
    FastJets conefinder(zfinder.remainingFinalState(), FastJets::D0ILCONE, 0.5, 20.0*GeV);
    addProjection(conefinder, "ConeFinder");
  }



  // Book histograms
  void D0_2008_S7863608::init() {

    /// @todo Dividing through by measured Z cross-section would be nice...
    _h_jet_pT_cross_section = bookHistogram1D(1, 1, 1, "Differential cross section in leading jet $p_\\perp$",
                                              "$p_{\\perp}$(1st jet) [GeV]", "$\\text{d}\\sigma/\\text{d}p_{\\perp}$(1st jet)");
    _h_jet_y_cross_section = bookHistogram1D(2, 1, 1, "Differential cross section in leading jet rapidity",
                                             "$|y|$(1st jet)", "$\\text{d}\\sigma/\\text{d}|y|$(1st jet)");
    _h_Z_pT_cross_section = bookHistogram1D(3, 1, 1, "Differential cross section in Z/$\\gamma*$ $p_\\perp$",
                                            "$p_{\\perp}$(Z) [GeV]", "$\\text{d}\\sigma/\\text{d}p_{\\perp}$(Z)");
    _h_Z_y_cross_section = bookHistogram1D(4, 1, 1, "Differential cross section in Z/$\\gamma*$ rapidity",
                                           "$|y|$(Z)", "$\\text{d}\\sigma/\\text{d}|y|$(Z)");
    _h_total_cross_section = bookHistogram1D(5, 1, 1, "Total Z + jet cross section", "$\\sqrt{s}$", "$\\sigma$");
    
  }



  // Do the analysis 
  void D0_2008_S7863608::analyze(const Event & e) {
    double weight = e.weight();

    const ZFinder& zfinder = applyProjection<ZFinder>(e, "ZFinder");
    if (zfinder.particles().size()==1) {
      const JetAlg& jetpro = applyProjection<JetAlg>(e, "ConeFinder");
      const Jets& jets = jetpro.jetsByPt(20.0*GeV);
      Jets jets_cut;
      foreach (const Jet& j, jets) {
        if (fabs(j.momentum().pseudorapidity()) < 2.8) {
          jets_cut.push_back(j);
        }
      }
      
      // Return if there are no jets:
      if(jets_cut.size()<1) {
        getLog() << Log::DEBUG << "Skipping event " << e.genEvent().event_number()
            << " because no jets pass cuts " << endl;
        vetoEvent;
      }
      
      // cut on Delta R between jet and muons
      foreach (const Jet& j, jets_cut) {
        foreach (const Particle& mu, zfinder.constituentsFinalState().particles()) {
          if (deltaR(mu.momentum().pseudorapidity(), mu.momentum().azimuthalAngle(),
                     j.momentum().pseudorapidity(), j.momentum().azimuthalAngle()) < 0.5) {
            vetoEvent;
          }
        }
      }
      
      const FourMomentum Zmom = zfinder.particles()[0].momentum();

      // In jet pT
      _h_jet_pT_cross_section->fill( jets_cut[0].momentum().pT(), weight);
      _h_jet_y_cross_section->fill( fabs(jets_cut[0].momentum().rapidity()), weight);
  
      // In Z pT
      _h_Z_pT_cross_section->fill(Zmom.pT(), weight);
      _h_Z_y_cross_section->fill(fabs(Zmom.rapidity()), weight);
    }
  }



  // Finalize
  void D0_2008_S7863608::finalize() {
    /// @todo Use the generator cross-section
    //_h_total_cross_section->fill(crossSection());
    normalize(_h_jet_pT_cross_section, 18.7);
    normalize(_h_jet_y_cross_section, 18.7);
    normalize(_h_Z_pT_cross_section, 18.7);
    normalize(_h_Z_y_cross_section, 18.7);
  }

}
