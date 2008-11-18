// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S7863608.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2008_S7863608::D0_2008_S7863608() : _events(0) 
  {
    setBeams(PROTON, ANTIPROTON);
    //full final state
    FinalState fs(-5.0, 5.0);
    addProjection(fs, "FS");
 
    //leading muons in tracking acceptance
    LeadingParticlesFinalState lpfs(fs, -1.7, 1.7);
    lpfs.addParticleId(13).addParticleId(-13);
    addProjection(lpfs, "LeadingMuons");

    //Vetoed fs for jets
    VetoedFinalState vfs(fs);
    // veto the muons from Z decay  
    vfs.addVetoOnThisFinalState(lpfs);

    //also need to exclude photons in dR<0.2 around each muon

    //jets
    D0ILConeJets jets(vfs);
    addProjection(jets, "Jets");

  } 



  // Book histograms
  void D0_2008_S7863608::init() {

    /// @todo Dividing through by measured Z cross-section would be nice...
    _h_jet_pT_cross_section = bookHistogram1D(1, 1, 1, "Differential cross section in leading jet pT");
    _h_jet_y_cross_section = bookHistogram1D(2, 1, 1, "Differential cross section in leading jet rapidity");
    _h_Z_pT_cross_section = bookHistogram1D(3, 1, 1, "Differential cross section in leading Z/gamma* pT");
    _h_Z_y_cross_section = bookHistogram1D(4, 1, 1, "Differential cross section in leading Z/gamma* rapidity");
    _h_total_cross_section = bookHistogram1D(5, 1, 1, "Total Z + jet cross section");
    
  }



  // Do the analysis 
  void D0_2008_S7863608::analyze(const Event & event) {

    Log & log = getLog();
    log << Log::DEBUG << "Starting analyzing" << endl;

    double weight = event.weight();
    log << Log::DEBUG << "Event weight is " << weight << endl;

    //skip if the event is empty
    const FinalState & fs = applyProjection<FinalState>(event, "FS");
    if (fs.isEmpty()) {
      log << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
        << " because no final state pair found " << endl;
      return;
    }


    //find the Z candidates
    const FinalState & muonfs = applyProjection<FinalState>(event, "LeadingMuons");
    //if there are no muons in the FinalState skip event,
    if (muonfs.particles().size()!=2) {
      log << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
        << " because no muon pair found " << endl;
      return;
    }


    
    //now build the jets on a fs without the muons from Z
    const D0ILConeJets & jetpro = applyProjection<D0ILConeJets>(event, "Jets");
    
    /// @todo Can we make this neater, using the JetAlg interface and the built-in sorting?

    //take the leading jet with pt > 20, |y|<2.8:
    const Jets &jets = jetpro.jets();
    Jets jets_cut;
    foreach (const Jet& j, jets) {
      if (j.momentum().pT()/GeV > 20 && fabs(j.momentum().rapidity()) < 2.8)
	jets_cut.push_back(j);
    }
    log << Log::DEBUG << "jets above 20 GeV size = " << jets_cut.size() << endl;
    //return if there are no jets:
    if(jets_cut.size()<1) {
      log << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
	  << " because no jets pass cuts " << endl;
      return;
    }

    //sort by pT:
    sort(jets_cut.begin(), jets_cut.end(), cmpJetsByPt);




    //calculate the Z pT, rapidity:
    const ParticleVector muons = muonfs.particles();
    const FourMomentum Zmom = muons[0].momentum() + muons[1].momentum();

    if(Zmom.mass()/GeV<65 || Zmom.mass()/GeV>115) {
       log << Log::DEBUG << "Skipping event " << event.genEvent().event_number()
	  << " because failing Z mass window " << endl;
      return;
    }



    
    //now plot:
    _events += weight;
    
    _h_jet_pT_cross_section->fill( jets_cut[0].momentum().pT(), weight);
    _h_jet_y_cross_section->fill( fabs(jets_cut[0].momentum().rapidity()), weight);
    
    
    _h_Z_pT_cross_section->fill(Zmom.pT(), weight);
    _h_Z_y_cross_section->fill(fabs(Zmom.rapidity()), weight);

    //    _h_total_cross_section = bookHistogram1D     _crossSectionRatio->fill(1, weight);

  }







  // Finalize
  void D0_2008_S7863608::finalize() {
    normalize(_h_jet_pT_cross_section, 18.7);
    normalize(_h_jet_y_cross_section, 18.7);
    normalize(_h_Z_pT_cross_section, 18.7);
    normalize(_h_Z_y_cross_section, 18.7);
  }

}
