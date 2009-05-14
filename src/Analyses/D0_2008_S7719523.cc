// -*- C++ -*-
#include "Rivet/Analyses/D0_2008_S7719523.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  D0_2008_S7719523::D0_2008_S7719523()
  {
    setBeams(PROTON, ANTIPROTON);
    
    /// @todo Use cross-section from generator
    //setNeedsCrossSection(true);

    // General FS
    FinalState fs(-5.0, 5.0);
    addProjection(fs, "FS");
 
    // Get leading photon
    LeadingParticlesFinalState photonfs(fs, -1.0, 1.0);
    photonfs.addParticleId(PHOTON);
    addProjection(photonfs, "LeadingPhoton");

    // FS for jets excludes the leading photon
    VetoedFinalState vfs(fs);
    vfs.addVetoOnThisFinalState(photonfs);
    addProjection(vfs, "JetFS");
  } 



  // Book histograms
  void D0_2008_S7719523::init() {
    const string xlabel = "$p_\\perp(\\gamma_\\text{lead})$ / GeV";
    /// @todo Cross-section units in label
    const string ylabel = "$\\mathrm{d}{\\sigma}/\\mathrm{d}{p_\\perp(\\gamma_\\text{lead})}$";
    const string basetitle = "Leading photon $p_\\perp$ ";

    _h_central_same_cross_section = 
      bookHistogram1D("d01-x01-y01", basetitle + "(central jets, same-sign rapidity)", xlabel, ylabel);
    _h_central_opp_cross_section  = 
      bookHistogram1D("d02-x01-y01", basetitle + "(central jets, opp-sign rapidity)", xlabel, ylabel);
    _h_forward_same_cross_section = 
      bookHistogram1D("d03-x01-y01", basetitle + "(forward jets, same-sign rapidity)", xlabel, ylabel);
    _h_forward_opp_cross_section  = 
      bookHistogram1D("d04-x01-y01", basetitle + "(forward jets, opp-sign rapidity)", xlabel, ylabel); 
  }



  // Do the analysis 
  void D0_2008_S7719523::analyze(const Event& event) {
    const double weight = event.weight();

    // Get the photon
    const FinalState& photonfs = applyProjection<FinalState>(event, "LeadingPhoton");
    if (photonfs.particles().size() != 1) {
      getLog() << Log::DEBUG << "No photon found" << endl;
      vetoEvent(event);
    }
    const FourMomentum photon = photonfs.particles().front().momentum();
    if (photon.pT()/GeV < 30) {
      getLog() << Log::DEBUG << "Leading photon has pT < 30 GeV: " << photon.pT()/GeV << endl;
      vetoEvent(event);
    }

    // Get all charged particles
    const FinalState& fs = applyProjection<FinalState>(event, "JetFS");
    if (fs.isEmpty()) {
      vetoEvent(event);
    }

    // Isolate photon by ensuring that a 0.4 cone around it contains less than 7% of the photon's energy
    const double egamma = photon.E();
    double econe = 0.0;
    foreach (const Particle& p, fs.particles()) {
      const double dr = deltaR(photon.pseudorapidity(), photon.azimuthalAngle(),
                               p.momentum().pseudorapidity(), p.momentum().azimuthalAngle());
      if (dr < 0.2) {
        econe += p.momentum().E();
        // Veto as soon as E_cone gets larger
        if (econe/egamma > 0.07) {
          getLog() << Log::DEBUG << "Vetoing event because photon is insufficiently isolated" << endl;
          vetoEvent(event);
        }
      }
    }


    /// @todo Allow proj creation w/o FS as ctor arg, so that calc can be used more easily.
    FastJets jetpro(fs, FastJets::D0ILCONE, 0.7); //< @todo This fs arg makes no sense!
    jetpro.calc(fs.particles());
    Jets isolated_jets;
    foreach (const Jet& j, jetpro.jets()) {
      const FourMomentum pjet = j.momentum();
      const double dr = deltaR(photon.pseudorapidity(), photon.azimuthalAngle(),
                               pjet.pseudorapidity(), pjet.azimuthalAngle());
      if (dr > 0.7 && pjet.pT()/GeV > 15) {
        isolated_jets.push_back(j);
      }
    }
    
    getLog() << Log::DEBUG << "Num jets after isolation and pT cuts = " 
             << isolated_jets.size() << endl;
    if (isolated_jets.empty()) {
      getLog() << Log::DEBUG << "No jets pass cuts" << endl;
      vetoEvent(event);
    }

    // Sort by pT and get leading jet
    sort(isolated_jets.begin(), isolated_jets.end(), cmpJetsByPt);
    const FourMomentum leadingJet = isolated_jets.front().momentum();
    int photon_jet_sign = sign( leadingJet.rapidity() * photon.rapidity() );

    // Veto if leading jet is outside plotted rapidity regions
    const double abs_y1 = fabs(leadingJet.rapidity());
    if (inRange(abs_y1, 0.8, 1.5) || abs_y1 > 2.5) {
      getLog() << Log::DEBUG << "Leading jet falls outside acceptance range; |y1| = " 
               << abs_y1 << endl;
      vetoEvent(event);
    }

    // Fill histos
    if (fabs(leadingJet.rapidity()) < 0.8) { 
      if (photon_jet_sign >= 1) {
        _h_central_same_cross_section->fill(photon.pT(), weight);
      } else {
        _h_central_opp_cross_section->fill(photon.pT(), weight);
      }
    } else if (inRange( fabs(leadingJet.rapidity()), 1.5, 2.5)) {
      if (photon_jet_sign >= 1) {
        _h_forward_same_cross_section->fill(photon.pT(), weight);
      } else {
        _h_forward_opp_cross_section->fill(photon.pT(), weight); 
      }
    }

  }



  // Finalize
  void D0_2008_S7719523::finalize() {
    // Cross-section ratios (6 plots)
    // Central/central and forward/forward ratios
    AIDA::IHistogramFactory& hf = histogramFactory();
    const string dir = histoDir();

    hf.divide(dir + "/d05-x01-y01", *_h_central_opp_cross_section, *_h_central_same_cross_section);
    hf.divide(dir + "/d08-x01-y01", *_h_forward_opp_cross_section, *_h_forward_same_cross_section);
    // Central/forward ratio combinations
    /// @todo Bins don't match
    // hf.divide(dir + "/d06-x01-y01", *_h_central_same_cross_section, *_h_forward_same_cross_section);
    // hf.divide(dir + "/d07-x01-y01", *_h_central_opp_cross_section,  *_h_forward_same_cross_section);
    // hf.divide(dir + "/d09-x01-y01", *_h_central_same_cross_section, *_h_forward_opp_cross_section);
    // hf.divide(dir + "/d10-x01-y01", *_h_central_opp_cross_section,  *_h_forward_opp_cross_section);

    /// @todo Use the generator cross-section
    // Must happen *after* the divs, since otherwise the pointers are null!
    //_h_total_cross_section->fill(crossSection());
    normalize(_h_central_same_cross_section, 347.4);
    normalize(_h_central_opp_cross_section,  281.8);
    normalize(_h_forward_same_cross_section, 164.8);
    normalize(_h_forward_opp_cross_section,   81.5);
  }


}
