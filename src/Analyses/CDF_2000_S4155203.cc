// -*- C++ -*-
// CDF Z pT analysis

#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2000_S4155203.hh"
#include "Rivet/RivetAIDA.hh"

namespace Rivet {


  // Book histograms
  void CDF_2000_S4155203::init() {
    _hist_zpt = bookHistogram1D(1, 1, 1, "$Z p_\\perp$");
  }


  // Do the analysis
  void CDF_2000_S4155203::analyze(const Event& e) {
    const double weight = e.weight();

    // Get the leptons
    const ParticleVector& leptons = applyProjection<ChargedLeptons>(e, "CL").chargedLeptons();
    ParticleVector electrons;
    electrons.clear();

    bool has_central_electron = false;
    for (size_t i = 0 ; i < leptons.size() ; ++i) {
      // Only take electrons
      if (abs(leptons[i].pdgId()) != ELECTRON) continue;
      
      // Check for electrons in the fiducial regions
      if (fabs(leptons[i].momentum().pseudorapidity()) < 1.1 &&
          leptons[i].momentum().pT()/GeV > 20) {
        has_central_electron = true;
        electrons.push_back(leptons[i]);
      }
      else if (fabs(leptons[i].momentum().pseudorapidity()) > 1.1 &&
               fabs(leptons[i].momentum().pseudorapidity()) < 4.2 &&
               leptons[i].momentum().pT()/GeV > 15) {
        electrons.push_back(leptons[i]);
      }
    }

    // One electron needs to be in the central region
    if (! has_central_electron) {
      vetoEvent(e);
    }

    // We want exactly two electrons with opposite charge
    getLog() << Log::DEBUG << "electron multiplicity = " << electrons.size() << endl;
    if (electrons.size() != 2 || electrons[0].pdgId() != -electrons[1].pdgId() ) {
      vetoEvent(e);
    }

    FourVector dilepton = electrons[0].momentum() + electrons[1].momentum();
    const double mZ = mass(dilepton);
    const double pTZ = pT(dilepton);

    // Lepton pair should have an invariant mass between 66 and 116
    if (mZ/GeV < 66.0 || mZ/GeV > 116.0) {
      vetoEvent(e);
    }

    getLog() << Log::DEBUG << "dilepton mass = " << mZ/GeV  << endl;
    getLog() << Log::DEBUG << "dilepton pT   = " << pTZ/GeV << endl;

    // Fill the histograms
    _hist_zpt->fill(pTZ/GeV, weight);
  }


  void CDF_2000_S4155203::finalize() {
    // Normalize to the experimental cross-section
    normalize(_hist_zpt, 247.4);
  }


}
