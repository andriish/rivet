// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/D0_1996_S3324664.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  D0_1996_S3324664::D0_1996_S3324664() : Analysis("D0_1996_S3324664") {
    setBeams(PROTON, ANTIPROTON);
    setNeedsCrossSection(false);

    const FinalState fs(-3.2, 3.2);
    addProjection(fs, "FS");
    /// @todo Use correct jet algorithm
    addProjection(FastJets(fs, FastJets::D0ILCONE, 0.7, 20.0*GeV), "ConeJets");
  }


  void D0_1996_S3324664::init() {

    _h_deta = bookHistogram1D(1, 1, 1);
    _h_dphi.addHistogram(0.0, 2.0, bookHistogram1D(2, 1, 1));
    _h_dphi.addHistogram(2.0, 4.0, bookHistogram1D(2, 1, 2));
    _h_dphi.addHistogram(4.0, 6.0, bookHistogram1D(2, 1, 3));
    _h_cosdphi_deta = bookProfile1D(3, 1, 1);
  }


  void D0_1996_S3324664::analyze(const Event& event) {

    const double weight = event.weight();
    
    Jets jets;
    foreach (const Jet& jet, applyProjection<FastJets>(event, "ConeJets").jets()) {
      if (fabs(jet.momentum().eta())<3.0) {
        jets.push_back(jet);
      }
    }
    
    if (jets.size()<2) {
      vetoEvent;
    }
    
    FourMomentum minjet=jets[0].momentum();
    FourMomentum maxjet=jets[1].momentum();
    double mineta = minjet.eta();
    double maxeta = maxjet.eta();
    
    foreach(const Jet& jet, jets) {
      double eta = jet.momentum().eta();
      if (eta < mineta) {
        minjet = jet.momentum();
        mineta = eta;
      }
      else if (eta > maxeta) {
        maxjet = jet.momentum();
        maxeta = eta;
      }
    }
    
    if (minjet.Et()<50*GeV && maxjet.Et()<50.0*GeV) {
      vetoEvent;
    }
    
    double deta = maxjet.eta()-minjet.eta();
    double dphi = mapAngle0To2Pi(maxjet.phi()-minjet.phi());
    
    _h_deta->fill(deta, weight);
    _h_dphi.fill(deta, 1.0-dphi/M_PI, weight);
    _h_cosdphi_deta->fill(deta, cos(M_PI-dphi), weight);

  }


  void D0_1996_S3324664::finalize() {

    normalize(_h_deta, 8830.0); // not normalised to cross section, but to #events
    
    // I have no idea what this is normalised to, in the paper it says unity!?!
    foreach (IHistogram1D* histo, _h_dphi.getHistograms()) {
      normalize(histo, 0.0798);
    }

  }


}
