// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2005_S6080774.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  CDF_2005_S6080774::CDF_2005_S6080774() : Analysis("CDF_2005_S6080774") {
    setBeams(PROTON, ANTIPROTON);
    setNeedsCrossSection(true);
    
    FinalState fs;
    addProjection(fs, "FS");
    
    IdentifiedFinalState ifs(-0.9, 0.9, 13.0*GeV);
    ifs.acceptId(PHOTON);
    addProjection(ifs, "IFS");
  }


  void CDF_2005_S6080774::init() {

    for (size_t yAxisId=1; yAxisId<5; ++yAxisId) {
      _h_m_PP.push_back(bookHistogram1D(1, 1, yAxisId));
      _h_pT_PP.push_back(bookHistogram1D(2, 1, yAxisId));
      _h_dphi_PP.push_back(bookHistogram1D(3, 1, yAxisId));
    }

  }


  void CDF_2005_S6080774::analyze(const Event& event) {

    const double weight = event.weight();
    
    ParticleVector photons = applyProjection<IdentifiedFinalState>(event, "IFS").particles();
    
    if (photons.size() < 2 ||
        (photons[0].momentum().pT() < 14.0*GeV && photons[1].momentum().pT() < 14.0*GeV)) {
      vetoEvent;
    }
    
    // isolate photons with ET_sum in cone
    ParticleVector isolated_photons;
    ParticleVector fs = applyProjection<FinalState>(event, "FS").particles();
    foreach (const Particle& photon, photons) {
      FourMomentum mom_in_cone;
      double eta_P = photon.momentum().eta();
      double phi_P = photon.momentum().phi();
      foreach (const Particle& p, fs) {
        if (deltaR(eta_P, phi_P, p.momentum().eta(), p.momentum().phi()) < 0.4) {
          mom_in_cone += p.momentum();
        }
      }
      if (mom_in_cone.Et()-photon.momentum().Et() < 1.0*GeV) {
        isolated_photons.push_back(photon);
      }
    }
    
    if (isolated_photons.size()!=2) {
      vetoEvent;
    }
    
    FourMomentum mom_PP = isolated_photons[0].momentum() + isolated_photons[1].momentum();
    
    for (size_t i=0; i<4; ++i) {
      _h_m_PP[i]->fill(mom_PP.mass(), weight);
      _h_pT_PP[i]->fill(mom_PP.pT(), weight);
      _h_dphi_PP[i]->fill(mapAngle0ToPi(isolated_photons[0].momentum().phi()-
                                        isolated_photons[1].momentum().phi())/M_PI,
                          weight);
    }
    
  }


  void CDF_2005_S6080774::finalize() {

    for (size_t i=0; i<4; ++i) {
      scale(_h_m_PP[i], crossSection()/sumOfWeights());
      scale(_h_pT_PP[i], crossSection()/sumOfWeights());
      scale(_h_dphi_PP[i], crossSection()/sumOfWeights());
    }

  }


}
