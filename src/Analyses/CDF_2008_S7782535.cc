// -*- C++ -*-
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2008_S7782535.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  CDF_2008_S7782535::CDF_2008_S7782535()
    : _Rjet(0.7) , _NpTbins(4)
  {
    setBeams(PROTON, ANTIPROTON);
    
    const FinalState fs(-3.6, 3.6);
    addProjection(fs, "FS");
    // Veto (anti)neutrinos, and muons with pT above 1.0 GeV
    VetoedFinalState vfs(fs);
    vfs
      .addVetoPairId(NU_E)
      .addVetoPairId(NU_MU)
      .addVetoPairId(NU_TAU)
      .addVetoDetail(MUON, 1.0*GeV, MAXDOUBLE);
    addProjection(vfs, "VFS");
    addProjection(FastJets(vfs, FastJets::CDFMIDPOINT, 0.7), "Jets");
    addProjection(JetShape(vfs, _jetaxes, 0.0, 0.7, 0.1, 0.3), "JetShape");

    // We don't attempt to model the following cuts:
    //  * missing ET significance
    //  * veto on additional vertices
    //  * Z_vtx < 50 cm
  }


  void CDF_2008_S7782535::init() {
    _pTbins += 52, 80, 104, 142, 300;
    // Book histograms
    for (int i = 0; i < _NpTbins; ++i) {
      stringstream title;
      title << "Integral jet shape $\\Psi$ for $" << _pTbins[i] 
            << " < p_\\perp < " << _pTbins[i+1] << "$"; 
      _h_Psi_pT[i] = bookProfile1D(i+1, 2, 1, title.str(), 
                                   "r/R", "$\\Psi$(r/R)");
    }
    _h_OneMinusPsi_vs_pT = 
      bookDataPointSet(5, 1, 1, "$1 - \\Psi$ vs jet $p_\\perp$", 
                       "$p_\\perp$ [GeV]", "1-$\\Psi$(0.2/R)");
  }  

  
  // Do the analysis
  void CDF_2008_S7782535::analyze(const Event& event) {
    // Put all b-quarks in a vector
    ParticleVector bquarks;
    /// @todo Provide nicer looping over HepMC event contents
    for (GenEvent::particle_const_iterator p = event.genEvent().particles_begin(); 
         p != event.genEvent().particles_end(); ++p) {
      if (abs((*p)->pdg_id()) == BQUARK) {
        bquarks.push_back(Particle(**p));
      }
    }
    if (bquarks.empty()) { 
      getLog() << Log::DEBUG << "No b-quarks, exiting" << endl;
      vetoEvent;
    }

    // Get jets     
    const Jets& jets = applyProjection<FastJets>(event, "Jets").jetsByPt();
    getLog() << Log::DEBUG << "Jet multiplicity before any pT cut = " << jets.size() << endl;

    // Determine the central jet axes
    _jetaxes.clear();
    foreach (const Jet& jt, jets) {
      // Only use central calorimeter jets
      FourMomentum pjet = jt.momentum();
      /// @todo Really y rather than eta?
      if (pjet.pT()/GeV > _pTbins[0] && fabs(pjet.rapidity()) < 0.7) {
        /// @todo Only push_back jets which contain particles with a b quark parent?
        _jetaxes.push_back(pjet);
      }
    }
    if (_jetaxes.empty())  {
      getLog() << Log::DEBUG << "No jet axes in acceptance" << endl;
      vetoEvent;
    }

    // Determine jet shapes
    const JetShape& js = applyProjection<JetShape>(event, "JetShape");

    /// @todo Replace with foreach
    for (size_t jind = 0; jind < _jetaxes.size(); ++jind) {

      /// @todo Replace with only running over jets containing b quark descendants
      bool bjet = false;
      foreach (const Particle& bquark,  bquarks) {
        if (deltaR(_jetaxes[jind], bquark.momentum()) < _Rjet) {
          bjet = true;
          break;
        }
      }

      if (bjet) {	
        // Put jet in correct pT bin
        int jet_pt_bin = -1;
        for (size_t i = 0; i < 4; ++i) {
          if (_jetaxes[jind].pT() > _pTbins[i] && 
              _jetaxes[jind].pT() <= _pTbins[i+1]) jet_pt_bin = i;
        }
        if (jet_pt_bin > -1) {
          // Fill each entry in profile
          for (size_t rbin = 0; rbin < js.numBins(); ++rbin) {
            const double rad_Psi = js.rMin() +(rbin+1.0) * js.interval();
            /// @todo Yuck... JetShape's interface sucks
            _h_Psi_pT[jet_pt_bin]->fill(rad_Psi/_Rjet, js.intJetShape(jind, rbin), event.weight() );
          }
        }
      }
    }
    
  }

  

  // Finalize
  void CDF_2008_S7782535::finalize() {  
    std::vector<double> y, ey;
    for (size_t i = 0; i < _pTbins.size()-1; ++i) {
      // Get entry for rad_Psi = 0.2 bin
      AIDA::IProfile1D* ph_i = _h_Psi_pT[i];
      y.push_back(1.0 - ph_i->binHeight(1));
      ey.push_back(ph_i->binError(1));
    }
    _h_OneMinusPsi_vs_pT->setCoordinate(1, y, ey); 
  }


}
