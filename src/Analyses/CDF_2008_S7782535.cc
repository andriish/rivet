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


  // Constructor.
  // jet cuts: |eta| <= 0.7
  // Don't attempt to model the following cuts :
  //   Missing ET significance
  //   veto on additional vertices
  //   Zvtx < 50
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
  }


  void CDF_2008_S7782535::init() {
    _pTbins += 52, 80, 104, 142, 300;
    // Book histograms
    for (int i = 0; i < _NpTbins; ++i) {
      stringstream title;
      title << "Integral jet shape $\\Psi$ for $" << _pTbins[i] << " < p_\\perp < " << _pTbins[i+1] << "$"; 
      _Psi_pT[i] = bookProfile1D(i+1, 2, 1, title.str(),"r/R","$\\Psi$(r/R)");
    }
    _OneMinusPsi_vs_pT = bookDataPointSet(5, 1, 1, "$1 - \\Psi$ vs jet $p_\\perp$","$p_\\perp$ [GeV]","1-$\\Psi$(0.2/R)");
  }  


  
  // Do the analysis
  void CDF_2008_S7782535::analyze(const Event& event) {
    // Put all b-quarks in a vector
    ParticleVector bquarks;
    /// @todo Provide nicer looping
    for (GenEvent::particle_const_iterator p = event.genEvent().particles_begin(); 
         p != event.genEvent().particles_end(); ++p) {
      if ( fabs((*p)->pdg_id()) == BQUARK ) {
        bquarks.push_back(Particle(**p));
      }
    }
    
    if (bquarks.empty()) { 
      getLog() << Log::DEBUG << "No b-quarks, exiting" << endl;
      vetoEvent(event);
    }

    // Get jets 
    const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
    getLog() << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.size() << endl;
    
    const Jets& jets = jetpro.jetsByPt();
    getLog() << Log::DEBUG << "Jet list size = " << jets.size() << endl;

    // Determine the central jet axes
    _jetaxes.clear();
    foreach (const Jet& jt, jets) {
      // Only use central calorimeter jets
      FourMomentum pjet = jt.momentum();
      if (pjet.pT()/GeV > _pTbins[0] && fabs(pjet.rapidity()) < 0.7) {
        _jetaxes.push_back(pjet);
      }
    }
    if (_jetaxes.empty())  { 
      getLog() << Log::DEBUG << "No jet axes" << endl;
      vetoEvent(event);
    }

    // Determine jet shapes
    const JetShape& js = applyProjection<JetShape>(event, "JetShape");

    /// @todo Replace with foreach
    for (size_t jind = 0; jind < _jetaxes.size(); ++jind) {

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
        if      (_jetaxes[jind].pT() > _pTbins[0] && _jetaxes[jind].pT() <= _pTbins[1]) jet_pt_bin = 0;
        else if (_jetaxes[jind].pT() > _pTbins[1] && _jetaxes[jind].pT() <= _pTbins[2]) jet_pt_bin = 1;
        else if (_jetaxes[jind].pT() > _pTbins[2] && _jetaxes[jind].pT() <= _pTbins[3]) jet_pt_bin = 2;
        else if (_jetaxes[jind].pT() > _pTbins[3] && _jetaxes[jind].pT() <= _pTbins[4]) jet_pt_bin = 3;
        if (jet_pt_bin > -1) {
          // Fill each entry in profile
          for (size_t rbin = 0; rbin < js.numBins(); ++rbin) {
            const double rad_Psi = js.rMin() +(rbin+1.0) * js.interval();
            /// @todo What is this getIntJetShape(jind, rbin) index arg? Yuck...
            _Psi_pT[jet_pt_bin]->fill(rad_Psi/_Rjet, js.intJetShape(jind, rbin), event.weight() );
          }
        }
      }
    }
    
  }

  

  // Finalize
  void CDF_2008_S7782535::finalize() {  
    std::vector<double> x, y, ex, ey;
    for (size_t i = 0; i < _pTbins.size()-1; ++i) {
      x.push_back((_pTbins[i]+_pTbins[i+1])/2.0);
      /// @todo Need this to be set to match the ref histo: autobook the DPS?
      ex.push_back(0.0);
      // get entry for rad_Psi = 0.2 bin
      y.push_back(1.0 - _Psi_pT[i]->binHeight(1));
      ey.push_back(_Psi_pT[i]->binError(1));
    }
    _OneMinusPsi_vs_pT->setCoordinate(0, x, ex);
    _OneMinusPsi_vs_pT->setCoordinate(1, y, ey); 
  }


}
