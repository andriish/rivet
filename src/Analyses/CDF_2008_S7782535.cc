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
    addProjection(JetShape(vfs, _jetaxes, 0.0, 0.7, 0.1, 0.3, ENERGY), "JetShape");
  }


  void CDF_2008_S7782535::init() {
    _pTbins.push_back(52.);
    _pTbins.push_back(80.);
    _pTbins.push_back(104.);
    _pTbins.push_back(142.);
    _pTbins.push_back(300.);
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
    
    /// @todo Don't expose FastJet objects in Rivet analyses
    const PseudoJets& jets = jetpro.pseudoJetsByPt();
    getLog() << Log::DEBUG << "jetlist size = " << jets.size() << endl;
    // Determine the central jet axes
    FourMomentum jetaxis;
    _jetaxes.clear();
    for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      // Only Central Calorimeter jets
      /// @todo Declare this cut
      //      cout<<jt->perp() << " " << jt->rapidity() << endl;
      if (jt->perp() > _pTbins[0] && fabs(jt->rapidity()) <= 0.7) {
        jetaxis.px(jt->px());
        jetaxis.py(jt->py());
        jetaxis.pz(jt->pz());
        jetaxis.E(jt->E());
        _jetaxes.push_back(jetaxis);
      }
    }
    // Determine jet shapes
    if (_jetaxes.empty())  { 
      getLog() << Log::DEBUG << "No jet axes" << endl;
      vetoEvent(event);
    }
      
    const JetShape& jetShape = applyProjection<JetShape>(event, "JetShape");
    for (size_t jind = 0; jind < _jetaxes.size(); ++jind) { 

      /// @todo Replace this with Jet::containsParticleId
      bool bjet = false;
      foreach (const Particle& bquark,  bquarks) {
        // double dr = deltaR(_jetaxes[jind].rapidity(), _jetaxes[jind].azimuthalAngle(),
        //                    bquark->momentum().rapidity(), bquark->momentum().azimuthalAngle())
        //if (dr <= _Rjet ) {
        if (deltaR(_jetaxes[jind], bquark.momentum()) <= _Rjet ) {
          bjet = true;
          break;
        }
      }

      if(bjet) {	
        // Put jet in correct pT bin
        int jet_pt_bin = -1;
        if      (_jetaxes[jind].pT() > _pTbins[0] && _jetaxes[jind].pT() <= _pTbins[1]) jet_pt_bin = 0;
        else if (_jetaxes[jind].pT() > _pTbins[1] && _jetaxes[jind].pT() <= _pTbins[2]) jet_pt_bin = 1;
        else if (_jetaxes[jind].pT() > _pTbins[2] && _jetaxes[jind].pT() <= _pTbins[3]) jet_pt_bin = 2;
        else if (_jetaxes[jind].pT() > _pTbins[3] && _jetaxes[jind].pT() <= _pTbins[4]) jet_pt_bin = 3;
        if (jet_pt_bin > -1) {
          // Fill each entry in profile
          for (size_t rbin = 0; rbin < jetShape.getNbins(); ++rbin) {
            const double rad_Psi = jetShape.getRmin() +(rbin+1.0)*jetShape.getInterval();
            _Psi_pT[jet_pt_bin]->fill(rad_Psi/_Rjet, jetShape.getIntJetShape(jind, rbin), event.weight() );
          }
        } // end valid jet_pt_bin
      } // end bjet
    } // end loop round jets
    
  }

  

  // Finalize
  void CDF_2008_S7782535::finalize() {  
    std::vector<double> x, y, ex, ey;
    for (unsigned int i = 0; i < _pTbins.size()-1; i++) {
      x.push_back((_pTbins[i]+_pTbins[i+1])/2.);
      ex.push_back(0.);
      // get entry for  rad_Psi = 0.2 bin
      y.push_back(1.0 - _Psi_pT[i]->binHeight(1));
      ey.push_back(_Psi_pT[i]->binError(1)); 
      
    }
    _OneMinusPsi_vs_pT->setCoordinate(0,x,ex);
    _OneMinusPsi_vs_pT->setCoordinate(1,y,ey); 
  }


}
