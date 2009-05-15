// -*- C++ -*-
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2005_S6217184.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// Constructor
  CDF_2005_S6217184::CDF_2005_S6217184()
  { 
    setBeams(PROTON, ANTIPROTON);

    const FinalState fs(-2.0, 2.0);
    addProjection(fs, "FS");
    addProjection(FastJets(fs, FastJets::CDFMIDPOINT, 0.7), "Jets"); 
    addProjection(TotalVisibleMomentum(fs), "CalMET");
    addProjection(PVertex(), "PV");
    // Veto (anti)neutrinos, and muons with pT above 1.0 GeV
    VetoedFinalState vfs(fs);
    vfs
      .addVetoPairId(NU_E)
      .addVetoPairId(NU_MU)
      .addVetoPairId(NU_TAU)
      .addVetoDetail(MUON, 1.0*GeV, MAXDOUBLE);
    addProjection(vfs, "VFS");
    addProjection(JetShape(vfs, _jetaxes, 0.0, 0.7, 0.1, 0.3), "JetShape");

    // Specify pT bins and initialise weight entries
    /// @todo Get these numbers from bundled data files
    _pTbins.resize(19);
    double ptbins[] = { 37.0, 45.0, 55.0, 63.0, 73.0, 84.0, 97.0, 112.0, 128.0, 148.0, 
                        166.0, 186.0, 208.0, 229.0, 250.0, 277.0, 304.0, 340.0, 380.0 };
    for (size_t i = 0; i <= 18; ++i) {
      _pTbins[i]  =  ptbins[i];
      _shapeWeights[i%18] = 0.0;
    }
  }

  
  // Book histograms
  void CDF_2005_S6217184::init() {
    // 18 = 6x3 pT bins, one histogram each
    for (size_t i = 0; i < 6; ++i) { 
      for (size_t j = 0; j < 3; ++j) {
        size_t k = i*3 + j;
        stringstream ss;
        ss << "Differential jet shape $\\rho$, $p_\\perp$ bin " << k+1;
        _profhistRho_pT[k] = 
          bookProfile1D(i+1, 1, j+1, ss.str(), "$r/R$", "$\\rho(r/R)$");
        ss.str("");
        ss << "Integral jet shape $\\psi$, $p_\\perp$ bin " << k+1;
        _profhistPsi_pT[k] = 
          bookProfile1D(6+i+1, 1, j+1, ss.str(), "$r/R$", "$\\psi(r/R)$");
      }
    }    
    _profhistPsi = 
      bookProfile1D(13, 1, 1, "$\\Psi$(0.3 over $R$)",
                    "$p_\\perp^\\text{jet}$ / GeV/$c$", "$\\psi(0.3/R)$");
  }
  

  
  // Do the analysis
  void CDF_2005_S6217184::analyze(const Event& event) {
    // Find primary vertex and veto on its separation from the nominal IP
    const PVertex& pv = applyProjection<PVertex>(event, "PV");
    if (fabs(pv.position().z())/mm > 600) {
      vetoEvent(event);
    }

    // Get jets and require at least one to pass pT and y cuts
    const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt();
    getLog() << Log::DEBUG << "Jet multiplicity before cuts = " << jets.size() << endl;
    bool jetcutpass = false;
    foreach (const Jet& jt, jets) {
      const FourMomentum pj = jt.momentum();
      if (pj.pT()/GeV > 37.0 && inRange(fabs(pj.rapidity()), 0.1, 0.7)) {
        jetcutpass = true;
      }
    }
    if (!jetcutpass) vetoEvent(event);

    // Check there's not too much missing Et
    const TotalVisibleMomentum& caloMissEt = applyProjection<TotalVisibleMomentum>(event, "CalMET");
    getLog() << Log::DEBUG << "CaloMissEt.momentum().pT() = " << caloMissEt.momentum().pT() << endl;
    if ((caloMissEt.momentum().pT()/GeV) / sqrt(caloMissEt.scalarET()/GeV) > 3.5) {
      vetoEvent(event);
    }
    
    // Determine the central jet axes
    _jetaxes.clear();
    foreach (const Jet& jt, jets) {
      const FourMomentum pj = jt.momentum();
      if (fabs(pj.rapidity()) < 1.1) {
        _jetaxes.push_back(pj);
      }
    }
    
    // Calculate and histogram jet shapes
    if (!_jetaxes.empty()) {
      const double weight = event.weight();
      const JetShape& js = applyProjection<JetShape>(event, "JetShape");
      const double R_JET = 0.7;

      for (size_t jind = 0; jind < _jetaxes.size(); ++jind) {
        for (size_t ipT = 0; ipT < 18; ++ipT) {
          if (_jetaxes[jind].pT() > _pTbins[ipT] && _jetaxes[jind].pT() <= _pTbins[ipT+1]) {
            /// @todo Is this not being used?
            _shapeWeights[ipT] += weight;
            for (size_t rbin = 0; rbin < js.numBins(); ++rbin) {
              const double rad_Rho = js.rMin() + (rbin+0.5)*js.interval();
              _profhistRho_pT[ipT]->fill(rad_Rho/R_JET, js.diffJetShape(jind, rbin), weight);
              const double rad_Psi = js.rMin() +(rbin+1.0)*js.interval();
              _profhistPsi_pT[ipT]->fill(rad_Psi/R_JET, js.intJetShape(jind, rbin), weight);
            }
            _profhistPsi->fill((_pTbins[ipT] + _pTbins[ipT+1])/2.0, js.psi(jind), weight);
          }
        }
      }

    }

  }
  

  // Finalize
  void CDF_2005_S6217184::finalize() {  
    /// @todo Do the shape weighting?
  }
  
  
}
