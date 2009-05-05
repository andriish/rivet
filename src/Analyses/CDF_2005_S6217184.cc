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
    : _pvzmax(600*mm), _Rjet(0.7)
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
    addProjection(JetShape(vfs, _jetaxes, 0.0, 0.7, 0.1, 0.3, ENERGY), 
                  "JetShape");

    // Specify pT bins and initialise weight entries
    /// @todo Get these numbers from bundled data files
    _pTbins.resize(19);
    double ptbins[] = { 37.0, 45.0, 55.0, 63.0, 73.0, 84.0, 97.0, 112.0, 128.0, 148.0, 
                        166.0, 186.0, 208.0, 229.0, 250.0, 277.0, 304.0, 340.0, 380.0 };
    for (size_t i = 0; i <= 18; ++i) {
      _pTbins[i]  =  ptbins[i];
      _ShapeWeights[i%18] = 0.0;
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
    /// @todo Improve title... "0.3 over R" means what?)
    _profhistPsi = 
      bookProfile1D(13, 1, 1, "$\\Psi$(0.3 over $R$)",
                    "$\\psi(0.3/R)$", "p_\\perp^\\text{jet} / GeV/$c$");
  }
  

  
  // Do the analysis
  void CDF_2005_S6217184::analyze(const Event& event) {
    // Find primary vertex and veto on its separation from the nominal IP
    const PVertex& pv = applyProjection<PVertex>(event, "PV");
    if (fabs(pv.position().z())/mm > _pvzmax) vetoEvent(event);

    // Analyse and print some info  
    const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
    getLog() << Log::DEBUG << "Jet multiplicity before any pT cut = " << jetpro.size() << endl;
    
    /// @todo Don't expose FastJet objects in Rivet analyses: the FastJets projection
    /// should convert them to Rivet 4-momentum classes (or similar).
    const PseudoJets& jets = jetpro.pseudoJetsByPt();
    getLog() << Log::DEBUG << "jetlist size = " << jets.size() << endl;
    size_t Njet = 0;
    bool jetcutpass = false;
    for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      getLog() << Log::DEBUG << "List item pT = " << jt->perp() << " E=" << jt->E() << " pz=" << jt->pz() << endl;
      /// @todo Declare these cuts, and use units
      if (jt->perp() > 37.0 && fabs(jt->rapidity()) > 0.1 && fabs(jt->rapidity()) < 0.7) jetcutpass = true;
      ++Njet;
      getLog() << Log::DEBUG << "Jet pT =" << jt->perp() << " y=" << jt->rapidity() << " phi=" << jt->phi() << endl; 
    }
    if (!jetcutpass) vetoEvent(event);

    const TotalVisibleMomentum& caloMissEt = applyProjection<TotalVisibleMomentum>(event, "CalMET");
    getLog() << Log::DEBUG << "CaloMissEt.momentum().pT() = " << caloMissEt.momentum().pT() << endl;
    /// @todo Declare this cut, and use units
    if (caloMissEt.momentum().pT()/sqrt(caloMissEt.scalarET()) > 3.5) {
      vetoEvent(event);
    }
    
    // Determine the central jet axes
    FourMomentum jetaxis;
    _jetaxes.clear();
    for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
      // Only Central Calorimeter jets
      /// @todo Declare this cut
      if (fabs(jt->rapidity()) < 1.1) {
        jetaxis.px(jt->px());
        jetaxis.py(jt->py());
        jetaxis.pz(jt->pz());
        jetaxis.E(jt->E());
        _jetaxes.push_back(jetaxis);
      }
    }
    
    // Determine jet shapes
    if (!_jetaxes.empty()) { 
      const JetShape& jetShape = applyProjection<JetShape>(event, "JetShape");
      for (size_t jind = 0; jind < _jetaxes.size(); ++jind) {
        for (size_t ipT = 0; ipT < 18; ++ipT) {
          if (_jetaxes[jind].pT() > _pTbins[ipT] && _jetaxes[jind].pT() <= _pTbins[ipT+1]) {
            _ShapeWeights[ipT] += event.weight(); 
            for (size_t rbin = 0; rbin < jetShape.getNbins(); ++rbin) {
              const double rad_Rho = jetShape.getRmin() + (rbin+0.5)*jetShape.getInterval();
              _profhistRho_pT[ipT]->fill(rad_Rho/_Rjet, jetShape.getDiffJetShape(jind, rbin), event.weight() );
              const double rad_Psi = jetShape.getRmin() +(rbin+1.0)*jetShape.getInterval();
              _profhistPsi_pT[ipT]->fill(rad_Psi/_Rjet, jetShape.getIntJetShape(jind, rbin), event.weight() );
            }
            _profhistPsi->fill((_pTbins[ipT]+_pTbins[ipT+1])/2., jetShape.getPsi(jind), event.weight());
          }
        }
      }
    }
  }
  

  // Finalize
  void CDF_2005_S6217184::finalize() {  }
  
  
}
