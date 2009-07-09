// -*- C++ -*-
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/CDF_2005_S6217184.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  // Constructor
  CDF_2005_S6217184::CDF_2005_S6217184()
    : Analysis("CDF_2005_S6217184")
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

    // Specify pT bins
    _pTbins += 37.0, 45.0, 55.0, 63.0, 73.0, 84.0, 97.0, 112.0, 128.0, 
      148.0, 166.0, 186.0, 208.0, 229.0, 250.0, 277.0, 304.0, 340.0, 380.0;
  }


  
  // Book histograms
  void CDF_2005_S6217184::init() {
    // 18 = 6x3 pT bins, one histogram each
    for (size_t i = 0; i < 6; ++i) { 
      for (size_t j = 0; j < 3; ++j) {
        size_t k = i*3 + j;
        stringstream ptrange;
        ptrange << "$" << _pTbins[k] << "\\,\\text{GeV}/c < p_\\perp^\\text{jet} < " 
                << _pTbins[k+1] << "\\,\\text{GeV}/c$";

        stringstream difftitle;
        difftitle << "Differential jet shape $\\rho$, " << ptrange.str();
        _profhistRho_pT[k] = 
          bookProfile1D(i+1, 1, j+1, difftitle.str(), "$r/R$", "$\\rho(r/R)$");

        stringstream inttitle;
        inttitle << "Integral jet shape $\\Psi$, " << ptrange.str();
        _profhistPsi_pT[k] = 
          bookProfile1D(6+i+1, 1, j+1, inttitle.str(), "$r/R$", "$\\Psi(r/R)$");
      }
    }    

    _profhistPsi = 
      bookProfile1D(13, 1, 1, "Integral jet shape, $\\Psi$(0.3/$R$), vs. $p_\\perp^\\text{jet}$",
                    "$p_\\perp^\\text{jet}$ / GeV/$c$", "$\\Psi(0.3/R)$");
  }
  
  
  
  // Do the analysis
  void CDF_2005_S6217184::analyze(const Event& event) {

    // Get jets and require at least one to pass pT and y cuts
    const Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt();
    getLog() << Log::DEBUG << "Jet multiplicity before cuts = " << jets.size() << endl;

    // Determine the central jet axes
    _jetaxes.clear();
    foreach (const Jet& jt, jets) {
      const FourMomentum pj = jt.momentum();
      if (inRange(pj.pT()/GeV, 37.0, 380.0) && inRange(fabs(pj.rapidity()), 0.1, 0.7)) {
        _jetaxes.push_back(jt.momentum());
      }
    }
    if (_jetaxes.empty()) vetoEvent;

    // Calculate and histogram jet shapes
    const double weight = event.weight();
    const JetShape& js = applyProjection<JetShape>(event, "JetShape");

    /// @todo Use BinnedHistogram, for collections of histos each for a range of values of an extra variable
    for (size_t jind = 0; jind < _jetaxes.size(); ++jind) {
      for (size_t ipT = 0; ipT < 18; ++ipT) {
        if (_jetaxes[jind].pT() > _pTbins[ipT] && _jetaxes[jind].pT() <= _pTbins[ipT+1]) {
          for (size_t rbin = 0; rbin < js.numBins(); ++rbin) {
            const double rad_Rho = js.rMin() + (rbin+0.5)*js.interval();
            _profhistRho_pT[ipT]->fill(rad_Rho/0.7, (0.7/1.0)*js.diffJetShape(jind, rbin), weight);
            /// @todo Calc int histos from diff histos
            const double rad_Psi = js.rMin() +(rbin+1.0)*js.interval();
            _profhistPsi_pT[ipT]->fill(rad_Psi/0.7, js.intJetShape(jind, rbin), weight);
          }
          /// @todo Calc int histos from diff histos
          _profhistPsi->fill((_pTbins[ipT] + _pTbins[ipT+1])/2.0, js.psi(jind), weight);
        }
      }
    }

  }
  

  // Finalize
  void CDF_2005_S6217184::finalize() {  
    //
  }
  
  
}
