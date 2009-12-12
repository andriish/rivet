// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// Implementation of CDF RunII b-jet shape paper
  class CDF_2008_S7782535 : public Analysis {
  public:

    /// Constructor
    CDF_2008_S7782535() : Analysis("CDF_2008_S7782535")
    {
      setBeams(PROTON, ANTIPROTON);
      _Rjet = 0.7;
      _NpTbins = 4;
    }
 

    /// @name Analysis methods
    //@{

    void init() {
      // Set up projections
      const FinalState fs(-3.6, 3.6);
      addProjection(fs, "FS");
      // Veto (anti)neutrinos, and muons with pT above 1.0 GeV
      VetoedFinalState vfs(fs);
      vfs.vetoNeutrinos();
      vfs.addVetoPairDetail(MUON, 1.0*GeV, MAXDOUBLE);
      addProjection(vfs, "VFS");
      addProjection(FastJets(vfs, FastJets::CDFMIDPOINT, 0.7), "Jets");
      addProjection(JetShape(vfs, _jetaxes, 0.0, 0.7, 0.1, 0.3), "JetShape");

      // Book histograms
      _pTbins += 52, 80, 104, 142, 300;
      for (int i = 0; i < _NpTbins; ++i) {
        _h_Psi_pT[i] = bookProfile1D(i+1, 2, 1);
      }
      _h_OneMinusPsi_vs_pT = bookDataPointSet(5, 1, 1);
    }
 
 
    // Do the analysis
    void analyze(const Event& event) {
      // Get jets
      const Jets& jets = applyProjection<FastJets>(event, "Jets").jetsByPt();
      getLog() << Log::DEBUG << "Jet multiplicity before any pT cut = " << jets.size() << endl;
   
      // Determine the central jet axes
      _jetaxes.clear();
      foreach (const Jet& j, jets) {
        if (j.containsBottom()) {
          // Only use central calorimeter jets
          FourMomentum pjet = j.momentum();
          if (pjet.pT()/GeV > _pTbins[0] && fabs(pjet.rapidity()) < 0.7) {
            _jetaxes.push_back(pjet);
          }
        }
      }
      if (_jetaxes.empty())  {
        getLog() << Log::DEBUG << "No b-jet axes in acceptance" << endl;
        vetoEvent;
      }
   
      // Determine jet shapes
      const JetShape& js = applyProjection<JetShape>(event, "JetShape");
   
      /// @todo Replace with foreach
      for (size_t jind = 0; jind < _jetaxes.size(); ++jind) {
        // Put jet in correct pT bin
        int jet_pt_bin = -1;
        for (size_t i = 0; i < 4; ++i) {
          if (inRange(_jetaxes[jind].pT(), _pTbins[i], _pTbins[i+1])) {
            jet_pt_bin = i;
            break;
          }
        }
        if (jet_pt_bin > -1) {
          // Fill each entry in profile
          for (size_t rbin = 0; rbin < js.numBins(); ++rbin) {
            const double rad_Psi = js.rMin() + (rbin+1.0)*js.interval();
            /// @todo Yuck... JetShape's interface sucks
            _h_Psi_pT[jet_pt_bin]->fill(rad_Psi/_Rjet, js.intJetShape(jind, rbin), event.weight() );
          }
        }
      }
   
    }
 
 
    /// Finalize
    void finalize() {
      vector<double> y, ey;
      for (size_t i = 0; i < _pTbins.size()-1; ++i) {
        // Get entry for rad_Psi = 0.2 bin
        AIDA::IProfile1D* ph_i = _h_Psi_pT[i];
        y.push_back(1.0 - ph_i->binHeight(1));
        ey.push_back(ph_i->binError(1));
      }
      _h_OneMinusPsi_vs_pT->setCoordinate(1, y, ey);
    }
 
    //@}


  private:

    /// @name Analysis data
    //@{
    vector<FourMomentum> _jetaxes;
    double _Rjet;
    vector<double> _pTbins;
    int _NpTbins;
    //@}


    /// @name Histograms
    //@{
    AIDA::IProfile1D* _h_Psi_pT[4];
    AIDA::IDataPointSet* _h_OneMinusPsi_vs_pT;
    //@}

  };


  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_2008_S7782535> plugin_CDF_2008_S7782535;

}
