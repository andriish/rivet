// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// @brief CDF Run II jet shape analysis
  /// @author Andy Buckley
  class CDF_2005_S6217184 : public Analysis {
  public:

    /// Constructor
    CDF_2005_S6217184()
      : Analysis("CDF_2005_S6217184")
    {
      setBeams(PROTON, ANTIPROTON);
    }


    /// @name Analysis methods
    //@{

    void init() {
      // Set up projections
      const FinalState fs(-2.0, 2.0);
      addProjection(fs, "FS");
      FastJets fj(fs, FastJets::CDFMIDPOINT, 0.7);
      fj.useInvisibles();
      addProjection(fj, "Jets");

      // Specify pT bins
      _ptedges += 37.0, 45.0, 55.0, 63.0, 73.0, 84.0, 97.0, 112.0, 128.0,
          148.0, 166.0, 186.0, 208.0, 229.0, 250.0, 277.0, 304.0, 340.0, 380.0;

      // Register a jet shape projection and histogram for each pT bin
      for (size_t i = 0; i < 6; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          const size_t k = i*3 + j;
          stringstream ss; ss << "JetShape" << k;
          const string pname = ss.str();
          _jsnames_pT[k] = pname;
          const JetShape jsp(fj, 0.0, 0.7, 6, _ptedges[k], _ptedges[k+1], 0.1, 0.7, RAPIDITY);
          addProjection(jsp, pname);
          _profhistRho_pT[k] = bookProfile1D(i+1, 1, j+1);
          _profhistPsi_pT[k] = bookProfile1D(6+i+1, 1, j+1);
        }
      }

      // Final histo
      _profhistPsi = bookProfile1D(13, 1, 1);
    }



    /// Do the analysis
    void analyze(const Event& evt) {

      // Get jets and require at least one to pass pT and y cuts
      const Jets jets = applyProjection<FastJets>(evt, "Jets").jetsByPt(37, 380, 0.1, 0.7);
      MSG_DEBUG("Jet multiplicity before cuts = " << jets.size());
      if (jets.size() == 0) {
        MSG_DEBUG("No jets found in required pT range");
        vetoEvent;
      }

      // Calculate and histogram jet shapes
      const double weight = evt.weight();

      for (size_t ipt = 0; ipt < 18; ++ipt) {
        const JetShape& jsipt = applyProjection<JetShape>(evt, _jsnames_pT[ipt]);
        for (size_t rbin = 0; rbin < jsipt.numBins(); ++rbin) {
          const double r_rho = jsipt.rBinMid(rbin);
          _profhistRho_pT[ipt]->fill(r_rho/0.7, (0.7/1.0)*jsipt.diffJetShape(rbin), weight);
          const double r_Psi = jsipt.rBinMax(rbin);
          _profhistPsi_pT[ipt]->fill(r_Psi/0.7, jsipt.intJetShape(rbin), weight);
        }

        // Final histo is Psi(0.3/R) as a function of jet pT bin
        /// @todo Can this actually be calculated event by event, or does it need to be assembled in a DPS in finalize? See CDF_2008_S7782535.
        const double ptmid = (_ptedges[ipt] + _ptedges[ipt+1])/2.0;
        _profhistPsi->fill(ptmid/GeV, jsipt.intJetShape(2), weight);
      }

    }


    // Finalize
    void finalize() {
      //
    }

    //@}


  private:

    /// @name Analysis data
    //@{

    /// Jet \f$ p_\perp\f$ bins.
    vector<double> _ptedges; // This can't be a raw array if we want to initialise it non-painfully

    /// JetShape projection name for each \f$p_\perp\f$ bin.
    string _jsnames_pT[18];

    //@}


    /// @name Histograms
    //@{
    AIDA::IProfile1D* _profhistRho_pT[18];
    AIDA::IProfile1D* _profhistPsi_pT[18];
    AIDA::IProfile1D* _profhistPsi;
    //@}

  };


  // This global object acts as a hook for the plugin system
  AnalysisBuilder<CDF_2005_S6217184> plugin_CDF_2005_S6217184;

}
