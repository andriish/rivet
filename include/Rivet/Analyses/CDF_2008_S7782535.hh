// -*- C++ -*-
#ifndef RIVET_CDF_2008_S7782535_HH
#define RIVET_CDF_2008_S7782535_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/Projections/JetShape.hh"

namespace Rivet {


  /// Implementation of CDF RunII b-jet shape paper
  /// @todo Test with Pythia
  class CDF_2008_S7782535 : public Analysis {

  public:

    /// Constructor.

    /// jet cuts: |eta| <= 0.7

    /// Don't attempt to model the following cuts :
    ///   Missing ET significance
    ///   veto on additional vertices
    ///   Zvtx < 50

    CDF_2008_S7782535()
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

  public:

    /// Factory method
    static Analysis* create() { 
      return new CDF_2008_S7782535(); 
    }


    /// @name Publication metadata
    //@{
    /// Get SPIRES ID code.
    string spiresId() const {
      return "7782535";
    }
    /// A short description of the analysis.
    string summary() const {
      return "CDF Run II b-jet shape paper";
    }

/// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "A measurement of the shapes of b-jets using 300 pb$^{-1}$ "
	 << "of data obtained with CDF II in $p\bar{p}$ collisions at "
	 << "$\sqrt{s}=1.96$ TeV. The measured quantity is the average "
	 << "integrated jet shape, which is computed over an ensemble of "
	 << "jets. This quantity is expressed as "
	 << "$\Psi(r/R) = \langle\frac{p_T(0 \rightarrow r)}"
         << "{p_T(0 \rightarrow R)}\rangle$, "
         << "where pT (0$\rightarrow$r) is the scalar sum of the "
	 << "transverse momenta of all objects inside a sub-cone of "
	 << "radius r around the jet axis. The integrated shapes are "
	 << "by definition normalized such that  $\Psi(r/R  =1) =1$. "
         << "\n "   
         << "The measurement is done in bins of jet pT in the range "
         << "52 to 300 GeV/c. The jets have $|\eta| < 0.7$. "
         << "The b-jets are expected to be broader than inclusive jets."
         << "Moreover, b-jets containing a single b-quark are expected "
         << "to be narrower than those containing a b bbar pair from "
         << "gluon splitting. ";
      return os.str();
    }
  /// Characteristics of events to be processed by this analysis
    string runInfo() const {
      ostringstream os;
      os << "Requires 2->2 QCD scattering processes."
         << "The minimum jet Et is 52 GeV, so a cut on kinematic pTmin"
         << "may be required for good statistics.";
      return os.str();
    }

    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }

  /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 2";
    }

    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2008";
    }

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn += "Alison Lister <alister@fnal.gov>";
      rtn += "Emily Nurse <nurse@hep.ucl.ac.uk>";
      return rtn;
    }

    /// Journal, and preprint references.
    virtual vector<string> references() const {
      vector<string> ret;
      ret.push_back("arXiv:0806.1699");
     ret.push_back("Phys.Rev.D78:072005,2008");
      return ret;
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}


  private:

    /// @name Analysis cuts

    /// @name Analysis data
    //@{
    vector<FourMomentum> _jetaxes;

    double _Rjet;
    vector<double> _pTbins;
    int _NpTbins;
    //@{
    /// Histograms
    AIDA::IProfile1D* _Psi_pT[4];
    AIDA::IDataPointSet* _OneMinusPsi_vs_pT;
    //@}

  };

}

#endif
