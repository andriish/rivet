// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"

namespace Rivet {


  /// @brief Just measures a few random things as an example.
  class EXAMPLE : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(EXAMPLE);

    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {
      // Projections
      MSG_TRACE(0);
      const FinalState cnfs(Cuts::abseta < 4 && Cuts::pT > 500*MeV);
      MSG_TRACE(1);
      const ChargedFinalState cfs(cnfs);
      MSG_TRACE(2);
      declare(cnfs, "FS");
      MSG_TRACE(3);
      declare(cfs, "CFS");
      declare(FastJets(cnfs, FastJets::KT, 0.7), "Jets");
      declare(Thrust(cfs), "Thrust");
      declare(Sphericity(cfs), "Sphericity");

      // Histograms
      _histTot         = bookHisto1D("TotalMult", 100, -0.5, 99.5);
      _histChTot       = bookHisto1D("TotalChMult", 50, -1.0, 99.0);
      _histHadrTot     = bookHisto1D("HadrTotalMult", 100, -0.5, 99.5);
      _histHadrChTot   = bookHisto1D("HadrTotalChMult", 50, -1.0, 99.0);
      _histMajor       = bookHisto1D("Major", 10, 0.0, 0.6);
      _histSphericity  = bookHisto1D("Sphericity", 10, 0.0, 0.8);
      _histAplanarity  = bookHisto1D("Aplanarity", 10, 0.0, 0.3);

      // Non-uniform binning example:
      double edges[11] = { 0.5, 0.6, 0.7, 0.80, 0.85, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0 };
      vector<double> vedges(edges, edges+11);
      _histThrust = bookHisto1D("Thrust", vedges);
    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Make sure to always include the event weight in histogram fills!
      const double weight = 1.0;

      const Particles& cnparticles = apply<FinalState>(event, "FS").particles();
      MSG_DEBUG("Total multiplicity = " << cnparticles.size());
      _histTot->fill(cnparticles.size(), weight);
      int cnhadronmult = 0;
      foreach (const Particle& p, cnparticles) if (isHadron(p)) cnhadronmult += 1;
      MSG_DEBUG("Hadron multiplicity = " << cnhadronmult);
      _histHadrTot->fill(cnhadronmult, weight);

      const Particles& cparticles = apply<FinalState>(event, "CFS").particles();
      MSG_DEBUG("Total charged multiplicity = " << cparticles.size());
      _histChTot->fill(cparticles.size(), weight);
      int chadronmult = 0;
      foreach (const Particle& p, cparticles) if (isHadron(p)) chadronmult += 1;
      MSG_DEBUG("Hadron charged multiplicity = " << chadronmult);
      _histHadrChTot->fill(chadronmult, weight);

      const Thrust& t = apply<Thrust>(event, "Thrust");
      MSG_DEBUG("Thrust = " << t.thrust());
      _histThrust->fill(t.thrust(), weight);
      _histMajor->fill(t.thrustMajor(), weight);

      const Sphericity& s = apply<Sphericity>(event, "Sphericity");
      MSG_DEBUG("Sphericity = " << s.sphericity());
      _histSphericity->fill(s.sphericity(), weight);
      MSG_DEBUG("Aplanarity = " << s.aplanarity());
      _histAplanarity->fill(s.aplanarity(), weight);

      const Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::pT > 5*GeV);
      size_t num_b_jets = count_if(jets.begin(), jets.end(), [](const Jet& j){ return j.bTagged(Cuts::pT > 500*MeV); });
      MSG_DEBUG("Num B-jets with pT > 5 GeV = " << num_b_jets);
    }


    /// Finalize
    void finalize() {
      normalize(_histTot);
      normalize(_histChTot);
      normalize(_histHadrTot);
      normalize(_histHadrChTot);
      normalize(_histThrust);
      normalize(_histMajor);
      normalize(_histSphericity);
      normalize(_histAplanarity);
    }

    //@}


  private:

    //@{
    /// Histograms
    Histo1DPtr _histTot;
    Histo1DPtr _histChTot;
    Histo1DPtr _histHadrTot;
    Histo1DPtr _histHadrChTot;
    Histo1DPtr _histThrust;
    Histo1DPtr _histMajor;
    Histo1DPtr _histSphericity;
    Histo1DPtr _histAplanarity;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EXAMPLE);

}
