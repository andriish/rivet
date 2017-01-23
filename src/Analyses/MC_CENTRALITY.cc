// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/CentralityHistogram.hh"

namespace Rivet {

/// Example of CentralityEstimator projection that uses summed Et in
/// the forward region.
class SumETFwdCentrality: public CentralityEstimator {

public:

  /// Constructor.
  SumETFwdCentrality() {
    declare(FinalState(3.2, 4.9, 100*MeV), "FSSumETFwdCentrality");
  }

protected:

  /// Perform the projection on the Event
  void project(const Event& e) {
    const FinalState & fsfwd =
      apply<FinalState>(e, "FSSumETFwdCentrality");
    _estimate = 0.0;
    for ( const Particle & p : fsfwd.particles() ) {
      _estimate += p.Et();
    }
  }
};
    

  /// Generic analysis looking at various distributions of final state particles
  class MC_CENTRALITY : public Analysis {
  public:

    /// Constructor
    MC_CENTRALITY() : Analysis("MC_CENTRALITY") {}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      declare(ChargedFinalState(-2.5, 2.5, 500*MeV), "CFScent");
      declare(FinalState(3.2, 4.9, 100*MeV), "FSfwd");
      

      // Histograms
      // The sum Et in the forward region.
      _histETfwd  = bookHisto1D("ETfwd", 100, 0.0, 500.0);
      
      // The overall charged multiplicity distribution as a function
      // of eta.
      _histEtaAll = bookHisto1D("EtaAll", 50, -2.5, 2.5);

      // Multiplicity distribution for different centrality intervals.
      // The notation is that the maximum centrality is 100% and the
      // lowest is 0%. THe _cent object will dynamically decide how to
      // cluster fills to the histogram together.
      _cent.add(_histEtaC95 = bookHisto1D("EtaC95", 50, -2.5, 2.5),
                95.0, 100.0);
      _cent.add(_histEtaC90 = bookHisto1D("EtaC90", 50, -2.5, 2.5),
                90.0,  95.0);
      _cent.add(_histEtaC80 = bookHisto1D("EtaC80", 50, -2.5, 2.5),
                80.0,  90.0);
      _cent.add(_histEtaC60 = bookHisto1D("EtaC60", 50, -2.5, 2.5),
                60.0,  80.0);
      _cent.add(_histEtaC00 = bookHisto1D("EtaC00", 50, -2.5, 2.5),
                 0.0,  60.0);

      // Try also the facility of using a CentralityEstimator
      // projection.
      declare(SumETFwdCentrality(), "SumETFwdCentrality");
      _cent2.setProjection("SumETFwdCentrality");
      _cent2.add(_histEta2C95 = bookHisto1D("Eta2C95", 50, -2.5, 2.5),
                95.0, 100.0);
      _cent2.add(_histEta2C90 = bookHisto1D("Eta2C90", 50, -2.5, 2.5),
                90.0,  95.0);
      _cent2.add(_histEta2C80 = bookHisto1D("Eta2C80", 50, -2.5, 2.5),
                80.0,  90.0);
      _cent2.add(_histEta2C60 = bookHisto1D("Eta2C60", 50, -2.5, 2.5),
                60.0,  80.0);
      _cent2.add(_histEta2C00 = bookHisto1D("Eta2C00", 50, -2.5, 2.5),
                 0.0,  60.0);

      /// Finally also try the case where we know the centrality bin
      /// edges beforehand.
      _cent3.setProjection("SumETFwdCentrality");
      _cent3.add(_histEta3C95 = bookHisto1D("Eta3C95", 50, -2.5, 2.5),
                 95.0, 100.0, 22.55, 1.0e20);
      _cent3.add(_histEta3C90 = bookHisto1D("Eta3C90", 50, -2.5, 2.5),
                 90.0,  95.0, 16.11, 22.55);
      _cent3.add(_histEta3C80 = bookHisto1D("Eta3C80", 50, -2.5, 2.5),
                 80.0,  90.0, 9.611, 16.11);
      _cent3.add(_histEta3C60 = bookHisto1D("Eta3C60", 50, -2.5, 2.5),
                 60.0,  80.0, 3.921, 9.611);
      _cent3.add(_histEta3C00 = bookHisto1D("Eta3C00", 50, -2.5, 2.5),
                 0.0,  60.0, 0.0, 3.921);



    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // First calculate the centrality estimator. In this case the
      // sumed transverse enegry in a forward region of the event, and
      // initialize the centrality bin.
      const FinalState & fsfwd = apply<FinalState>(event, "FSfwd");
      double sumet = 0.0;
      for ( const Particle & p : fsfwd.particles() ) {
        sumet += p.Et();
      }
      _cent.setup(event, sumet);
      _histETfwd->fill(sumet, weight);

      // Also setup the other centrality binner.
      _cent2.setup(event);

      // Then fill the selected centrality histogram.
      const ChargedFinalState & cfscent =
        apply<ChargedFinalState>(event, "CFScent");
      for ( const Particle & p : cfscent.particles() ) {
        _cent.fill(p.eta(), weight);
        _cent2.fill(p.eta(), weight);
        _cent3.fill(p.eta(), weight);
      }
    }


    /// Finalize
    void finalize() {
      normalize(_histETfwd, 1.0/sumOfWeights());
      _cent.finalize();
      _cent2.finalize();
      _cent3.finalize(); // This should be a noop.
      _cent.normalizePerEvent();
      _cent2.normalizePerEvent();
      _cent3.normalizePerEvent();
      _histEtaAll->normalize(1.0/sumOfWeights());

      map<double,double> edges1 = _cent.edges();
      map<double,double> edges2 = _cent2.edges();
      map<double,double> edges3 = _cent3.edges();
      MSG_INFO("Cross check of centalty edges:"
               << setw(10) << "%" << setw(10) << "1"
               << setw(10) << "2" << setw(10) << "3");
      for ( auto e : edges3 ) {
        MSG_INFO(""
                 << setw(10) << e.first << setw(10) << e.second << setw(10) 
                 << edges1[e.first] << setw(10) << edges2[e.first]);
      }
      

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histETfwd;
    Histo1DPtr _histEtaAll;
    Histo1DPtr _histEtaC95, _histEtaC90, _histEtaC80, _histEtaC60, _histEtaC00;
    CentralityHistogram _cent;
    Histo1DPtr _histEta2C95, _histEta2C90, _histEta2C80,
    _histEta2C60, _histEta2C00;
    CentralityHistogram _cent2;
    Histo1DPtr _histEta3C95, _histEta3C90, _histEta3C80,
    _histEta3C60, _histEta3C00;
    CentralityHistogram _cent3;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_CENTRALITY);

}
