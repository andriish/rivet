// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/CentralityBinner.hh"

namespace Rivet {

/// Example of CentralityEstimator projection that uses summed Et in
/// the forward region.
class SumETFwdCentrality: public CentralityEstimator {

public:

  /// Constructor.
  SumETFwdCentrality() {
    declare(FinalState(3.2, 4.9, 100*MeV), "FSSumETFwdCentrality");
  }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SumETFwdCentrality);

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
  
  /// Compare projections
  int compare(const Projection& p) const {
    return mkNamedPCmp(p, "FSSumETFwdCentrality");
  }

};
    

  /// Generic analysis looking at various distributions of final state particles
  class MC_Centrality : public Analysis {
  public:

    /// Constructor
    MC_Centrality() : Analysis("MC_Centrality"), _cent(200, 0.02) {}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Just testing
      CentralityBinner< tuple<Histo1DPtr,Histo2DPtr,Profile1DPtr,Profile2DPtr> >
        testbinner;

      // Projections
      declare(ChargedFinalState(-2.5, 2.5, 500*MeV), "CFScent");
      declare(FinalState(3.2, 4.9, 100*MeV), "FSfwd");
      

      // Histograms
      // The sum Et in the forward region.
      _histETfwd  = bookHisto1D("ETfwd", 200, 0.0, 100.0);
      
      // The overall charged multiplicity distribution as a function
      // of eta.
      _histEtaAll = bookHisto1D("EtaAll", 50, -2.5, 2.5);

      // Multiplicity distribution for different centrality intervals.
      // The notation is that the maximum centrality is 100% and the
      // lowest is 0%. The _cent1 object will dynamically decide how to
      // cluster fills to the histogram together.
      _histETfwdC95 = bookHisto1D("ETfwdC95", 200, 0.0, 100.0);
      _histETfwdC90 = bookHisto1D("ETfwdC90", 200, 0.0, 100.0);
      _histETfwdC80 = bookHisto1D("ETfwdC80", 200, 0.0, 100.0);
      _histETfwdC60 = bookHisto1D("ETfwdC60", 200, 0.0, 100.0);
      _histETfwdC00 = bookHisto1D("ETfwdC00", 200, 0.0, 100.0);
      _histEtaC95 = bookHisto1D("EtaC95", 50, -2.5, 2.5);
      _histEtaC90 = bookHisto1D("EtaC90", 50, -2.5, 2.5);
      _histEtaC80 = bookHisto1D("EtaC80", 50, -2.5, 2.5);
      _histEtaC60 = bookHisto1D("EtaC60", 50, -2.5, 2.5);
      _histEtaC00 = bookHisto1D("EtaC00", 50, -2.5, 2.5);
      _cent.add(make_tuple(_histETfwdC95, _histEtaC95), 95.0, 100.0);
      _cent.add(make_tuple(_histETfwdC90, _histEtaC90), 90.0,  95.0);
      _cent.add(make_tuple(_histETfwdC80, _histEtaC80), 80.0,  90.0);
      _cent.add(make_tuple(_histETfwdC60, _histEtaC60), 60.0,  80.0);
      _cent.add(make_tuple(_histETfwdC00, _histEtaC00),  0.0,  60.0);

      _centrue.clear();

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

      _histETfwd->fill(sumet, weight);
      _centrue.insert(make_pair(sumet, weight));


      // Also setup other centrality binners.
      auto ch = _cent.select(sumet, weight);
      std::get<0>(ch)->fill(sumet, weight);

     // Then fill the selected centrality histogram.
      const ChargedFinalState & cfscent =
        apply<ChargedFinalState>(event, "CFScent");
      for ( const Particle & p : cfscent.particles() ) {
        std::get<1>(ch)->fill(p.eta(), weight);
	_histEtaAll->fill(p.eta(), weight);
      }

    }


    /// Finalize
    void finalize() {
        
      normalize(_histETfwd, _histETfwd->sumW()/sumOfWeights());
      normalize(_histEtaAll, _histEtaAll->sumW()/sumOfWeights());

      _cent.finalize();
      _cent.normalizePerEvent();

      map<double,double> edges = _cent.edges();
      map<double,double> edges0 = edges;

      auto curr = edges0.begin();
      curr->second = 0.0;
      ++curr;
      pair<double,double> prev = *_centrue.begin();
      double acc = 0.0;
      for ( auto next: _centrue ) {
	double del = next.second/sumOfWeights();
	if ( acc + del > curr->first ) {
	  curr->second = prev.first + (curr->first - acc)*(next.first - prev.first)/del;
	  ++curr;
	}
	prev = next;
	acc += del;
      }


      MSG_INFO("Cross check of centalty edges:\n");
      MSG_INFO("" << setw(10) << "%" << setw(10) << "true"
               << setw(10) << "estimate");
      for ( auto e : edges0 ) {
        MSG_INFO(""
		 << setw(10) << e.first << setw(10) << e.second
		 << setw(10) << edges[e.first]);
      }
      

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histETfwd;
    Histo1DPtr _histEtaAll;
    Histo1DPtr _histETfwdC95, _histETfwdC90, _histETfwdC80,
      _histETfwdC60, _histETfwdC00;
    Histo1DPtr _histEtaC95, _histEtaC90, _histEtaC80,
      _histEtaC60, _histEtaC00;
    CentralityBinner< tuple<Histo1DPtr,Histo1DPtr> > _cent;
    //@}

    /// Keep track of the actually generated centralities.
    multimap<double, double> _centrue;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Centrality);

}
