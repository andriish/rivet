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

      // Projections
      declare(ChargedFinalState(-2.5, 2.5, 500*MeV), "CFScent");
      declare(FinalState(3.2, 4.9, 100*MeV), "FSfwd");

      // Histograms
      // The sum Et in the forward region.
      _hETfwd  = bookHisto1D("ETfwd", 200, 0.0, 100.0);
      
      // The overall charged multiplicity distribution as a function
      // of eta.
      _hEtaAll = bookHisto1D("EtaAll", 50, -2.5, 2.5);

      // Distribution for different centrality intervals.  The
      // notation is that the maximum centrality is 100% and the
      // lowest is 0%. The _cent object will dynamically decide how to
      // cluster fills to the histogram together.
      _hETfwdC95 = bookHisto1D("ETfwdC95", 200, 0.0, 100.0);
      _hETfwdC90 = bookHisto1D("ETfwdC90", 200, 0.0, 100.0);
      _hETfwdC80 = bookHisto1D("ETfwdC80", 200, 0.0, 100.0);
      _hETfwdC60 = bookHisto1D("ETfwdC60", 200, 0.0, 100.0);
      _hETfwdC00 = bookHisto1D("ETfwdC00", 200, 0.0, 100.0);
      _hEtaC95 = bookHisto1D("EtaC95", 50, -2.5, 2.5);
      _hEtaC90 = bookHisto1D("EtaC90", 50, -2.5, 2.5);
      _hEtaC80 = bookHisto1D("EtaC80", 50, -2.5, 2.5);
      _hEtaC60 = bookHisto1D("EtaC60", 50, -2.5, 2.5);
      _hEtaC00 = bookHisto1D("EtaC00", 50, -2.5, 2.5);
      _cent.add(make_tuple(_hETfwdC95, _hEtaC95), 95.0, 100.0);
      _cent.add(make_tuple(_hETfwdC90, _hEtaC90), 90.0,  95.0);
      _cent.add(make_tuple(_hETfwdC80, _hEtaC80), 80.0,  90.0);
      _cent.add(make_tuple(_hETfwdC60, _hEtaC60), 60.0,  80.0);
      _cent.add(make_tuple(_hETfwdC00, _hEtaC00),  0.0,  60.0);
      // Distribution for different centrality intervals as reported by HepMC3.  The
      // notation is that the maximum centrality is 100% and the
      // lowest is 0%.
      _hETfwdGC95 = bookHisto1D("ETfwdGC95", 200, 0.0, 100.0);
      _hETfwdGC90 = bookHisto1D("ETfwdGC90", 200, 0.0, 100.0);
      _hETfwdGC80 = bookHisto1D("ETfwdGC80", 200, 0.0, 100.0);
      _hETfwdGC60 = bookHisto1D("ETfwdGC60", 200, 0.0, 100.0);
      _hETfwdGC00 = bookHisto1D("ETfwdGC00", 200, 0.0, 100.0);
      _hEtaGC95 = bookHisto1D("EtaGC95", 50, -2.5, 2.5);
      _hEtaGC90 = bookHisto1D("EtaGC90", 50, -2.5, 2.5);
      _hEtaGC80 = bookHisto1D("EtaGC80", 50, -2.5, 2.5);
      _hEtaGC60 = bookHisto1D("EtaGC60", 50, -2.5, 2.5);
      _hEtaGC00 = bookHisto1D("EtaGC00", 50, -2.5, 2.5);
      _gencent.setProjection(GeneratedCentrality(), "GenCent");
      _gencent.add(make_tuple(_hETfwdGC95, _hEtaGC95), 95.0, 100.0, 95.0, 100.0);
      _gencent.add(make_tuple(_hETfwdGC90, _hEtaGC90), 90.0,  95.0, 90.0,  95.0);
      _gencent.add(make_tuple(_hETfwdGC80, _hEtaGC80), 80.0,  90.0, 80.0,  90.0);
      _gencent.add(make_tuple(_hETfwdGC60, _hEtaGC60), 60.0,  80.0, 60.0,  80.0);
      _gencent.add(make_tuple(_hETfwdGC00, _hEtaGC00),  0.0,  60.0,  0.0,  60.0);

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

      _hETfwd->fill(sumet, weight);
      _centrue.insert(make_pair(sumet, weight));


      // Setup centrality binners.
      auto ch = _cent.select(sumet, weight);
      auto gch = _gencent.select(event, weight);
      std::get<0>(ch)->fill(sumet, weight);
      std::get<0>(gch)->fill(sumet, weight);
      

     // Then fill the selected centrality histogram.
      const ChargedFinalState & cfscent =
        apply<ChargedFinalState>(event, "CFScent");
      for ( const Particle & p : cfscent.particles() ) {
        std::get<1>(ch)->fill(p.eta(), weight);
        std::get<1>(gch)->fill(p.eta(), weight);
	_hEtaAll->fill(p.eta(), weight);
      }

    }


    /// Finalize
    void finalize() {
        
      normalize(_hETfwd, _hETfwd->sumW()/sumOfWeights());
      normalize(_hEtaAll, _hEtaAll->sumW()/sumOfWeights());

      _cent.finalize();
      _cent.normalizePerEvent();
      _gencent.finalize();
      _gencent.normalizePerEvent();

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
    Histo1DPtr _hETfwd;
    Histo1DPtr _hEtaAll;
    Histo1DPtr _hETfwdC95, _hETfwdC90, _hETfwdC80, _hETfwdC60, _hETfwdC00;
    Histo1DPtr _hEtaC95, _hEtaC90, _hEtaC80, _hEtaC60, _hEtaC00;
    CentralityBinner< tuple<Histo1DPtr,Histo1DPtr> > _cent;
    Histo1DPtr _hETfwdGC95, _hETfwdGC90, _hETfwdGC80, _hETfwdGC60, _hETfwdGC00;
    Histo1DPtr _hEtaGC95, _hEtaGC90, _hEtaGC80, _hEtaGC60, _hEtaGC00;
    CentralityBinner< tuple<Histo1DPtr,Histo1DPtr> > _gencent;
    //@}

    /// Keep track of the actually generated centralities.
    multimap<double, double> _centrue;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Centrality);

}
