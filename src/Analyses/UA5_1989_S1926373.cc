// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {


  class UA5_1989_S1926373 : public Analysis {
  public:
 
    /// Constructor
    UA5_1989_S1926373() : Analysis("UA5_1989_S1926373") {
      setBeams(PROTON, ANTIPROTON);
      _sumWPassed = 0;
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and projections
    void init() {
      addProjection(TriggerUA5(), "Trigger");
      addProjection(Beam(), "Beams");
      addProjection(ChargedFinalState(-0.5, 0.5), "CFS05");
      addProjection(ChargedFinalState(-1.5, 1.5), "CFS15");
      addProjection(ChargedFinalState(-3.0, 3.0), "CFS30");
      addProjection(ChargedFinalState(-5.0, 5.0), "CFS50");

      // NB. _hist_nch{200,900} and _hist_nch{200,900}eta50 use the same data but different binning
      _hist_nch200       = bookHistogram1D(1, 1, 1);
      _hist_nch900       = bookHistogram1D(2, 1, 1);
      _hist_nch200eta05  = bookHistogram1D(3, 1, 1);
      _hist_nch200eta15  = bookHistogram1D(4, 1, 1);
      _hist_nch200eta30  = bookHistogram1D(5, 1, 1);
      _hist_nch200eta50  = bookHistogram1D(6, 1, 1);
      _hist_nch900eta05  = bookHistogram1D(7, 1, 1);
      _hist_nch900eta15  = bookHistogram1D(8, 1, 1);
      _hist_nch900eta30  = bookHistogram1D(9, 1, 1);
      _hist_nch900eta50  = bookHistogram1D(10, 1, 1);
      _hist_mean_nch_200 = bookHistogram1D(11, 1, 1);
      _hist_mean_nch_900 = bookHistogram1D(12, 1, 1);

      /// @todo Moments of distributions
    }
 
 
    /// Do the analysis
    void analyze(const Event& event) {
      // Trigger
      const TriggerUA5& trigger = applyProjection<TriggerUA5>(event, "Trigger");
      if (!trigger.nsdDecision()) vetoEvent;

      const double sqrtS = applyProjection<Beam>(event, "Beams").sqrtS();
      const double weight = event.weight();
      _sumWPassed += weight;
   
      // Count final state particles in several eta regions
      const int numP05 = applyProjection<ChargedFinalState>(event, "CFS05").size();
      const int numP15 = applyProjection<ChargedFinalState>(event, "CFS15").size();
      const int numP30 = applyProjection<ChargedFinalState>(event, "CFS30").size();
      const int numP50 = applyProjection<ChargedFinalState>(event, "CFS50").size();
   
      // Fill histograms
      if (fuzzyEquals(sqrtS/GeV, 200.0, 1E-4)) {
        _hist_nch200->fill(numP50, weight);
        _hist_nch200eta05->fill(numP05, weight);
        _hist_nch200eta15->fill(numP15, weight);
        _hist_nch200eta30->fill(numP30, weight);
        _hist_nch200eta50->fill(numP50, weight);
        _hist_mean_nch_200->fill(_hist_mean_nch_200->binMean(0), numP50);
      }
      else if (fuzzyEquals(sqrtS/GeV, 900.0, 1E-4)) {
        _hist_nch900->fill(numP50, weight);
        _hist_nch900eta05->fill(numP05, weight);
        _hist_nch900eta15->fill(numP15, weight);
        _hist_nch900eta30->fill(numP30, weight);
        _hist_nch900eta50->fill(numP50, weight);
        _hist_mean_nch_900->fill(_hist_mean_nch_900->binMean(0), numP50);
      }
    }
 
 
 
    void finalize() {
      scale(_hist_nch200, _sumWPassed);
      scale(_hist_nch900, _sumWPassed);
      scale(_hist_nch200eta05, _sumWPassed);
      scale(_hist_nch200eta15, _sumWPassed);
      scale(_hist_nch200eta30, _sumWPassed);
      scale(_hist_nch200eta50, _sumWPassed);
      scale(_hist_nch900eta05, _sumWPassed);
      scale(_hist_nch900eta15, _sumWPassed);
      scale(_hist_nch900eta30, _sumWPassed);
      scale(_hist_nch900eta50, _sumWPassed);
      scale(_hist_mean_nch_200, 1.0/_sumWPassed);
      scale(_hist_mean_nch_900, 1.0/_sumWPassed);
    }

    //@}


  private:

    /// @name Counters
    //@{
    double _sumWPassed;
    //@}

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _hist_nch200;
    AIDA::IHistogram1D* _hist_nch900;
    AIDA::IHistogram1D* _hist_nch200eta05;
    AIDA::IHistogram1D* _hist_nch200eta15;
    AIDA::IHistogram1D* _hist_nch200eta30;
    AIDA::IHistogram1D* _hist_nch200eta50;
    AIDA::IHistogram1D* _hist_nch900eta05;
    AIDA::IHistogram1D* _hist_nch900eta15;
    AIDA::IHistogram1D* _hist_nch900eta30;
    AIDA::IHistogram1D* _hist_nch900eta50;
    AIDA::IHistogram1D* _hist_mean_nch_200;
    AIDA::IHistogram1D* _hist_mean_nch_900;
    //@}

  };



  // This global object acts as a hook for the plugin system
  AnalysisBuilder<UA5_1989_S1926373> plugin_UA5_1989_S1926373;

}
