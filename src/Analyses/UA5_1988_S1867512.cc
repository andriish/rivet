// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {
  

  class UA5_1988_S1867512 : public Analysis {
  public:

    UA5_1988_S1867512() : Analysis("UA5_1988_S1867512")
    {
      setBeams(PROTON, ANTIPROTON);
    }
    
    
    /// @name Analysis methods
    //@{

    void init() {
      addProjection(TriggerUA5(), "Trigger");
      addProjection(Beam(), "Beams");
      
      // Symmetric eta interval
      addProjection(ChargedFinalState(-0.5, 0.5), "CFS05");

      // Asymmetric intervals first
      // Forward eta intervals
      addProjection(ChargedFinalState(0.0, 1.0), "CFS10F");
      addProjection(ChargedFinalState(0.5, 1.5), "CFS15F");
      addProjection(ChargedFinalState(1.0, 2.0), "CFS20F");
      addProjection(ChargedFinalState(1.5, 2.5), "CFS25F");
      addProjection(ChargedFinalState(2.0, 3.0), "CFS30F");
      addProjection(ChargedFinalState(2.5, 3.5), "CFS45F");
      addProjection(ChargedFinalState(3.0, 4.0), "CFS40F");
      
      // Backward eta intervals
      addProjection(ChargedFinalState(-1.0,  0.0), "CFS10B");
      addProjection(ChargedFinalState(-1.5, -0.5), "CFS15B");
      addProjection(ChargedFinalState(-2.0, -1.0), "CFS20B");
      addProjection(ChargedFinalState(-2.5, -1.5), "CFS25B");
      addProjection(ChargedFinalState(-3.0, -2.0), "CFS30B");
      addProjection(ChargedFinalState(-3.5, -2.5), "CFS45B");
      addProjection(ChargedFinalState(-4.0, -3.0), "CFS40B");
            

      // Histogram booking, we have sqrt(s) = 200, 546 and 900 GeV
      _hist_correl_10_200 = bookHistogram1D(1, 1, 1);
      _hist_correl_10_546 = bookHistogram1D(1, 1, 2);
      _hist_correl_10_900 = bookHistogram1D(1, 1, 3);
      
      _hist_correl_15_200 = bookHistogram1D(2, 1, 1);
      _hist_correl_15_546 = bookHistogram1D(2, 1, 2);
      _hist_correl_15_900 = bookHistogram1D(2, 1, 3);
      
      _hist_correl_20_200 = bookHistogram1D(3, 1, 1);
      _hist_correl_20_546 = bookHistogram1D(3, 1, 2);
      _hist_correl_20_900 = bookHistogram1D(3, 1, 3);
      
      _hist_correl_25_200 = bookHistogram1D(4, 1, 1);
      _hist_correl_25_546 = bookHistogram1D(4, 1, 2);
      _hist_correl_25_900 = bookHistogram1D(4, 1, 3);
      
      _hist_correl_30_200 = bookHistogram1D(5, 1, 1);
      _hist_correl_30_546 = bookHistogram1D(5, 1, 2);
      _hist_correl_30_900 = bookHistogram1D(5, 1, 3);
      
      _hist_correl_35_200 = bookHistogram1D(6, 1, 1);
      _hist_correl_35_546 = bookHistogram1D(6, 1, 2);
      _hist_correl_35_900 = bookHistogram1D(6, 1, 3);
      
      _hist_correl_40_200 = bookHistogram1D(7, 1, 1);
      _hist_correl_40_546 = bookHistogram1D(7, 1, 2);
      _hist_correl_40_900 = bookHistogram1D(7, 1, 3);
      
      _hist_correl_asym_15_200 = bookHistogram1D(8, 1, 1);
      _hist_correl_asym_15_546 = bookHistogram1D(8, 1, 2);
      _hist_correl_asym_15_900 = bookHistogram1D(8, 1, 3);
      
      _hist_correl_asym_20_200 = bookHistogram1D(9, 1, 1);
      _hist_correl_asym_20_546 = bookHistogram1D(9, 1, 2);
      _hist_correl_asym_20_900 = bookHistogram1D(9, 1, 3);
      
      _hist_correl_asym_25_200 = bookHistogram1D(10, 1, 1);
      _hist_correl_asym_25_546 = bookHistogram1D(10, 1, 2);
      _hist_correl_asym_25_900 = bookHistogram1D(10, 1, 3);
      
      _hist_correl_asym_30_200 = bookHistogram1D(11, 1, 1);
      _hist_correl_asym_30_546 = bookHistogram1D(11, 1, 2);
      _hist_correl_asym_30_900 = bookHistogram1D(11, 1, 3);
    }
    
    

    void analyze(const Event& event) {
      // Trigger
      const bool trigger = applyProjection<TriggerUA5>(event, "Trigger").nsdDecision();
      if (!trigger) vetoEvent;

      const double sqrtS = applyProjection<Beam>(event, "Beams").sqrtS();
      const double weight = event.weight();
            
      // Count forward/backward rates
      n_10f += applyProjection<ChargedFinalState>(event, "CFS10F").size();
      n_15f += applyProjection<ChargedFinalState>(event, "CFS15F").size();
      n_20f += applyProjection<ChargedFinalState>(event, "CFS20F").size();
      n_25f += applyProjection<ChargedFinalState>(event, "CFS25F").size();
      n_30f += applyProjection<ChargedFinalState>(event, "CFS30F").size();
      n_35f += applyProjection<ChargedFinalState>(event, "CFS35F").size();
      n_40f += applyProjection<ChargedFinalState>(event, "CFS40F").size();
      //
      n_10b += applyProjection<ChargedFinalState>(event, "CFS10B").size();
      n_15b += applyProjection<ChargedFinalState>(event, "CFS15B").size();
      n_20b += applyProjection<ChargedFinalState>(event, "CFS20B").size();
      n_25b += applyProjection<ChargedFinalState>(event, "CFS25B").size();
      n_30b += applyProjection<ChargedFinalState>(event, "CFS30B").size();
      n_35b += applyProjection<ChargedFinalState>(event, "CFS35B").size();
      n_40b += applyProjection<ChargedFinalState>(event, "CFS40B").size();
      //
      n_05 += applyProjection<ChargedFinalState>(event, "CFS05").size();
      
      // Dummy fill
      if (fuzzyEquals(sqrtS, 200.0, 1E-4)) {
        _hist_correl_10_200->fill(_hist_correl_10_200->binMean(0), weight);
        _hist_correl_15_200->fill(_hist_correl_15_200->binMean(0), weight);
        _hist_correl_20_200->fill(_hist_correl_20_200->binMean(0), weight);
        _hist_correl_25_200->fill(_hist_correl_25_200->binMean(0), weight);
        _hist_correl_30_200->fill(_hist_correl_30_200->binMean(0), weight);
        _hist_correl_35_200->fill(_hist_correl_35_200->binMean(0), weight);
        _hist_correl_40_200->fill(_hist_correl_40_200->binMean(0), weight);
        _hist_correl_asym_15_200->fill(_hist_correl_asym_15_200->binMean(0), weight);
        _hist_correl_asym_20_200->fill(_hist_correl_asym_20_200->binMean(0), weight);
        _hist_correl_asym_25_200->fill(_hist_correl_asym_25_200->binMean(0), weight);
        _hist_correl_asym_30_200->fill(_hist_correl_asym_30_200->binMean(0), weight);
      }
      
      else if (fuzzyEquals(sqrtS, 546.0, 1E-4)) {
        _hist_correl_10_546->fill(_hist_correl_10_546->binMean(0), weight);
        _hist_correl_15_546->fill(_hist_correl_15_546->binMean(0), weight);
        _hist_correl_20_546->fill(_hist_correl_20_546->binMean(0), weight);
        _hist_correl_25_546->fill(_hist_correl_25_546->binMean(0), weight);
        _hist_correl_30_546->fill(_hist_correl_30_546->binMean(0), weight);
        _hist_correl_35_546->fill(_hist_correl_35_546->binMean(0), weight);
        _hist_correl_40_546->fill(_hist_correl_40_546->binMean(0), weight);
        _hist_correl_asym_15_546->fill(_hist_correl_asym_15_546->binMean(0), weight);
        _hist_correl_asym_20_546->fill(_hist_correl_asym_20_546->binMean(0), weight);
        _hist_correl_asym_25_546->fill(_hist_correl_asym_25_546->binMean(0), weight);
        _hist_correl_asym_30_546->fill(_hist_correl_asym_30_546->binMean(0), weight);
      }
      
      else if (fuzzyEquals(sqrtS, 900.0, 1E-4)) {
        _hist_correl_10_900->fill(_hist_correl_10_900->binMean(0), weight);
        _hist_correl_15_900->fill(_hist_correl_15_900->binMean(0), weight);
        _hist_correl_20_900->fill(_hist_correl_20_900->binMean(0), weight);
        _hist_correl_25_900->fill(_hist_correl_25_900->binMean(0), weight);
        _hist_correl_30_900->fill(_hist_correl_30_900->binMean(0), weight);
        _hist_correl_35_900->fill(_hist_correl_35_900->binMean(0), weight);
        _hist_correl_40_900->fill(_hist_correl_40_900->binMean(0), weight);
        _hist_correl_asym_15_900->fill(_hist_correl_asym_15_900->binMean(0), weight);
        _hist_correl_asym_20_900->fill(_hist_correl_asym_20_900->binMean(0), weight);
        _hist_correl_asym_25_900->fill(_hist_correl_asym_25_900->binMean(0), weight);
        _hist_correl_asym_30_900->fill(_hist_correl_asym_30_900->binMean(0), weight);
      }
    }
    
    
    void finalize() {
      // Get the correlation coefficients
      //
      // Symmetric eta intervals first
      double correlation_cfs10 = correlation(n_10f, n_10b);
      double correlation_cfs15 = correlation(n_15f, n_15b);
      double correlation_cfs20 = correlation(n_20f, n_20b);
      double correlation_cfs25 = correlation(n_25f, n_25b);
      double correlation_cfs30 = correlation(n_30f, n_30b);
      double correlation_cfs35 = correlation(n_35f, n_35b);
      double correlation_cfs40 = correlation(n_40f, n_40b);

      // Assymetric eta intervals
      //  1.5 ... 2.5 & -1.5 ... -0.5
      double correlation_as_cfs15 = correlation(n_25f, n_15b);
      //  2.0 ... 3.0 & -1.0 ...  0.0
      double correlation_as_cfs20 = correlation(n_30f, n_10b);
      //  2.5 ... 3.5 & -0.5 ...  0.5
      double correlation_as_cfs25 = correlation(n_35f, n_05);
      //  3.0 ... 4.0 &  0.0 ...  1.0
      double correlation_as_cfs30 = correlation(n_40f, n_10f);

      normalize(_hist_correl_10_200, correlation_cfs10);
      normalize(_hist_correl_10_546, correlation_cfs10);
      normalize(_hist_correl_10_900, correlation_cfs10);
      
      normalize(_hist_correl_15_200, correlation_cfs15);
      normalize(_hist_correl_15_546, correlation_cfs15);
      normalize(_hist_correl_15_900, correlation_cfs15);
      
      normalize(_hist_correl_20_200, correlation_cfs20);
      normalize(_hist_correl_20_546, correlation_cfs20);
      normalize(_hist_correl_20_900, correlation_cfs20);
      
      normalize(_hist_correl_25_200, correlation_cfs25);
      normalize(_hist_correl_25_546, correlation_cfs25);
      normalize(_hist_correl_25_900, correlation_cfs25);
      
      normalize(_hist_correl_30_200, correlation_cfs30);
      normalize(_hist_correl_30_546, correlation_cfs30);
      normalize(_hist_correl_30_900, correlation_cfs30);
      
      normalize(_hist_correl_35_200, correlation_cfs35);
      normalize(_hist_correl_35_546, correlation_cfs35);
      normalize(_hist_correl_35_900, correlation_cfs35);
      
      normalize(_hist_correl_40_200, correlation_cfs40);
      normalize(_hist_correl_40_546, correlation_cfs40);
      normalize(_hist_correl_40_900, correlation_cfs40);
      
      normalize(_hist_correl_asym_15_200, correlation_as_cfs15);
      normalize(_hist_correl_asym_15_546, correlation_as_cfs15);
      normalize(_hist_correl_asym_15_900, correlation_as_cfs15);
      
      normalize(_hist_correl_asym_20_200, correlation_as_cfs20);
      normalize(_hist_correl_asym_20_546, correlation_as_cfs20);
      normalize(_hist_correl_asym_20_900, correlation_as_cfs20);
      
      normalize(_hist_correl_asym_25_200, correlation_as_cfs25);
      normalize(_hist_correl_asym_25_546, correlation_as_cfs25);
      normalize(_hist_correl_asym_25_900, correlation_as_cfs25);
      
      normalize(_hist_correl_asym_30_200, correlation_as_cfs30);
      normalize(_hist_correl_asym_30_546, correlation_as_cfs30);
      normalize(_hist_correl_asym_30_900, correlation_as_cfs30);      
    }
    
    //@}
    
    
  private:
    
    /// @name Vectors for storing the number of particles in the different eta intervals per event.
    /// @todo Is there a better way?
    //@{
    
    std::vector<int> n_10f;
    std::vector<int> n_15f;
    std::vector<int> n_20f;
    std::vector<int> n_25f;
    std::vector<int> n_30f;
    std::vector<int> n_35f;
    std::vector<int> n_40f;
                           
    std::vector<int> n_10b;
    std::vector<int> n_15b;
    std::vector<int> n_20b;
    std::vector<int> n_25b;
    std::vector<int> n_30b;
    std::vector<int> n_35b;
    std::vector<int> n_40b;
   
    std::vector<int> n_05;

    //@}


    /// @name Histograms
    //@{

    // Symmetric eta intervals
    AIDA::IHistogram1D *_hist_correl_10_200;
    AIDA::IHistogram1D *_hist_correl_10_546;
    AIDA::IHistogram1D *_hist_correl_10_900;

    AIDA::IHistogram1D *_hist_correl_15_200;
    AIDA::IHistogram1D *_hist_correl_15_546;
    AIDA::IHistogram1D *_hist_correl_15_900;

    AIDA::IHistogram1D *_hist_correl_20_200;
    AIDA::IHistogram1D *_hist_correl_20_546;
    AIDA::IHistogram1D *_hist_correl_20_900;
    
    AIDA::IHistogram1D *_hist_correl_25_200;
    AIDA::IHistogram1D *_hist_correl_25_546;
    AIDA::IHistogram1D *_hist_correl_25_900;
    
    AIDA::IHistogram1D *_hist_correl_30_200;
    AIDA::IHistogram1D *_hist_correl_30_546;
    AIDA::IHistogram1D *_hist_correl_30_900;

    AIDA::IHistogram1D *_hist_correl_35_200;
    AIDA::IHistogram1D *_hist_correl_35_900;
    AIDA::IHistogram1D *_hist_correl_35_546;
    
    AIDA::IHistogram1D *_hist_correl_40_200;
    AIDA::IHistogram1D *_hist_correl_40_546;
    AIDA::IHistogram1D *_hist_correl_40_900;
    
    // For asymmetric eta intervals
    AIDA::IHistogram1D *_hist_correl_asym_15_200;
    AIDA::IHistogram1D *_hist_correl_asym_15_546;
    AIDA::IHistogram1D *_hist_correl_asym_15_900;
                                      
    AIDA::IHistogram1D *_hist_correl_asym_20_200;
    AIDA::IHistogram1D *_hist_correl_asym_20_546;
    AIDA::IHistogram1D *_hist_correl_asym_20_900;
                                      
    AIDA::IHistogram1D *_hist_correl_asym_25_200;
    AIDA::IHistogram1D *_hist_correl_asym_25_546;
    AIDA::IHistogram1D *_hist_correl_asym_25_900;
                                      
    AIDA::IHistogram1D *_hist_correl_asym_30_200;
    AIDA::IHistogram1D *_hist_correl_asym_30_546;
    AIDA::IHistogram1D *_hist_correl_asym_30_900;
    //@}

  };


  
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<UA5_1988_S1867512> plugin_UA5_1988_S1867512;

}
