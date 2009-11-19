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
    
    
    inline double cov_w_mean(int m, double m_mean, int n, double n_mean) {
      return (m - m_mean)*(n - n_mean);
    }
  
  
    /// Calculate the correlation strength between two samples
    inline double c_str(int m, double m_mean, int n, double n_mean) {
      const double cov  = cov_w_mean(m, m_mean, n, n_mean);
      const double var1 = cov_w_mean(m, m_mean, m, m_mean);
      const double var2 = cov_w_mean(n, n_mean, n, n_mean);
      const double correlation = cov/sqrt(var1*var2);
      const double corr_strength = correlation*sqrt(var2/var1);
      return corr_strength;
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
      addProjection(ChargedFinalState(2.5, 3.5), "CFS35F");
      addProjection(ChargedFinalState(3.0, 4.0), "CFS40F");
      
      // Backward eta intervals
      addProjection(ChargedFinalState(-1.0,  0.0), "CFS10B");
      addProjection(ChargedFinalState(-1.5, -0.5), "CFS15B");
      addProjection(ChargedFinalState(-2.0, -1.0), "CFS20B");
      addProjection(ChargedFinalState(-2.5, -1.5), "CFS25B");
      addProjection(ChargedFinalState(-3.0, -2.0), "CFS30B");
      addProjection(ChargedFinalState(-3.5, -2.5), "CFS35B");
      addProjection(ChargedFinalState(-4.0, -3.0), "CFS40B");
            
      // Histogram booking, we have sqrt(s) = 200, 546 and 900 GeV
      _hist_correl_200 = bookProfile1D(2, 1, 1);
      _hist_correl_546 = bookProfile1D(2, 1, 2);
      _hist_correl_900 = bookProfile1D(2, 1, 3);
      
      _hist_correl_asym_200 = bookProfile1D(3, 1, 1);
      _hist_correl_asym_546 = bookProfile1D(3, 1, 2);
      _hist_correl_asym_900 = bookProfile1D(3, 1, 3);
    }
    
    

    void analyze(const Event& event) {
      sqrtS = applyProjection<Beam>(event, "Beams").sqrtS();
      
      // Trigger
      const bool trigger = applyProjection<TriggerUA5>(event, "Trigger").nsdDecision();
      if (!trigger) vetoEvent;

            
      // Count forward/backward particles
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
      
    }
    
    
    void finalize() {
      // Calculate mean number of particles in eta intervals
      double mean_n_10f = mean(n_10f);
      double mean_n_15f = mean(n_15f);
      double mean_n_20f = mean(n_20f);
      double mean_n_25f = mean(n_25f);
      double mean_n_30f = mean(n_30f);
      double mean_n_35f = mean(n_35f);
      double mean_n_40f = mean(n_40f);
                 
      double mean_n_10b = mean(n_10b);
      double mean_n_15b = mean(n_15b);
      double mean_n_20b = mean(n_20b);
      double mean_n_25b = mean(n_25b);
      double mean_n_30b = mean(n_30b);
      double mean_n_35b = mean(n_35b);
      double mean_n_40b = mean(n_40b);
                 
      double mean_n_05  = mean(n_05) ;


      // Fill histos
      if (fuzzyEquals(sqrtS, 200.0, 1E-4)) {
        for (size_t i = 0; i < n_10f.size(); i++) {
          // Fill gap size histo (Fig 14), iterate over central gap size
          _hist_correl_200->fill(0.0, c_str(n_10f[i], mean_n_10f, n_10b[i], mean_n_10b));
          _hist_correl_200->fill(1.0, c_str(n_15f[i], mean_n_15f, n_15b[i], mean_n_15b));
          _hist_correl_200->fill(2.0, c_str(n_20f[i], mean_n_20f, n_20b[i], mean_n_20b));
          _hist_correl_200->fill(3.0, c_str(n_25f[i], mean_n_25f, n_25b[i], mean_n_25b));
          _hist_correl_200->fill(4.0, c_str(n_30f[i], mean_n_30f, n_30b[i], mean_n_30b));
          _hist_correl_200->fill(5.0, c_str(n_35f[i], mean_n_35f, n_35b[i], mean_n_35b));
          _hist_correl_200->fill(6.0, c_str(n_40f[i], mean_n_40f, n_40b[i], mean_n_40b));

          // Fill gap-center histos (Fig 15), iterate over gap centers
          //
          // The first bin contains all the c_str strengths of
          // the gap size histo above
          _hist_correl_asym_200->fill(0.0, c_str(n_10f[i], mean_n_10f, n_10b[i], mean_n_10b));
          _hist_correl_asym_200->fill(0.0, c_str(n_15f[i], mean_n_15f, n_15b[i], mean_n_15b));
          _hist_correl_asym_200->fill(0.0, c_str(n_20f[i], mean_n_20f, n_20b[i], mean_n_20b));
          _hist_correl_asym_200->fill(0.0, c_str(n_25f[i], mean_n_25f, n_25b[i], mean_n_25b));
          _hist_correl_asym_200->fill(0.0, c_str(n_30f[i], mean_n_30f, n_30b[i], mean_n_30b));
          _hist_correl_asym_200->fill(0.0, c_str(n_35f[i], mean_n_35f, n_35b[i], mean_n_35b));
          _hist_correl_asym_200->fill(0.0, c_str(n_40f[i], mean_n_40f, n_40b[i], mean_n_40b));
          // Fill in c_str strength for assymetric intervals
          _hist_correl_asym_200->fill(0.5, c_str(n_25f[i], mean_n_25f, n_15b[i], mean_n_15b));
          _hist_correl_asym_200->fill(1.0, c_str(n_30f[i], mean_n_30f, n_10b[i], mean_n_10b));
          _hist_correl_asym_200->fill(1.5, c_str(n_35f[i], mean_n_35f, n_05[i] , mean_n_05 ));
          _hist_correl_asym_200->fill(2.0, c_str(n_40f[i], mean_n_40f, n_10f[i], mean_n_10f));
        }
      }
      
      else if (fuzzyEquals(sqrtS, 546.0, 1E-4)) {
        for (size_t i = 0; i < n_10f.size(); i++) {
          _hist_correl_546->fill(0.0, c_str(n_10f[i], mean_n_10f, n_10b[i], mean_n_10b));
          _hist_correl_546->fill(1.0, c_str(n_15f[i], mean_n_15f, n_15b[i], mean_n_15b));
          _hist_correl_546->fill(2.0, c_str(n_20f[i], mean_n_20f, n_20b[i], mean_n_20b));
          _hist_correl_546->fill(3.0, c_str(n_25f[i], mean_n_25f, n_25b[i], mean_n_25b));
          _hist_correl_546->fill(4.0, c_str(n_30f[i], mean_n_30f, n_30b[i], mean_n_30b));
          _hist_correl_546->fill(5.0, c_str(n_35f[i], mean_n_35f, n_35b[i], mean_n_35b));
          _hist_correl_546->fill(6.0, c_str(n_40f[i], mean_n_40f, n_40b[i], mean_n_40b));

          _hist_correl_asym_546->fill(0.0, c_str(n_10f[i], mean_n_10f, n_10b[i], mean_n_10b));
          _hist_correl_asym_546->fill(0.0, c_str(n_15f[i], mean_n_15f, n_15b[i], mean_n_15b));
          _hist_correl_asym_546->fill(0.0, c_str(n_20f[i], mean_n_20f, n_20b[i], mean_n_20b));
          _hist_correl_asym_546->fill(0.0, c_str(n_25f[i], mean_n_25f, n_25b[i], mean_n_25b));
          _hist_correl_asym_546->fill(0.0, c_str(n_30f[i], mean_n_30f, n_30b[i], mean_n_30b));
          _hist_correl_asym_546->fill(0.0, c_str(n_35f[i], mean_n_35f, n_35b[i], mean_n_35b));
          _hist_correl_asym_546->fill(0.0, c_str(n_40f[i], mean_n_40f, n_40b[i], mean_n_40b));
          
          _hist_correl_asym_546->fill(0.5, c_str(n_25f[i], mean_n_25f, n_15b[i], mean_n_15b));
          _hist_correl_asym_546->fill(1.0, c_str(n_30f[i], mean_n_30f, n_10b[i], mean_n_10b));
          _hist_correl_asym_546->fill(1.5, c_str(n_35f[i], mean_n_35f, n_05[i] , mean_n_05 ));
          _hist_correl_asym_546->fill(2.0, c_str(n_40f[i], mean_n_40f, n_10f[i], mean_n_10f));
        }
      }
      
      else if (fuzzyEquals(sqrtS, 900.0, 1E-4)) {
        for (size_t i = 0; i < n_10f.size(); i++) {
          _hist_correl_900->fill(0.0, c_str(n_10f[i], mean_n_10f, n_10b[i], mean_n_10b));
          _hist_correl_900->fill(1.0, c_str(n_15f[i], mean_n_15f, n_15b[i], mean_n_15b));
          _hist_correl_900->fill(2.0, c_str(n_20f[i], mean_n_20f, n_20b[i], mean_n_20b));
          _hist_correl_900->fill(3.0, c_str(n_25f[i], mean_n_25f, n_25b[i], mean_n_25b));
          _hist_correl_900->fill(4.0, c_str(n_30f[i], mean_n_30f, n_30b[i], mean_n_30b));
          _hist_correl_900->fill(5.0, c_str(n_35f[i], mean_n_35f, n_35b[i], mean_n_35b));
          _hist_correl_900->fill(6.0, c_str(n_40f[i], mean_n_40f, n_40b[i], mean_n_40b));

          _hist_correl_asym_900->fill(0.0, c_str(n_10f[i], mean_n_10f, n_10b[i], mean_n_10b));
          _hist_correl_asym_900->fill(0.0, c_str(n_15f[i], mean_n_15f, n_15b[i], mean_n_15b));
          _hist_correl_asym_900->fill(0.0, c_str(n_20f[i], mean_n_20f, n_20b[i], mean_n_20b));
          _hist_correl_asym_900->fill(0.0, c_str(n_25f[i], mean_n_25f, n_25b[i], mean_n_25b));
          _hist_correl_asym_900->fill(0.0, c_str(n_30f[i], mean_n_30f, n_30b[i], mean_n_30b));
          _hist_correl_asym_900->fill(0.0, c_str(n_35f[i], mean_n_35f, n_35b[i], mean_n_35b));
          _hist_correl_asym_900->fill(0.0, c_str(n_40f[i], mean_n_40f, n_40b[i], mean_n_40b));
          
          _hist_correl_asym_900->fill(0.5, c_str(n_25f[i], mean_n_25f, n_15b[i], mean_n_15b));
          _hist_correl_asym_900->fill(1.0, c_str(n_30f[i], mean_n_30f, n_10b[i], mean_n_10b));
          _hist_correl_asym_900->fill(1.5, c_str(n_35f[i], mean_n_35f, n_05[i] , mean_n_05 ));
          _hist_correl_asym_900->fill(2.0, c_str(n_40f[i], mean_n_40f, n_10f[i], mean_n_10f));
        }
      }

      
    }
    
    //@}
    
    
  private:
   
    // CoM energy
    double sqrtS;

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
    AIDA::IProfile1D *_hist_correl_200;
    AIDA::IProfile1D *_hist_correl_546;
    AIDA::IProfile1D *_hist_correl_900;
    
    // For asymmetric eta intervals
    AIDA::IProfile1D *_hist_correl_asym_200;
    AIDA::IProfile1D *_hist_correl_asym_546;
    AIDA::IProfile1D *_hist_correl_asym_900;
    //@}

  };


  
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<UA5_1988_S1867512> plugin_UA5_1988_S1867512;

}
