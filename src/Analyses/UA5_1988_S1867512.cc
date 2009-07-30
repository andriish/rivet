// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Analyses/UA5_1988_S1867512.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/AnalysisLoader.hh"

namespace Rivet {


  UA5_1988_S1867512::UA5_1988_S1867512()
      : Analysis("UA5_1988_S1867512")
  {
    /// @todo Set approriate for your analysis
    setBeams(PROTON, ANTIPROTON);
    addProjection(Beam(), "Beams");
   
    // All charged final state particles, needed for trigger implementation only
    const ChargedFinalState cfs;
    addProjection(cfs,   "CFSAll");
    
    // Symmetric intervals first
    // Maybe its possible to define symmetric eta intervals with gaps
    // Forward eta intervals
    const ChargedFinalState cfs10f(0.0, 1.0);
    const ChargedFinalState cfs15f(0.5, 1.5);
    const ChargedFinalState cfs20f(1.0, 2.0);
    const ChargedFinalState cfs25f(1.5, 2.5);
    const ChargedFinalState cfs30f(2.0, 3.0);
    const ChargedFinalState cfs35f(2.5, 3.5);
    const ChargedFinalState cfs40f(3.0, 4.0);
      
    // Backward eta intervals
    const ChargedFinalState cfs10b(-1.0,  0.0);
    const ChargedFinalState cfs15b(-1.5, -0.5);
    const ChargedFinalState cfs20b(-2.0, -1.0);
    const ChargedFinalState cfs25b(-2.5, -1.5);
    const ChargedFinalState cfs30b(-3.0, -2.0);
    const ChargedFinalState cfs35b(-3.5, -2.5);
    const ChargedFinalState cfs40b(-4.0, -3.0);

    // Symmetric eta interval
    const ChargedFinalState cfs05(-0.5,  0.5);

    addProjection(cfs10f, "CFS10F");
    addProjection(cfs15f, "CFS15F");
    addProjection(cfs20f, "CFS20F");
    addProjection(cfs25f, "CFS25F");
    addProjection(cfs30f, "CFS30F");
    addProjection(cfs35f, "CFS35F");
    addProjection(cfs40f, "CFS40F");

    addProjection(cfs10b, "CFS10B");
    addProjection(cfs15b, "CFS15B");
    addProjection(cfs20b, "CFS20B");
    addProjection(cfs25b, "CFS25B");
    addProjection(cfs30b, "CFS30B");
    addProjection(cfs35b, "CFS35B");
    addProjection(cfs40b, "CFS40B");

    addProjection(cfs05, "CFS05");
  }


  void UA5_1988_S1867512::init() {
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

  double UA5_1988_S1867512::calc_mean(std::vector<int> sample) {
      // A simple function to calculate the mean of a sample
      double mean = 0.0;
      foreach (const int& i, sample) {
          mean += i;
      }
      return mean/sample.size();
  }
  
  double UA5_1988_S1867512::calc_covariance(std::vector<int> sample1, std::vector<int> sample2) {
      // A function to calculate the covariance (variance) of one quality (n_particles)
      // between two samples
      double mean1 = UA5_1988_S1867512::calc_mean(sample1);
      double mean2 = UA5_1988_S1867512::calc_mean(sample2);
      int N = sample1.size();
      double cov = 0.0;
      for (int i = 0; i < N; i++) {
          double cov_i = (sample1[i] - mean1)*(sample2[i] - mean2);
            cov += cov_i;
            }
      if ( N > 1 ) return cov/(N-1);
      else return 0.0;
  }

  double UA5_1988_S1867512::calc_correlation(std::vector<int> sample1, std::vector<int> sample2) {
      // A function to calculate the correlation strength of one quality
      // (n_particles) between two samples
      double cov = UA5_1988_S1867512::calc_covariance(sample1, sample2);
      double var1 = UA5_1988_S1867512::calc_covariance(sample1, sample1);
      double var2 = UA5_1988_S1867512::calc_covariance(sample2, sample2);
      double correlation = cov/sqrt(var1*var2);
      double corr_strength = correlation*sqrt(var2/var1);
      return corr_strength;
  }


  void UA5_1988_S1867512::analyze(const Event& event) {
      Log log = getLog();

      const double sqrtS = applyProjection<Beam>(event, "Beams").sqrtS();
      const double weight = event.weight();

      // Minimum Bias trigger requirements from the hodoscopes
      int n_trig_1 = 0;
      int n_trig_2 = 0;
      
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFSAll");
      foreach (const Particle& p, cfs.particles()) {
           double eta = p.momentum().pseudorapidity();
           if ( ( -5.6 < eta ) && ( eta < -2.0 ) ) n_trig_1++;
           else if ( ( 2.0 < eta ) && ( eta < 5.6 ) ) n_trig_2++;
      }
      
      // Require at least one coincidence hit in trigger hodoscopes
      if ( n_trig_1* n_trig_2 < 1. ) vetoEvent; 
      getLog() << Log::DEBUG << "Trigger 1: " << n_trig_1 << " Trigger 2: " << n_trig_2 << endl;
      
      // Declare final states in several eta regions
      const ChargedFinalState& cfs10f = applyProjection<ChargedFinalState>(event, "CFS10F");
      const ChargedFinalState& cfs15f = applyProjection<ChargedFinalState>(event, "CFS15F");
      const ChargedFinalState& cfs20f = applyProjection<ChargedFinalState>(event, "CFS20F");
      const ChargedFinalState& cfs25f = applyProjection<ChargedFinalState>(event, "CFS25F");
      const ChargedFinalState& cfs30f = applyProjection<ChargedFinalState>(event, "CFS30F");
      const ChargedFinalState& cfs35f = applyProjection<ChargedFinalState>(event, "CFS35F");
      const ChargedFinalState& cfs40f = applyProjection<ChargedFinalState>(event, "CFS40F");
      
      const ChargedFinalState& cfs10b = applyProjection<ChargedFinalState>(event, "CFS10B");
      const ChargedFinalState& cfs15b = applyProjection<ChargedFinalState>(event, "CFS15B");
      const ChargedFinalState& cfs20b = applyProjection<ChargedFinalState>(event, "CFS20B");
      const ChargedFinalState& cfs25b = applyProjection<ChargedFinalState>(event, "CFS25B");
      const ChargedFinalState& cfs30b = applyProjection<ChargedFinalState>(event, "CFS30B");
      const ChargedFinalState& cfs35b = applyProjection<ChargedFinalState>(event, "CFS35B");
      const ChargedFinalState& cfs40b = applyProjection<ChargedFinalState>(event, "CFS40B");

      const ChargedFinalState& cfs05 = applyProjection<ChargedFinalState>(event, "CFS05");
      // Push back the number of particles in the different regions
      n_10f.push_back(cfs10f.particles().size());
      n_15f.push_back(cfs15f.particles().size());
      n_20f.push_back(cfs20f.particles().size());
      n_25f.push_back(cfs25f.particles().size());
      n_30f.push_back(cfs30f.particles().size());
      n_35f.push_back(cfs35f.particles().size());
      n_40f.push_back(cfs40f.particles().size());
                            
      n_10b.push_back(cfs10b.particles().size());
      n_15b.push_back(cfs15b.particles().size());
      n_20b.push_back(cfs20b.particles().size());
      n_25b.push_back(cfs25b.particles().size());
      n_30b.push_back(cfs30b.particles().size());
      n_35b.push_back(cfs35b.particles().size());
      n_40b.push_back(cfs40b.particles().size());

      n_05.push_back(cfs05.particles().size());

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


  void UA5_1988_S1867512::finalize() {
      // Get the correlation coefficients
      //
      // Symmetric eta intervals first
      double correlation_cfs10 = UA5_1988_S1867512::calc_correlation(n_10f, n_10b);
      double correlation_cfs15 = UA5_1988_S1867512::calc_correlation(n_15f, n_15b);
      double correlation_cfs20 = UA5_1988_S1867512::calc_correlation(n_20f, n_20b);
      double correlation_cfs25 = UA5_1988_S1867512::calc_correlation(n_25f, n_25b);
      double correlation_cfs30 = UA5_1988_S1867512::calc_correlation(n_30f, n_30b);
      double correlation_cfs35 = UA5_1988_S1867512::calc_correlation(n_35f, n_35b);
      double correlation_cfs40 = UA5_1988_S1867512::calc_correlation(n_40f, n_40b);

      // Assymetric eta intervals
      //  1.5 ... 2.5 & -1.5 ... -0.5
      double correlation_as_cfs15 = UA5_1988_S1867512::calc_correlation(n_25f, n_15b);
      //  2.0 ... 3.0 & -1.0 ...  0.0
      double correlation_as_cfs20 = UA5_1988_S1867512::calc_correlation(n_30f, n_10b);
      //  2.5 ... 3.5 & -0.5 ...  0.5
      double correlation_as_cfs25 = UA5_1988_S1867512::calc_correlation(n_35f, n_05);
      //  3.0 ... 4.0 &  0.0 ...  1.0
      double correlation_as_cfs30 = UA5_1988_S1867512::calc_correlation(n_40f, n_10f);

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


}
