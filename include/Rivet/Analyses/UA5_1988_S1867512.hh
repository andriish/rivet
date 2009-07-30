// -*- C++ -*-
#ifndef RIVET_UA5_1988_S1867512_HH
#define RIVET_UA5_1988_S1867512_HH

#include "Rivet/Analysis.hh"

namespace Rivet {


  class UA5_1988_S1867512 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    UA5_1988_S1867512();

    /// Factory method
    static Analysis* create() {
      return new UA5_1988_S1867512();
    }
    //@}


  public:

    /// @name Analysis methods
    //@{
     //string name() const {
      //return "UA5_1988_S1867512";
    //}
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    // Vectors for storing the number of particles in the
    // different  eta intervals per event. Memory leak?
    // Is there a better way???
  private:
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

    double calc_mean(std::vector<int> sample);
    double calc_covariance(std::vector<int> sample1, std::vector<int> sample2);
    double calc_correlation(std::vector<int> sample1, std::vector<int> sample2);

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
  

}

#endif

