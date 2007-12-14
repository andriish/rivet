// -*- C++ -*-
#ifndef RIVET_D0_2001_S4674421_HH
#define RIVET_D0_2001_S4674421_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/D0ILConeJets.hh"
#include "Rivet/Projections/PVertex.hh"
#include "Rivet/Projections/TotalVisibleMomentum.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {  

  /// Implementation of D0 Run I, differential W/Z boson cross 
  /// section paper hep-ex/0107012
  /// this is the version from Oleg just with adjustments
  /// like e.g. replacing Root histograms by our LWH implementation
  /// It is NOT ready for release!!!
  class D0_2001_S4674421 : public Analysis {

  public:

    /// Default constructor.
    D0_2001_S4674421()
      // NB. eta in [-3,3] cut specified via FinalState constructor
      : _fsproj(-3.0, 3.0) 
    { 

      setBeams(PROTON, ANTIPROTON);
      
      addProjection(_fsproj);

      /// @todo Fix this!      
      _bins_pt_w.resize(23);
      _bins_pt_w[0] = 0.;
      _bins_pt_w[1] = 2.;
      _bins_pt_w[2] = 4.;
      _bins_pt_w[3] = 6.;
      _bins_pt_w[4] = 8.;
      _bins_pt_w[5] = 10.;
      _bins_pt_w[6] = 12.;
      _bins_pt_w[7] = 14.;
      _bins_pt_w[8] = 16.;
      _bins_pt_w[9] = 18.;
      _bins_pt_w[10] = 20.;
      _bins_pt_w[11] = 25.;
      _bins_pt_w[12] = 30.;
      _bins_pt_w[13] = 35.;
      _bins_pt_w[14] = 40.;
      _bins_pt_w[15] = 50.;
      _bins_pt_w[16] = 60.;
      _bins_pt_w[17] = 70.;
      _bins_pt_w[18] = 80.;
      _bins_pt_w[19] = 100.;
      _bins_pt_w[20] = 120.;
      _bins_pt_w[21] = 160.;
      _bins_pt_w[22] = 200.;
      
    }    


    /// Factory method
    static Analysis* create() { return new D0_2001_S4674421(); }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "4674421";
    }
    /// Get a description of the analysis.
    // string getDescription() const {
    //   return "";
    // }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "D0";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2001";
    }
    //@}


    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// The final state projector used by this analysis.
    FinalState _fsproj;


    /// pT bins to be distinguished during analysis
    /// @todo Remove in favour of autobinning.
    vector<double> _bins_pt_w;


    D0_2001_S4674421& operator=(const D0_2001_S4674421& x);


    //@{
    /// Histograms
    AIDA::IHistogram1D* _h_pt_e;
    AIDA::IHistogram1D* _h_pt_miss;
    AIDA::IHistogram1D* _h_pt_w;
    AIDA::IHistogram1D* _h_pt_w_true;
    AIDA::IHistogram1D* _h_w_rec_eff_factor;
    //@}    


  };

}

#endif
