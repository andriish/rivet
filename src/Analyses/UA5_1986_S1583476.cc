// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {

  class UA5_1986_S1583476 : public Analysis {
  public:

    /// Constructor
    UA5_1986_S1583476() 
      : Analysis("UA5_1986_S1583476") 
    {
      setBeams(PROTON, ANTIPROTON);
      addProjection(Beam(), "Beams");
      
      // All charged final state particles, needed for triggers
      const ChargedFinalState cfs;
      addProjection(cfs, "CFSAll");
      // Charged particles in |eta| < 5.0 
      const ChargedFinalState cfs50(-5.0, 5.0);
      addProjection(cfs50, "CFS50");
    }
    


    /// @name Analysis methods
    //@{
    
    void init() {
      _hist_eta_nsd_200      = bookHistogram1D(1,1,1);
      _hist_eta_inelastic_200  = bookHistogram1D(1,1,2);
      _hist_eta_nsd_900      = bookHistogram1D(1,1,3);
      _hist_eta_inelastic_900  = bookHistogram1D(1,1,4);
      
      _hist_eta_nsd_n_2_200  = bookHistogram1D(2,1,1);
      _hist_eta_nsd_n_12_200 = bookHistogram1D(2,1,2);
      _hist_eta_nsd_n_22_200 = bookHistogram1D(2,1,3);
      _hist_eta_nsd_n_32_200 = bookHistogram1D(2,1,4);
      _hist_eta_nsd_n_42_200 = bookHistogram1D(2,1,5);
      _hist_eta_nsd_n_52_200 = bookHistogram1D(2,1,6);
      
      _hist_eta_nsd_n_2_900  = bookHistogram1D(3,1,1);
      _hist_eta_nsd_n_12_900 = bookHistogram1D(3,1,2);
      _hist_eta_nsd_n_22_900 = bookHistogram1D(3,1,3);
      _hist_eta_nsd_n_32_900 = bookHistogram1D(3,1,4);
      _hist_eta_nsd_n_42_900 = bookHistogram1D(3,1,5);
      _hist_eta_nsd_n_52_900 = bookHistogram1D(3,1,6);
      _hist_eta_nsd_n_62_900 = bookHistogram1D(3,1,7);
      _hist_eta_nsd_n_72_900 = bookHistogram1D(3,1,8);
      _hist_eta_nsd_n_82_900 = bookHistogram1D(3,1,9);
    }
    
    
    void analyze(const Event& event) {
      const double sqrtS = applyProjection<Beam>(event, "Beams").sqrtS();
      const double weight = event.weight();
      
      // Trigger requirements from the hodoscopes (1 arm (elastic) and 2 arms (NSD))
      int n_trig_1(0), n_trig_2(0);
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFSAll");
      foreach (const Particle& p, cfs.particles()) {
        double eta = p.momentum().pseudorapidity();
        if (inRange(eta, -5.6, -2.0)) n_trig_1 += 1;
        else if (inRange(eta, 2.0, 5.6)) n_trig_2 += 1;
      }
      getLog() << Log::DEBUG << "Trigger -: " << n_trig_1 << ", Trigger +: " << n_trig_2 << endl;
      
      // Check if we have a coincidence hit in hodoscopes == NSD
      if (n_trig_1 == 0 && n_trig_2 == 0) vetoEvent
      // Require at least one hit in trigger hodoscopes
      bool isNSD = true;
      if (n_trig_1 == 0 || n_trig_2 == 0) isNSD = false;
      
      // Declare final state for |eta| < 5.0
      const ChargedFinalState& cfs50 = applyProjection<ChargedFinalState>(event, "CFS50");
      int numP = cfs50.particles().size();


      // Iterate over particles in |eta| < 5.0 and fill histos with |eta| 
      foreach (const Particle& p, cfs.particles()) {
        double eta = fabs(p.momentum().pseudorapidity());
        
        // Fill 200 GeV histos
        if (fuzzyEquals(sqrtS, 200.0, 1E-4)) {
          // Fill histos that don't require a certain multiplicity
          _hist_eta_inelastic_200->fill(eta, weight);
          if ( isNSD ) {
            // Fill histos that require a certain multiplicity
            _hist_eta_nsd_200->fill(eta, weight);
            if ( ( 2 <= numP ) && ( numP <= 10 ) ) _hist_eta_nsd_n_2_200->fill(eta, weight);
            else if ( ( 12 <= numP ) && ( numP <= 20 ) ) _hist_eta_nsd_n_12_200->fill(eta, weight);
            else if ( ( 22 <= numP ) && ( numP <= 30 ) ) _hist_eta_nsd_n_22_200->fill(eta, weight);
            else if ( ( 32 <= numP ) && ( numP <= 40 ) ) _hist_eta_nsd_n_32_200->fill(eta, weight);
            else if ( ( 42 <= numP ) && ( numP <= 50 ) ) _hist_eta_nsd_n_42_200->fill(eta, weight);
            else if ( numP >= 52 ) _hist_eta_nsd_n_52_200->fill(eta, weight);
          }
        }
        
        // Fill 900 GeV histos
        else if (fuzzyEquals(sqrtS, 900.0, 1E-4)) {
          // Fill histos that don't require a certain multiplicity
          _hist_eta_inelastic_900->fill(eta, weight);
          if ( isNSD ) {
            // Fill histos that require a certain multiplicity
            _hist_eta_nsd_900->fill(eta, weight);
            if ( ( 2 <= numP ) && ( numP <= 10 ) ) _hist_eta_nsd_n_2_900->fill(eta, weight);
            else if ( ( 12 <= numP ) && ( numP <= 20 ) ) _hist_eta_nsd_n_12_900->fill(eta, weight);
            else if ( ( 22 <= numP ) && ( numP <= 30 ) ) _hist_eta_nsd_n_22_900->fill(eta, weight);
            else if ( ( 32 <= numP ) && ( numP <= 40 ) ) _hist_eta_nsd_n_32_900->fill(eta, weight);
            else if ( ( 42 <= numP ) && ( numP <= 50 ) ) _hist_eta_nsd_n_42_900->fill(eta, weight);
            else if ( ( 52 <= numP ) && ( numP <= 60 ) ) _hist_eta_nsd_n_52_900->fill(eta, weight);
            else if ( ( 62 <= numP ) && ( numP <= 70 ) ) _hist_eta_nsd_n_62_900->fill(eta, weight);
            else if ( ( 72 <= numP ) && ( numP <= 80 ) ) _hist_eta_nsd_n_72_900->fill(eta, weight);
            else if ( numP >= 82 ) _hist_eta_nsd_n_82_900->fill(eta, weight);
          }
        } 
      }        
  }
  
  
    void finalize() {
      // Scale histos to the area of the corresponding reference histos
      normalize(_hist_eta_nsd_200, 10.2225);    
      normalize(_hist_eta_inelastic_200, 9.255);
      normalize(_hist_eta_nsd_900, 15.285);
      normalize(_hist_eta_inelastic_900, 13.9725);
      
      normalize(_hist_eta_nsd_n_2_200, 3.285);
      normalize(_hist_eta_nsd_n_12_200, 7.34);
      normalize(_hist_eta_nsd_n_22_200, 12.02);
      normalize(_hist_eta_nsd_n_32_200, 17.2);
      normalize(_hist_eta_nsd_n_42_200, 21.99);
      normalize(_hist_eta_nsd_n_52_200, 27.8);
      
      normalize(_hist_eta_nsd_n_2_900, 2.7);
      normalize(_hist_eta_nsd_n_12_900, 6.425);
      normalize(_hist_eta_nsd_n_22_900, 10.54);
      normalize(_hist_eta_nsd_n_32_900, 15.225);
      normalize(_hist_eta_nsd_n_42_900, 19.885);
      normalize(_hist_eta_nsd_n_52_900, 25.13);
      normalize(_hist_eta_nsd_n_62_900, 29.235);
      normalize(_hist_eta_nsd_n_72_900, 33.81);
      normalize(_hist_eta_nsd_n_82_900, 41.75);      
    }
    

  private:

    /// @name Histograms
    //@{
    // Histos of Figure 1 (HepData Table 1)
    AIDA::IHistogram1D *_hist_eta_nsd_200;
    AIDA::IHistogram1D *_hist_eta_inelastic_200;
    AIDA::IHistogram1D *_hist_eta_nsd_900;
    AIDA::IHistogram1D *_hist_eta_inelastic_900;

    // Histos of Figure 3a (HepData Table 2)
    AIDA::IHistogram1D *_hist_eta_nsd_n_2_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_12_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_22_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_32_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_42_200;
    AIDA::IHistogram1D *_hist_eta_nsd_n_52_200;

    // Histos of Figure 3b (HepData Table 3)
    AIDA::IHistogram1D *_hist_eta_nsd_n_2_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_12_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_22_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_32_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_42_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_52_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_62_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_72_900;
    AIDA::IHistogram1D *_hist_eta_nsd_n_82_900;
    //@}

  };
 

  
  // This global object acts as a hook for the plugin system
  AnalysisBuilder<UA5_1986_S1583476> plugin_UA5_1986_S1583476;
  
}
