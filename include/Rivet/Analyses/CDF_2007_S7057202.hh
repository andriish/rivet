// -*- C++ -*-
#ifndef RIVET_CDF_2007_S7057202_HH
#define RIVET_CDF_2007_S7057202_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief CDF Run II inclusive jet cross-section using the kT algorithm.
  /// @author James Monk
  class CDF_2007_S7057202 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2007_S7057202()
      : _minY(0.1), _maxY(0.7), _jetMinPT(54.0*GeV)
   {
      setBeams(PROTON, ANTIPROTON);
      //setSqrtS(1960*GeV);
      const FinalState fs;
      addProjection(FastJets(fs, FastJets::KT, 0.5), "JetsD05");
      addProjection(FastJets(fs, FastJets::KT, 0.7), "JetsD07");
      addProjection(FastJets(fs, FastJets::KT, 1.0), "JetsD10");
      setNeedsCrossSection(true);
    }
    
    /// Factory method
    static Analysis* create() { 
      return new CDF_2007_S7057202(); 
    }
    //@}
    

    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "7057202";
    }
    /// A short description of the analysis.
    string summary() const {
      return "CDF Run II inclusive jet cross-section using the kT algorithm";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "CDF";
    }
    
    /// Collider on which the experiment ran.
    string collider() const {
      return "Tevatron Run 2";
    }
    
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2007";
    }
    //@}

    /// Names & emails of paper/analysis authors.
    vector<string> authors() const {
      vector<string> rtn;
      rtn.push_back( "David Voong" );
      rtn.push_back( "James Monk <jmonk@hep.ucl.ac.uk>" );
      return rtn;
    }
    
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "CDF Run II measurement of inclusive jet cross sections at a p-pbar collision energy of 1.96~TeV"
         << "Jets are reconstructed in the central part of the detector ($|y|<2.1$) using the kT algorithm with an R parameter of 0.7."
         << "The minimum jet pT considered is 54~GeV, with a maximum around 700~GeV."
         << "The inclusive jet pT is plotted in bins of rapidity $|y|<0.1$, $0.1<|y|<0.7$, $0.7<|y|<1.1$, $1.1<|y|<1.6$ and $1.6<|y|<2.1$.";
      return os.str();
    }
    
    /// Information about the events needed as input for this analysis.
    string runInfo() const{
      ostringstream os;
      os << "Standard Tevatron Run II: p-pbar collisions at 1960~GeV."
         << " Jet pT bins from 54~GeV to 700~GeV."  
         << " Jet rapidity $< |2.1|$.";
      return os.str();
    }
    
    /// Status of this routine (VALIDATED or UNVALIDATED)
    string status() const{
      return "UNVALIDATED";
    }
    
    /// Journal, and preprint references.
    vector<string> references() const{
      vector<string> refs;
      refs.push_back( "Phys.Rev.D75:092006,2007");
      refs.push_back( "Erratum-ibid.D75:119901,2007");
      refs.push_back( "FNAL-PUB 07/026-E" );
      refs.push_back( "hep-ex/0701051" );
      return refs;
    }
    
    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event & event);
    void finalize();
    //@}
    
  private:

    /// Rapidity range of histograms for R=0.05 and R=1 kt jets
    const double _minY, _maxY;
        
    /// Min jet \f$ p_T \f$ cut.
    /// @todo Make static const and UPPERCASE?
    const double _jetMinPT;
    
    /// Counter for the number of events analysed (actually the sum of weights, hence double).
    double _eventsTried;

    /// @name Histograms
    //@{
    /// The number of events in each histogram
    map<AIDA::IHistogram1D*, double> _eventsPassed;

    /// The y bin width of each histogram
    map<AIDA::IHistogram1D*, double> _yBinWidths;

    /// The y bin edge values
    static const double _ybins[6];

    /// Histograms in different eta regions
    BinnedHistogram<double> _binnedHistosD07;

    // Single histogram for the \f$R=0.5\f$ \f$k_\perp\f$ jets
    AIDA::IHistogram1D* _histoD05;

    // Single histogram for the \f$R=1.0\f$ \f$k_\perp\f$ jets
    AIDA::IHistogram1D* _histoD10;
    //@}

  };


}

#endif
