// -*- C++ -*-
#ifndef RIVET_ZEUS_2001_S4815815_HH
#define RIVET_ZEUS_2001_S4815815_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This class is a reproduction of the HZTool routine for the ZEUS 
  /// dijet photoproduction paper which was used in the ZEUS Jets PDF fit.  
  ///
  /// @author Jon Butterworth
  class ZEUS_2001_S4815815 : public Analysis {

  public:

    /// Default constructor.
    ZEUS_2001_S4815815()
      : _fsproj(), _jetsproj(_fsproj) 
    { 
      setBeams(POSITRON, PROTON);
      addProjection(_fsproj);
      addProjection(_jetsproj);
    }

    /// Factory method.
    static Analysis* create() { 
      return new ZEUS_2001_S4815815(); 
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "4815815";
    }
    /// Get a description of the analysis.
    // string getDescription() const {
    //   return "ZEUS";
    // }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "ZEUS";
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

    /// The final state projection used by this analysis.
    FinalState _fsproj;

    /// The jets projection used by this analysis.
    FastJets _jetsproj;

  private:

    /// Hide the assignment operator
    ZEUS_2001_S4815815& operator=(const ZEUS_2001_S4815815&);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histJetEt1;
    //@}

  };

}

#endif
