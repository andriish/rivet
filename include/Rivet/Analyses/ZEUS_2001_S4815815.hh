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
    ZEUS_2001_S4815815() { 
      setBeams(POSITRON, PROTON);
      addProjection(*new FastJets(*new FinalState()), "Jets");
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
    string getDescription() const {
      return "Dijet photoproduction analysis";
    }
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

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histJetEt1;
    //@}

    /// Hide the assignment operator
    ZEUS_2001_S4815815& operator=(const ZEUS_2001_S4815815&);
  };

}

#endif
