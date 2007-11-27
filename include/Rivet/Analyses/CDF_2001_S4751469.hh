// -*- C++ -*-
#ifndef RIVET_CDF_2001_S4751469_HH
#define RIVET_CDF_2001_S4751469_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/TrackJet.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  class CDF_2001_S4751469 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Default constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.5 \f$ GeV.
    CDF_2001_S4751469()
      : _fsproj(-1.0, 1.0, 0.5), _trackjetproj(_fsproj),
        _dataToward(0) ,_dataTrans(0), _dataAway(0)
    { 
      setBeams(PROTON, ANTIPROTON);
      addProjection(_fsproj);
      addProjection(_trackjetproj);
    }

    /// Factory method
    static Analysis* create() { return new CDF_2001_S4751469(); }

    //@}

  public:

    /// @name Metadata
    //@{

    /// Get the name of the analysis.
    string getName() const { 
      return getExpt() + "_" + getYear() + "_S" + getSpiresId();
    }
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "4751469";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Field & Stuart underlying event analysis at CDF.";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "CDF";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "2001";
    }
    /// Journal, and preprint references.
    vector<string> getReferences() const {
      vector<string> ret;
      ret.push_back("Phys.Rev.D65:092002,2002");
      ret.push_back("FNAL-PUB 01/211-E");
      // No arXiv code
      return ret;
    }

    //@}

  public:

    /// @name Analysis methods
    //@{
    void init();
    void analyze(const Event& event);
    void finalize();
    //@}

  private:

    /// @name Internal projections
    //@{
    /// The FinalState projection used.
    ChargedFinalState _fsproj;
    /// The TrackJet projection used by this analysis.
    TrackJet _trackjetproj;
    //@}

  private:

    /// @name Histogram collections
    //@{
    /// Profile histograms, binned in the \f$ p_T \f$ of the
    /// leading jet, for the \f$ p_T \f$ sum in the toward, 
    /// transverse and away regions.
    AIDA::IProfile1D* _dataToward;
    AIDA::IProfile1D* _dataTrans;
    AIDA::IProfile1D* _dataAway;
    //@}


  private:

    /// Hide the assignment operator.
    CDF_2001_S4751469& operator=(const CDF_2001_S4751469&);

  };

}

#endif
