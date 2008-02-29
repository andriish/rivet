// -*- C++ -*-
#ifndef RIVET_CDF_2001_S4751469_HH
#define RIVET_CDF_2001_S4751469_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/TrackJet.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/LossyFinalState.hh"


namespace Rivet {

  class CDF_2001_S4751469 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor: cuts on final state are \f$ -1 < \eta < 1 \f$ 
    /// and \f$ p_T > 0.5 \f$ GeV. Use a lossy charged FS projection, which
    /// randomly discards 8% of charged particles, as a kind of hacky detector 
    /// correction.
    CDF_2001_S4751469()
      : _cfsproj(-1.0, 1.0, 0.5*GeV), _fsproj(_cfsproj, 0.08), _trackjetproj(_fsproj),
        _ptsumToward(0) ,_ptsumTrans(0), _ptsumAway(0),
        _numToward(0) ,_numTrans(0), _numAway(0)
    { 
      setBeams(PROTON, ANTIPROTON);
      addProjection(_cfsproj);
      addProjection(_fsproj);
      addProjection(_trackjetproj);
    }

    /// Factory method
    static Analysis* create() { return new CDF_2001_S4751469(); }

    //@}

  public:

    /// @name Publication metadata
    //@{
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
    /// Use only charged tracks
    ChargedFinalState _cfsproj;
    /// Lose 8% of charged tracks randomly.
    LossyFinalState _fsproj;
    /// The TrackJet projection used by this analysis.
    TrackJet _trackjetproj;
    //@}

  private:

    /// @name Histogram collections
    //@{
    /// Profile histograms, binned in the \f$ p_T \f$ of the leading jet, for
    /// the \f$ p_T \f$ sum in the toward, transverse and away regions.
    AIDA::IProfile1D* _ptsumToward;
    AIDA::IProfile1D* _ptsumTrans;
    AIDA::IProfile1D* _ptsumAway;

    /// Profile histograms, binned in the \f$ p_T \f$ of the leading jet, for
    /// the number of charged particles per jet in the toward, transverse and
    /// away regions.
    AIDA::IProfile1D* _numToward;
    AIDA::IProfile1D* _numTrans;
    AIDA::IProfile1D* _numAway;
    //@}


  private:

    /// Hide the assignment operator.
    CDF_2001_S4751469& operator=(const CDF_2001_S4751469&);

  };

}

#endif
