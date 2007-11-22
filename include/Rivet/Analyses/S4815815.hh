// -*- C++ -*-
#ifndef RIVET_S4815815_HH
#define RIVET_S4815815_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/KtJets.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {

  /// This class is a reproduction of the HZTool routine for the ZEUS 
  /// dijet photoproduction paper which was used in the ZEUS Jets PDF fit.  
  ///
  /// @author Jon Butterworth
  class S4815815 : public Analysis {

  public:

    /// Default constructor.
    inline S4815815()
      : _fsproj(), _ktjetsproj(_fsproj) 
    { 
      setBeams(POSITRON, PROTON);
      addProjection(_fsproj);
      addProjection(_ktjetsproj);
    }

    /// Factory method
    static Analysis* create() { return new S4815815(); }

    /// Get the name of this analysis.
    inline string getName() const {
      return "S4815815";
    }

  public:

    void init();
    
    void analyze(const Event& event);
    
    void finalize();

  private:

    /// The FinalState projector used by this analysis.
    FinalState _fsproj;

    /// The KtJets projector used by this analysis.
    KtJets _ktjetsproj;

  private:

    /// Hide the assignment operator
    S4815815& operator=(const S4815815&);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histJetEt1;
    //@}

  };

}

#endif
