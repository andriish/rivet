// -*- C++ -*-
#ifndef RIVET_HepEx0112029_H
#define RIVET_HepEx0112029_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/KtJets.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {

  /// This class is a reproduction of the HZTool routine for the ZEUS 
  /// dijet photoproduction paper which was used in the ZEUS Jets PDF fit.  
  ///
  /// @author Jon Butterworth
  class HepEx0112029 : public Analysis {

  public:

    /// Default constructor.
    inline HepEx0112029()
      : _fsproj(), _ktjetsproj(_fsproj) 
    { 
      setBeams(POSITRON, PROTON);
      addProjection(_fsproj);
      addProjection(_ktjetsproj);
    }

    /// Factory method
    static Analysis* create() { return new HepEx0112029(); }

    /// Get the name of this analysis.
    inline string getName() const {
      return "HepEx0112029";
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
    HepEx0112029& operator=(const HepEx0112029&);

    /// @name Histograms
    //@{
    AIDA::IHistogram1D* _histJetEt1;
    //@}

  };

}

#endif
