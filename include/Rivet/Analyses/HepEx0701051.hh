// -*- C++ -*-
#ifndef RIVET_HepEx0701051_HH
#define RIVET_HepEx0701051_HH

#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.fhh"

#include "Rivet/Projections/KtJets.hh"
#ifdef HAVE_FASTJET	
#include "Rivet/Projections/FastJets.hh"
#endif


namespace Rivet {


  /// @todo Needs full Doxygen commenting
  class HepEx0701051 : public Analysis {

  public:

    /// Default constructor
    HepEx0701051() {}

    /// Destructor
    ~HepEx0701051() {}

  public:

    /// Factory method
    static Analysis* create() { 
      return new HepEx0701051(); 
    }

    /// Get the name of this analysis.
    inline string getName() const {
      return "HepEx0701051";
    }

  public:

    // Declaration of initialization method
    void init();

    // Declaration of analyzing method
    void analyze(const Event & event);
    
    void finalize();

  private:

    /// Hide the assignment operator
    HepEx0701051& operator=(const HepEx0701051&);

    //@{
    /// Histograms
    /// @todo Private members should have an "_" prefix. Convention for histograms is _histName.
    AIDA::IHistogram1D* pHistogramObject1;
    AIDA::IHistogram1D* pHistogramObject1N;
    AIDA::IHistogram1D* pHistogramObject2;
    AIDA::IHistogram1D* pHistogramObject2N;
    AIDA::IHistogram1D* pHistogramObject3;
    AIDA::IHistogram1D* pHistogramObject3N;
    AIDA::IHistogram1D* pHistogramObject4;
    AIDA::IHistogram1D* pHistogramObject4N;
    AIDA::IHistogram1D* pHistogramObject5;
    AIDA::IHistogram1D* pHistogramObject5N;
    //@}

  };

}

#endif
