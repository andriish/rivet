// -*- C++ -*-
#ifndef RIVET_HepEx0112029_H
#define RIVET_HepEx0112029_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/KtJets.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {

  /**
   *  This class is under construction. It will be a reproduction of the HZTool routine
   *  for the ZEUS dijet photoproduction paper which was used in the ZEUS Jets PDF fit.
   *  And hopefully it will teach me how to write a rivet analysis.
   *  @author Jon Butterworth
   */
  class HepEx0112029 : public Analysis {

  public:

    /// Default constructor.
    inline HepEx0112029()
      : fs(), ktjets(fs) 
    { }

    /// Copy constructor.
    inline HepEx0112029(const HepEx0112029& x) 
      : Analysis(x), ktjets(x.ktjets) 
    { }

    /// Destructor
    inline ~HepEx0112029()
    { }

    /// The name of this analysis is "HepEx0112029"
    inline std::string name() const {
      return "HepEx0112029";
    }

  public:

    void init();
    
    void analyze(const Event & event);
    
    void finalize();

    /// Return the RivetInfo object of this analysis object.
    RivetInfo getInfo() const;

  private:

    /// The FinalState projector used by this analysis.
    FinalState fs;

    /// The KtJets projector used by this analysis.
    KtJets ktjets;


  private:

    /// Hide the assignment operator
    HepEx0112029 & operator=(const HepEx0112029& x);

    //@{
    /// Histograms
    AIDA::IHistogram1D* histJetEt1_;
    //@}

  };

}

#endif
