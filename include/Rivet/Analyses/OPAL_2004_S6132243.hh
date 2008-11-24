// -*- C++ -*-
#ifndef RIVET_OPAL_2004_S6132243_HH
#define RIVET_OPAL_2004_S6132243_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

//   /// 
//   class OPAL_2004_S6132243 : public Analysis {

//   public:

//     /// Default constructor.
//     OPAL_2004_S6132243()
//       : mult(fsproj), spher(fsproj) 
//     { 
//       setBeams(ELECTRON, POSITRON); 
//       addProjection(fsproj);
//       addProjection(mult);
//       addProjection(spher);
//     }

//   public:

  //    /// Factory method
  //    static Analysis* create() { return new TestAnalysis(); }

//     /// The name of this analysis is "OPAL_2004_S6132243"
//     string name() const {
//       return "OPAL_2004_S6132243";
//     }

//     virtual void init();

//     virtual void analyze(const Event & event);

//     virtual void finalize();

//     /// Return the RivetInfo object of this analysis object.
//     //virtual RivetInfo getInfo() const;

//   private:

//     /// The projectors used by this analysis.
//     FinalState fsproj;

//     Multiplicity mult;

//     Sphericity spher;

//   private:

//     /// Hide the assignment operator
//     OPAL_2004_S6132243 & operator=(const OPAL_2004_S6132243&);

//     /// @name Histograms
//     //@{
//     AIDA::IHistogram1D* histChTot_;
//     AIDA::IHistogram1D* histSphericity_;
//     AIDA::IHistogram1D* histPlanarity_;
//     AIDA::IHistogram1D* histAplanarity_;
//     //@}

//   };

}


#endif
