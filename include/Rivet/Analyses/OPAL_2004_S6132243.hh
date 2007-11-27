// -*- C++ -*-
#ifndef RIVET_EURPHYS40C287_HH
#define RIVET_EURPHYS40C287_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Multiplicity.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/RivetAIDA.fhh"

namespace Rivet {

//   /// This class...
//   class S2435284 : public Analysis {

//   public:

//     /// Default constructor.
//     inline S2435284()
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

//     /// The name of this analysis is "S2435284"
//     inline string getName() const {
//       return "S2435284";
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
//     S2435284 & operator=(const S2435284 &);

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
