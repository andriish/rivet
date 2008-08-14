// -*- C++ -*-
#ifndef RIVET_TrackJet_HH
#define RIVET_TrackJet_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Jet.hh"

namespace Rivet {


  /// Build jets using the non-IR-safe algorithm of the CDF Field/Stuart UE analysis.
  class TrackJet : public JetAlg {
  public:

    /// @name Standard constructors and destructors.
    //@{

    /// Constructor. The specified FinalState projection is assumed to live
    /// throughout the run and should be used to specify the max and min \f$
    /// \eta \f$ values and the min \f$ p_T \f$ (in GeV).
    TrackJet(const FinalState& fsp, const double Rmax=0.7)
      : _Rmax(Rmax)
    {
      setName("TrackJet");
      addProjection(fsp, "FS");
    }


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new TrackJet(*this);
    }
    //@}

    /// Typedef for the tracks (a list so that elements can be consistently removed
    typedef list<FourMomentum> Tracks;

    /// Get the computed jets.
    Jets getJets(double ptmin=0.0) const {
      if (ptmin == 0.0) return _jets;
      Jets rtn;
      for (Jets::const_iterator j = _jets.begin(); j != _jets.end(); ++j) {
        if (j->ptSum() >= ptmin) {
          rtn.push_back(*j);
        }
      }
      return rtn;
    }

  protected:   

    /// Perform the projection on the Event. The collection of jets that results
    /// will be sorted in order of decreasing \f$ p_T \f$.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const; 

  private:
    
    /// The computed jets
    Jets _jets;

    /// \f$ R = \sqrt{\eta^2 + \phi^2} \f$ cut in jet definition.
    double _Rmax;
  };
  
 
}

#endif
