// -*- C++ -*-
#ifndef RIVET_DISKinematics_H
#define RIVET_DISKinematics_H

#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/RivetCLHEP.hh"

namespace Rivet {

  /// This class projects out the DIS kinematic variables and relevant
  /// boosts for an event.
  class DISKinematics: public Projection {

  public:
        
    /// The default constructor. Must specify Beam and DISLepton
    /// projection objects which are guaranteed to live throughout the
    /// run. Also the PDG code of the incoming hadron (\a hadid) must be
    /// specified.
    inline DISKinematics(Beam& beamp, DISLepton& leptonp, const ParticleName& hadid)
      : _beams(&beamp), _lepton(&leptonp), _idhad(hadid), 
        _theQ2(-1.0), _theW2(-1.0), _theX(-1.0) 
    {
      addBeamPair(ANY, hadid);
      addProjection(beamp);
      addProjection(leptonp);
    }
    
  public:
    /// Return the name of the projection
    inline string getName() const {
      return "DISKinematics";
    }
    
  protected:
    
    /// Perform the projection operation on the supplied event.
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual int compare(const Projection& p) const;

  public:

    /// The \f$Q^2\f$.
    inline double Q2() const { return _theQ2; }

    /// The \f$W^2\f$.
    inline double W2() const { return _theW2; }

    /// The Bjorken \f$x\f$.
    inline double x() const { return _theX; }

    /// The LorentzRotation needed to boost a particle to the hadronic CM frame.
    inline const LorentzRotation& boostHCM() const {
      return _hcm; 
    }

    /// The LorentzRotation needed to boost a particle to the hadronic Breit frame.
    inline const LorentzRotation& boostBreit() const {
      return _breit;
    }

  private:

    /// The Beam projector object defining the incoming beam particles.
    Beam* _beams;

    /// The projector for the scattered lepton.
    DISLepton* _lepton;

    /// The PDG id of the incoming hadron.
    long _idhad;

    /// The \f$Q^2\f$.
    double _theQ2;

    /// The \f$W^2\f$.
    double _theW2;

    /// The Bjorken \f$x\f$.
    double _theX;

    /// The LorentzRotation needed to boost a particle to the hadronic CM frame.
    LorentzRotation _hcm;

    /// The LorentzRotation needed to boost a particle to the hadronic Breit frame.
    LorentzRotation _breit;

//   private:

//     /// Hide the assignment operator.
//     DISKinematics & operator=(const DISKinematics &);

  };

}

#endif
