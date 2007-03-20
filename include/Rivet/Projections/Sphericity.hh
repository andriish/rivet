// -*- C++ -*-
#ifndef RIVET_Sphericity_HH
#define RIVET_Sphericity_HH
// Declaration of the Sphericity class.

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/Event/Particle.hh"
#include "Rivet/Tools/Event/Event.hh"


namespace Rivet {

  /**
   * Project out the event shape.
   */
  class Sphericity : public Projection {

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// Default constructor. Must specify a FinalState projection which is
    //assumed to live throughout the run.
    inline Sphericity(FinalState& fsp, double rparam=2)
      : _sphericity(0), _planarity(0), _aplanarity(0), _regparam(rparam), 
        _fsproj(&fsp)
    { }

    /// Copy constructor.
    inline Sphericity(const Sphericity& sp)
      : Projection(sp), 
        _sphericity(sp._sphericity), _planarity(sp._planarity),
        _aplanarity(sp._aplanarity), _regparam(sp._regparam),
        _fsproj(sp._fsproj)
    { }

    /// Destructor.
    virtual ~Sphericity() { }
    //@}

  protected:

    /// Perform the projection on the Event: only to be called by 
    /// Event::applyProjection(Projection &).
    void project(const Event& e);

    /// This function defines a unique ordering between different 
    /// Projection objects of the same class. Should only be called from 
    /// operator<(const Projection &).
    int compare(const Projection& p) const;

  public:

    /// name access the event shapes, 
    /// Sphericity, Planarity and APlanarity
    inline const double sphericity() const { return _sphericity; }
    inline const double planarity() const { return _planarity; }
    inline const double aplanarity() const { return _aplanarity; }

  /**
   * Return the RivetInfo object of this Projection. Derived classes
   * should re-implement this function to return the combined
   * RivetInfo object of this and of any other Projection upon which
   * this depends.
   */
  virtual RivetInfo getInfo() const;

  private:

    /// The event shapes
    double _sphericity, _planarity, _aplanarity;

    /// Regularizing parameter, used to force infra-red safety.
    const double _regparam;

    /// The FinalState projection used by this projection.
    FinalState* _fsproj;

  private:

    /**
     * The assignment operator is private and must never be called.
     * In fact, it shouldn't even be implemented.
     */
    Sphericity & operator=(const Sphericity &);

  };

}


#endif /* RIVET_Sphericity_HH */
