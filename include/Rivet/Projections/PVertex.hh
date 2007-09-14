// -*- C++ -*-
#ifndef RIVET_PVertex_HH
#define RIVET_PVertex_HH

#include "Rivet/Projection.hh"
#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"


namespace Rivet {
  
  
  /// Returns the primary vertex of an event
  class PVertex: public Projection {

  public:

    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor.
    inline PVertex() { 
    }

    /// The destructor.
    virtual ~PVertex() { }
    //@}

    public:
      /// Return the name of the projection
      inline string getName() const {
        return "PVertex";
      }

  protected:

    /// Do the projection.
    virtual void project(const Event & e);

    /// Compare projections.
    inline virtual int compare(const Projection & p) const {
      return 0;
    }

  public:

    /// Get the primary vertex.
    inline const GenVertex& getPrimaryVertex() const {
      return *_thePVertex;
    }

  private:

    /// The Primary Vertex in the current collision.
    GenVertex* _thePVertex;

  };

}


#endif
