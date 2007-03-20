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
     @brief Calculate the sphericity event shape.
     
     The sphericity tensor (or quadratic momentum tensor) is defined as 
     \f[ 
     S^{\alpha \beta} = \frac{\sum_i p_i^\alpha p_i^\beta}{\sum_i |\mathbf{p}_i|^2} 
     \f],
     from which the sphericity, aplanarity and planarity can be
     calculated by combinations of eigenvalues.
     
     Defining the three eigenvalues
     \f$ \lambda_1 \ge \lambda_2 \ge \lambda_3 \f$, with \f$ \lambda_1 + \lambda_2 + \lambda_3 = 1 \f$, 
     the sphericity is
     \f[ 
     S = \frac{3}{2} (\lambda_2 + \lambda_3) 
     \f].
     
     The aplanarity is \f$ A = \frac{3}{2}\lambda_3 \f$ and the planarity
     can be obtained through the trace constraint. The eigenvectors define a 
     set of spatial axes comparable with the thrust axes, but more sensitive to 
     high momentum particles due to the quadratic sensitivity of the tensor to
     the particle momenta.
     
     Since the sphericity is quadratic in the particle momenta, it is not an
     infrared safe observable in perturbative QCD. This can be fixed by adding
     a regularizing power of \f$r\f$ to the definition:
     \f[ 
     S^{\alpha \beta} = 
     \frac{\sum_i |\mathbf{p}_i|^{r-2} \sum_i p_i^\alpha p_i^\beta}
     {\sum_i |\mathbf{p}_i|^r} 
     \f].
     
     \f$r\f$ is available as a constructor argument on this class and will be
     taken into account by the Cmp<Projection> operation, so a single analysis
     can use several sphericity projections with different \f$r\f$ values without
     fear of a clash.
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

    inline const double lambda1() const { return _lambdas[0]; }
    inline const double lambda2() const { return _lambdas[1]; }
    inline const double lambda3() const { return _lambdas[2]; }

  /**
   * Return the RivetInfo object of this Projection. Derived classes
   * should re-implement this function to return the combined
   * RivetInfo object of this and of any other Projection upon which
   * this depends.
   */
  virtual RivetInfo getInfo() const;

  private:

    /// Eigenvalues
    double _lambdas[3];

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
