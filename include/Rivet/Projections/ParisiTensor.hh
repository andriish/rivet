// -*- C++ -*-
#ifndef RIVET_ParisiTensor_HH
#define RIVET_ParisiTensor_HH

#include "Rivet/Projections/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Event.hh"


namespace Rivet {

  /**
     @brief Calculate the Parisi event shape tensor (or linear momentum tensor).
     
     The Parisi event shape C and D variables are derived from the eigenvalues of
     the linear momentum tensor
     \f[ 
     \theta^{\alpha \beta} = 
     \frac{\sum_i \frac{p_i^\alpha p_i^\beta}{|\mathbf{p}_i|}}
          {\sum_i |\mathbf{p}_i|} 
     \f]
     which is actually a linearized (and hence infra-red safe) version of the 
     {@link Sphericity} tensor.

     Defining the three eigenvalues of \f$\theta\f$
     \f$ \lambda_1 \ge \lambda_2 \ge \lambda_3 \f$, with \f$ \lambda_1 + \lambda_2 + \lambda_3 = 1 \f$, 
     the C and D parameters are defined as
     \f[ 
     C = 3(\lambda_1\lambda_2 + \lambda_1\lambda_3 + \lambda_2\lambda_3)
     \f]
     and
     \f[ 
     D = 27 \lambda_1\lambda_2\lambda_3
     \f]

     Internally, this Projection uses the Sphericity projection with the generalising
     \f$r\f$ parameter set to 1.
  */
  class ParisiTensor : public Projection {

  public:

    /// Constructor. The provided FinalState projection must live throughout the run.
    inline ParisiTensor(FinalState& fsp)
      : _C(0), _D(0), _sphproj(Sphericity(fsp, 1.0))
    { }

  public:
    /// Return the name of the projection
    inline string name() const {
      return "ParisiTensor";
    }

  protected:

    /// Perform the projection on the Event.
    void project(const Event& e);

    /// Compare with other projections.
    int compare(const Projection& p) const;

  public:

    /// @name Access the C and D params.
    ///@{
    inline const double C() const { return _C; }
    inline const double D() const { return _D; }
    ///@}

    /// @name Access the eigenvalues of \f$\theta\f$.
    ///@{
    inline const double lambda1() const { return _lambda[0]; }
    inline const double lambda2() const { return _lambda[1]; }
    inline const double lambda3() const { return _lambda[2]; }
    ///@}
    
//     /// Return the RivetInfo object of this Projection.
//     virtual RivetInfo getInfo() const;
    
  private:
    
    /// The event shapes
    double _C, _D;

    /// Eigenvalues
    double _lambda[3];

    /// The Sphericity projection which this projection is really just a wrapper for.
    Sphericity _sphproj;

  private:

    /// The assignment operator is private and must never be called.
    /// In fact, it shouldn't even be implemented.
    ParisiTensor & operator=(const ParisiTensor &);

  };

}


#endif /* RIVET_ParisiTensor_HH */
