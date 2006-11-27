#ifndef MATH_LORENTZVECTOR_HH 
#define MATH_LORENTZVECTOR_HH 1

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

namespace Rivet {

  /// Representation of a Lorentz vector.
  ///
  /// This vector is simpler than a CLHEP (Hep)LorentzVector and performs 
  /// some internal magic to manage light cone type and Cartesian type
  /// internal representations in an equal-handed way. It's based on a class
  /// coded up for KtJet by Andy Buckley.
  ///
  /// @todo Re-implement, using vectors based on expression templates cf.
  /// Boost::uBLAS or Blitz++.
  ///
  /// @todo Define operator+ and operator- (and unary minus?), *(1,2) and (1,2).
  ///
  /// @todo Make external manipulators for Minkowski dot product and tensor product.
  ///
  /// @todo Make ThreeVector class ... similar base functionality to this plus a cross product.
  ///
  /// @author Andy Buckley
  /// @date   2006-11-27
  /// $Id: $
  ///
  class LorentzVector {
  public:
    /// Default Constructor: create null vector
    LorentzVector();

    /// Cartesian constructor
    LorentzVector(double ct, double x, double y, double z);

    /// Destructor
    ~LorentzVector();

    /// @brief Set components in Cartesian representation.
    /// Note that \a c = 1, so \a ct can be treated as just \a t.
    inline void setCartesian(double ct, double x, double y, double z);

    /// Set components in Snowmass light cone representation.
    inline void setSnowmass(double rapidity, double phi, double perp, double m = 0);

    /// Clears the vector
    inline double lorentzInvariant() const;
 
    /// Clears the vector
    inline void clear();

    /// Compare vectors: reasonable, but should this just mean an address comparison?
    //inline bool operator== (const LorentzVector&) const;

    /// Compare vectors: reasonable, but should this just mean an address comparison?
    //inline bool operator!= (const LorentzVector&) const;

    /// Compare vectors: not obvious
    //inline bool operator< (const LorentzVector&) const;

    /// Compare vectors: not obvious
    //inline bool operator> (const LorentzVector&) const;

    inline LorentzVector& operator+=(const LorentzVector& other);

    inline LorentzVector& operator-=(const LorentzVector& other);

    inline std::string toString() const;

    static const LorentzVector X;
    static const LorentzVector Y;
    static const LorentzVector Z;

    /// Getters for Cartesian variables
    inline double ct() const;
    inline double x() const;
    inline double y() const;
    inline double z() const;

    /// @todo Setters for Cartesian variables

    /// Getters for Snowmass variables
    inline double rapidity() const;
    inline double phi() const;
    inline double perp() const;
    inline double perp2() const;
    inline double m() const;

    /// @todo Setters for Snowmass variables

  protected:    
    /// Calculate Snowmass co-ords from Cartesian ones
    inline void calcSnow() const;

    /// Calculate Cartesian co-ords from Snowmass ones
    inline void calcCart() const;


  private:
    mutable double ct_;
    mutable double x_;
    mutable double y_;
    mutable double z_;

    mutable double rapidity_;
    mutable double phi_;
    mutable double perp_;
    mutable double perp2_;
    mutable double m_;

    enum ConsistencyStatus { CONSISTENT, BADSNOW, BADCART };
    mutable ConsistencyStatus status_;

  };

  LorentzVector operator*(double num, const LorentzVector& vec);

  LorentzVector operator*(const LorentzVector& vec, double num);
  
  LorentzVector operator/(const LorentzVector& vec, double num);

  std::ostream& operator<<(std::ostream& os, LorentzVector vec);

} // end of namespace


#include "Rivet/Tools/Math/LorentzVector.icc"


#endif // MATH_LORENTZVECTOR_HH
