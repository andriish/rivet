#ifndef RIVET_MATH_LORENTZTRANS
#define RIVET_MATH_LORENTZTRANS

#include "Rivet/Math/MathHeader.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/MatrixN.hh"
#include "Rivet/Math/Matrix3.hh"
#include "Rivet/Math/Vector4.hh"
#include <iostream>

namespace Rivet {


  /// @name Simple Lorentz factor conversions
  //@{

  /// Calculate the gamma factor from beta
  inline double lorentzBeta2Gamma(double beta) {
    return 1.0 / sqrt(1 - sqr(beta));
  }

  /// Calculate beta from the gamma factor
  inline double lorentzGamma2Beta(double gamma) {
    return sqrt(1 - sqr(1/gamma));
  }

  //@}


  /// @brief Object implementing Lorentz transform calculations and boosts.
  ///
  /// @note These boosts are defined actively, i.e. as modifications of vectors
  /// rather than frame transformations. So the boost vector is the opposite of
  /// what you might expect.
  ///
  /// @todo Review the active/passive convention choice. Seems counterintuitive now...
  class LorentzTransform {
  public:

    // friend string toString(const LorentzTransform& lt); ///< @todo Ick!

    LorentzTransform() {
      _boostMatrix = Matrix<4>::mkIdentity();
    }

    LorentzTransform(const Vector3& boost) {
      setBoost(boost);
    }

    LorentzTransform(double betaX, double betaY, double betaZ) {
      setBoost(betaX, betaY, betaZ);
    }

    /// Set the \f$ \vec\beta \f$ vector for an active Lorentz boost
    LorentzTransform& setBoost(const Vector3& boost) {
      assert(boost.mod2() < 1);
      const double beta = boost.mod();
      const double gamma = lorentzBeta2Gamma(beta);
      _boostMatrix = Matrix<4>::mkIdentity();
      _boostMatrix.set(0, 0, gamma);
      _boostMatrix.set(1, 1, gamma);
      // Positive coeff since these are active boosts
      _boostMatrix.set(0, 1, +beta*gamma);
      _boostMatrix.set(1, 0, +beta*gamma);
      if (beta > 0) _boostMatrix = rotate(Vector3::mkX(), boost)._boostMatrix;
      return *this;
    }

    /// Set the \f$ \vec\beta \f$ vector for an active Lorentz boost, as 3-components
    LorentzTransform& setBoost(double betaX, double betaY, double betaZ) {
      return setBoost(Vector3(betaX, betaY, betaZ));
    }

    /// Get the \f$ \vec\beta \f$ vector for an active Lorentz boost
    Vector3 boost() const {
      FourMomentum boost(_boostMatrix.getColumn(0));
      //cout << "!!!" << boost << endl;
      if (boost.isZero()) return boost;
      assert(boost.E() > 0);
      const double beta = boost.p3().mod() / boost.E();
      return boost.p3().unit() * beta;
    }

    /// Get the \f$ \beta \f$ factor
    double beta() const {
      return boost().mod();
    }

    /// Get the \f$ \gamma \f$ factor
    double gamma() const {
      return lorentzBeta2Gamma(beta());
    }

    LorentzTransform rotate(const Vector3& from, const Vector3& to) const {
      return rotate(Matrix3(from, to));
    }

    LorentzTransform rotate(const Vector3& axis, double angle) const {
      return rotate(Matrix3(axis, angle));
    }

    LorentzTransform rotate(const Matrix3& rot) const {
      LorentzTransform lt = *this;
      const Matrix4 rot4 = _mkMatrix4(rot);
      const Matrix4 newlt = rot4 * _boostMatrix * rot4.inverse();
      lt._boostMatrix = newlt;
      return lt;
    }

    /// Apply this transformation to the given 4-vector
    FourVector transform(const FourVector& v4) const {
      return multiply(_boostMatrix, v4);
    }

    /// Calculate the inverse transform
    LorentzTransform inverse() const {
      LorentzTransform rtn;
      rtn._boostMatrix = _boostMatrix.inverse();
      return rtn;
    }


    /// Combine LTs, treating @a this as the LH matrix.
    LorentzTransform combine(const LorentzTransform& lt) const {
      LorentzTransform rtn;
      rtn._boostMatrix = _boostMatrix * lt._boostMatrix;
      return rtn;
    }

    /// Return the matrix form
    Matrix4 toMatrix() const {
      return _boostMatrix;
    }

    /// Operator combination of two LTs
    LorentzTransform operator*(const LorentzTransform& lt) const {
      return combine(lt);
    }

    /// Pre-multiply m3 by this LT
    LorentzTransform preMult(const Matrix3& m3) {
      _boostMatrix = multiply(_mkMatrix4(m3),_boostMatrix);
      return *this;
    }

    /// Post-multiply m3 by this LT
    LorentzTransform postMult(const Matrix3& m3) {
      _boostMatrix *= _mkMatrix4(m3);
      return *this;
    }


  private:

    Matrix4 _mkMatrix4(const Matrix3& m3) const {
      Matrix4 m4 = Matrix4::mkIdentity();
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          m4.set(i+1, j+1, m3.get(i, j));
        }
      }
      return m4;
    }

    Matrix4 _boostMatrix;

  };



  inline LorentzTransform inverse(const LorentzTransform& lt) {
    return lt.inverse();
  }

  inline LorentzTransform combine(const LorentzTransform& a, const LorentzTransform& b) {
    return a.combine(b);
  }

  inline FourVector transform(const LorentzTransform& lt, const FourVector& v4) {
      return lt.transform(v4);
  }


  //////////////////////////


  inline string toString(const LorentzTransform& lt) {
    return toString(lt.toMatrix());
  }

  inline ostream& operator<<(std::ostream& out, const LorentzTransform& lt) {
    out << toString(lt);
    return out;
  }


}

#endif
