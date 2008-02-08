// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Tools/Logging.hh"


namespace Rivet {

  void Thrust::calcT(const vector<Vector3>& p, double& t, Vector3& taxis) const {
    t = 0.0;
    Vector3 tv, ptot;
    vector<Vector3> cpm(4);
    const size_t psize = p.size();
    for (size_t k = 1; k < psize; ++k) {
      for (size_t j = 0; j < k; ++j) {
        tv = p[j].cross(p[k]);
        ptot = Vector3();
        for (size_t l = 0; l < psize; ++l) {
          if (l != j && l != k) {
            if (p[l].dot(tv) > 0.0) { 
              ptot += p[l];
            } else {
              ptot -= p[l];
            }
          }
        }
        cpm[0] = ptot - p[j] - p[k];
        cpm[1] = ptot - p[j] + p[k];
        cpm[2] = ptot + p[j] - p[k];
        cpm[3] = ptot + p[j] + p[k];
        for (size_t i = 0; i < 4; ++i) {
          double tval = mod2(cpm[i]);
          if (tval > t) {
            t = tval;
            taxis = cpm[i];
          }
        }
      } // j loop
    } // k loop
  }



  void Thrust::calcM(const vector<Vector3>& p, double& m, Vector3& maxis) const {
    double mval;
    m = 0.0;
    Vector3 tv, ptot;
    vector<Vector3> cpm;
    for (unsigned int j = 0; j < p.size(); ++j) {
      tv = p[j];
      ptot = Vector3();
      for (unsigned int l = 0; l < p.size(); ++l) {
        if (l != j) {
          if (p[l].dot(tv) > 0.0) { 
            ptot += p[l];
          } else {
            ptot -= p[l];
          }
        }
      }
      cpm.clear();
      cpm.push_back(ptot - p[j]);
      cpm.push_back(ptot + p[j]);
      for (vector<Vector3>::iterator it = cpm.begin(); it != cpm.end(); ++it) {
        mval = it->mod2();
        if (mval > m) {
          m = mval;
          maxis = *it;
        }
      }
    } // j loop
  }

  inline bool mod2ReverseCmp(const Vector3& a, const Vector3& b) {
    return a.mod2() < b.mod2();
  }

  void Thrust::calcThrust(const FinalState& fs) {
    // Make a vector of the three-momenta in the final state
    vector<Vector3> threeMomenta;
    double momentumSum(0.0);
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      Vector3 p3 = p->getMomentum().vector3();
      threeMomenta.push_back(p3);
      momentumSum += mod(threeMomenta.back());
    }


    // Clear the caches
    _thrusts.clear();
    _thrustAxes.clear(); 


    // If there are fewer than 2 visible particles, we can't do much
    if (threeMomenta.size() < 2) {
      for (int i = 0; i < 3; ++i) {
        _thrusts.push_back(-1);
        _thrustAxes.push_back(Vector3());
      }
      return;
    }


    // Handle special case of thrust = 1 if there are only 2 particles
    if (threeMomenta.size() == 2) {
      Vector3 axis;
      _thrusts.push_back(1.0);
      _thrusts.push_back(0.0);
      _thrusts.push_back(0.0);
      axis = threeMomenta[0].unit();
      if (axis.z() < 0) axis = -axis;
      _thrustAxes.push_back(axis);
      /// @todo Improve this --- special directions bad...
      /// (a,b,c) _|_ 1/(a^2+b^2) (b,-a,0) etc., but which combination minimises error?
      _thrustAxes.push_back( axis.cross(Vector3(0,0,1)) );
      _thrustAxes.push_back( _thrustAxes[0].cross(_thrustAxes[1]) );
      return;
    }


    // Handle special case where thrust = max{x1,x2,x3} if there are 3 particles
    if (threeMomenta.size() == 3) {
      Vector3 axis;
      // Order by magnitude
      std::sort(threeMomenta.begin(), threeMomenta.end(), mod2ReverseCmp);
      //std::reverse(threeMomenta.begin(), threeMomenta.end());
      // Thrust
      axis = threeMomenta[0].unit();
      if (axis.z() < 0) axis = -axis;
      _thrusts.push_back(2.0 * mod(threeMomenta[0]) / momentumSum);
      _thrustAxes.push_back(axis);
      // Thrust major (independent part of next largest momentum)
      axis = ( threeMomenta[1] - dot(axis, threeMomenta[1]) * axis ).unit();
      if (axis.x() < 0) axis = -axis;
      _thrusts.push_back( (fabs(threeMomenta[1].dot(axis)) + 
                           fabs(threeMomenta[2].dot(axis)))  / momentumSum);
      _thrustAxes.push_back(axis);
      // Thrust minor
      _thrusts.push_back(0.0);
      _thrustAxes.push_back( _thrustAxes[0].cross(_thrustAxes[1]) );
      return;
    }


    // If the special cases don't apply, we have to use a general method. Here
    // we explicitly calculate in units of MeV using an algorithm based on 
    // Brandt/Dahmen Z Phys C1 (1978) and 'tasso' code from HERWIG. Re-coded
    // from Herwig++ implementation by Stefan Gieseke by Andy Buckley.
    // NB. special case with >= 4 coplanar particles will still fail.
    // probably not too important... 

    // Temporary variables for calcs
    Vector3 axis; double val;

    // Get thrust
    calcT(threeMomenta, val, axis);
    getLog() << Log::DEBUG << "Mom sum = " << momentumSum << endl;
    _thrusts.push_back(sqrt(val) / momentumSum);
    // Make sure that thrust always points along the +ve z-axis.
    if (axis.z() < 0) axis = -axis;
    axis = axis.unit();
    getLog() << Log::DEBUG << "Axis = " << axis << endl;
    _thrustAxes.push_back(axis);

    // Get thrust major 
    threeMomenta.clear();
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      // Get the part of each 3-momentum which is perpendicular to the thrust axis
      const Vector3 v = p->getMomentum().vector3();
      const Vector3 vpar = dot(v, axis.unit()) * axis.unit();
      threeMomenta.push_back(v - vpar);
    }
    calcM(threeMomenta, val, axis);
    _thrusts.push_back(sqrt(val) / momentumSum);
    if (axis.x() < 0) axis = -axis;
    axis = axis.unit();
    _thrustAxes.push_back(axis); 

    // Get thrust minor
    if (_thrustAxes[0].dot(_thrustAxes[1]) < 1e-10) {
      axis = _thrustAxes[0].cross(_thrustAxes[1]);
      _thrustAxes.push_back(axis);
      val = 0.0;
      for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
        val += fabs(dot(axis, p->getMomentum().vector3()));
      }
      _thrusts.push_back(val / momentumSum);
    } else {
      _thrusts.push_back(-1.0);
      _thrustAxes.push_back(Vector3());
    }

  }


  void Thrust::project(const Event& e) {
    const FinalState& fs = e.applyProjection(_fsproj);
    calcThrust(fs);
  }


}
