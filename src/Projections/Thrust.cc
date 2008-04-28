// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Tools/Logging.hh"


namespace Rivet {

  inline bool mod2Cmp(const Vector3& a, const Vector3& b) {
    return a.mod2() > b.mod2();
  }


  void Thrust::calcT(vector<Vector3>& p, double& t, Vector3& taxis) const {
    /* This function implements the iterative algorithm as described in the
     * Pythia manual. We take eight (four) different starting vectors
     * constructed from the four (three) leading particles to make sure that
     * we don't find a local maximum.
     */
    assert(p.size()>=3);
    unsigned int n;
    p.size()==3 ? n=3 : n=4;
    vector<Vector3> tvec;
    vector<double> tval;
    std::sort(p.begin(), p.end(), mod2Cmp);
    for (unsigned int i=0 ; i<pow(2,n-1) ; i++) {
      // Create an initial vector from the leading four jets
      Vector3 foo(0,0,0);
      int sign=i;
      for (unsigned int k=0 ; k<n ; k++) {
        (sign%2)==1 ? foo+=p[k] : foo-=p[k];
        sign/=2;
      }
      foo=foo.unit();

      // Iterate
      double diff=999.;
      while (diff>1e-5) {
        Vector3 foobar(0,0,0);
        for (unsigned int k=0 ; k<p.size() ; k++)
          foo.dot(p[k])>0 ? foobar+=p[k] : foobar-=p[k];
        diff=(foo-foobar.unit()).mod();
        foo=foobar.unit();
      }

      // Calculate the thrust value for the vector we found
      t=0.;
      for (unsigned int k=0 ; k<p.size() ; k++)
        t+=fabs(foo.dot(p[k]));

      // Store everything
      tval.push_back(t);
      tvec.push_back(foo);
    }

    // Pick the solution with the largest thrust
    t=0.;
    for (unsigned int i=0 ; i<tvec.size() ; i++)
      if (tval[i]>t){
        t=tval[i];
        taxis=tvec[i];
      }
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
    getLog() << Log::DEBUG << "number of particles = " << threeMomenta.size() << endl;


    // Clear the caches
    _thrusts.clear();
    _thrustAxes.clear(); 


    // If there are fewer than 2 visible particles, we can't do much
    if (threeMomenta.size() < 2) {
      for (int i = 0; i < 3; ++i) {
        _thrusts.push_back(-1);
        _thrustAxes.push_back(Vector3(0,0,0));
      }
      return;
    }


    // Handle special case of thrust = 1 if there are only 2 particles
    if (threeMomenta.size() == 2) {
      Vector3 axis(0,0,0);
      _thrusts.push_back(1.0);
      _thrusts.push_back(0.0);
      _thrusts.push_back(0.0);
      axis = threeMomenta[0].unit();
      if (axis.z() < 0) axis = -axis;
      _thrustAxes.push_back(axis);
      /// @todo Improve this --- special directions bad...
      /// (a,b,c) _|_ 1/(a^2+b^2) (b,-a,0) etc., but which combination minimises error?
      if (axis.z() < 0.75)
        _thrustAxes.push_back( (axis.cross(Vector3(0,0,1))).unit() );
      else
        _thrustAxes.push_back( (axis.cross(Vector3(0,1,0))).unit() );
      _thrustAxes.push_back( _thrustAxes[0].cross(_thrustAxes[1]) );
      return;
    }



    // Temporary variables for calcs
    Vector3 axis(0,0,0); double val=0.;

    // Get thrust
    calcT(threeMomenta, val, axis);
    getLog() << Log::DEBUG << "Mom sum = " << momentumSum << endl;
    _thrusts.push_back(val / momentumSum);
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
    calcT(threeMomenta, val, axis);
    _thrusts.push_back(val / momentumSum);
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
      _thrustAxes.push_back(Vector3(0,0,0));
    }

  }


}
