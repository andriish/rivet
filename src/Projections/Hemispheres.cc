// -*- C++ -*-
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {

  void Hemispheres::project(const Event& e) {
    // Get thrust axes.
    const AxesDefinition& ax = e.applyProjection(*_axproj);
    const Vector3 n = ax.axis1();

    FourMomentum p4With, p4Against;
    double Evis(0), broadWith(0), broadAgainst(0), broadDenom(0);
    const ParticleVector particles = e.applyProjection(_fsproj).particles();  
    for (ParticleVector::const_iterator p = particles.begin(); p != particles.end(); ++p) {
      const FourMomentum p4 = p->getMomentum();
      const Vector3 p3 = p4.vector3();
      const double p3Mag = mod(p3);
      const double p3Para = dot(p3, n);
      const double p3Trans = mod(p3 - p3Para * n);

      // Update normalisations
      Evis += p4.t();
      broadDenom += p3Mag;

      // Update the mass and broadening variables
      if (p3Para > 0) {
        p4With += p4;
        broadWith += p3Trans;
      } else if (p3Para < 0) {
        p4Against += p4;
        broadAgainst += p3Trans;
      } else {
        // In the incredibly unlikely event that a particle goes exactly along the
        // thrust plane, add half to each hemisphere.
        p4With += 0.5 * p4;
        p4Against += 0.5 * p4;
        broadWith += 0.5 * p3Trans;
        broadAgainst += 0.5 * p3Trans;
      }
    }

    // Visible energy squared.
    _E2vis = Evis * Evis;

    // Calculate masses.
    const double mass2With = mod2(p4With);
    const double mass2Against = mod2(p4Against);
    const bool withIsMaxMass2 = (mass2With > mass2Against);
    _M2high = (withIsMaxMass2) ? mass2With : mass2Against;
    _M2low = (withIsMaxMass2) ? mass2Against : mass2With;

    // Calculate broadenings.
    broadWith /= 2.0 * broadDenom;
    broadAgainst /= 2.0 * broadDenom;
    const bool withIsMaxBroad = (broadWith > broadAgainst);
    _Bmax = (withIsMaxBroad) ? broadWith : broadAgainst;
    _Bmin = (withIsMaxBroad) ? broadAgainst : broadWith;

    // Calculate high-max correlation flag.
    _highMassEqMaxBroad = (withIsMaxMass2 && withIsMaxBroad || !withIsMaxMass2 && !withIsMaxBroad);

  }

}
