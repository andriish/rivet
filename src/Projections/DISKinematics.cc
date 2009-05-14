// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Cmp.hh"
#include "Rivet/Tools/ParticleIDMethods.hh"

namespace Rivet {


  void DISKinematics::project(const Event& e) {
    const DISLepton& dislep = applyProjection<DISLepton>(e, "Lepton");
    const ParticlePair& inc = applyProjection<Beam>(e, "Beam").beams();

    bool firstIsHadron  = PID::isHadron(inc.first.pdgId());
    bool secondIsHadron = PID::isHadron(inc.second.pdgId());
    
    if (firstIsHadron && !secondIsHadron) {
      _inHadron = inc.first;
    } else if (!firstIsHadron && secondIsHadron) {
      _inHadron = inc.second;
    } else {
      //help!
      throw Error("DISKinematics projector could not find the correct beam hadron");
    }

    const FourMomentum pLepIn = dislep.in().momentum();
    const FourMomentum pLepOut = dislep.out().momentum();
    const FourMomentum pHad = _inHadron.momentum();
    const FourMomentum pGamma = pLepIn - pLepOut;
    const FourMomentum tothad = pGamma + pHad;
    _theQ2 = -pGamma.mass2();
    _theW2 = tothad.mass2();
    _theX = Q2()/(2.0 * pGamma * pHad);
    _theY = (pGamma * pHad) / (pLepIn * pHad);
    _theS = invariant(pLepIn + pHad);

    // Calculate boost vector for boost into HCM-system
    _hcm.setBoost(-tothad.boostVector());

    // Boost the gamma and lepton
    FourMomentum pGammaHCM = _hcm.transform(pGamma);

    // Rotate so the photon is in x-z plane
    _hcm.preMult(Matrix3(Vector3::mkZ(), -pGammaHCM.azimuthalAngle()));
    pGammaHCM = _hcm.transform(pGamma);
    assert(isZero(dot(pGammaHCM.vector3(), Vector3::mkY())));

    // Rotate so the photon is along the negative z-axis
    _hcm.preMult(Matrix3(Vector3::mkY(), pi - pGammaHCM.polarAngle()));

    // Check that final HCM photon lies along -ve z as expected
    pGammaHCM = _hcm.transform(pGamma);
    assert(isZero(dot(pGammaHCM.vector3(), Vector3::mkX())) &&
           isZero(dot(pGammaHCM.vector3(), Vector3::mkY())));
    assert(isZero(angle(pGammaHCM.vector3(), -Vector3::mkZ())));

    // Finally rotate so outgoing lepton at phi = 0
    FourMomentum pLepOutHCM = _hcm.transform(pLepOut);
    _hcm.preMult(Matrix3(Vector3::mkZ(), -pLepOutHCM.azimuthalAngle()));
    pLepOutHCM = _hcm.transform(pLepOut);
    assert(isZero(pLepOutHCM.azimuthalAngle()));

    // Boost to Breit frame    
    const double bz = 1 - 2*x();
    _breit = LorentzTransform(Vector3::mkZ() * bz).combine(_hcm);
  }


  int DISKinematics::compare(const Projection & p) const {
    const DISKinematics& other = pcast<DISKinematics>(p);
    return mkNamedPCmp(other, "Lepton"); 
  }


}
