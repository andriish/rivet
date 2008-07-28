// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {

  void DISKinematics::project(const Event& e) {
    //Log& log = getLog();

    const DISLepton& dislep = applyProjection<DISLepton>(e, "Lepton");
    const ParticlePair& inc = applyProjection<Beam>(e, "Beam").getBeams();
    Particle hadron;

    if ( inc.second.getPdgId() == _idhad ) hadron = inc.second;
    else if ( inc.first.getPdgId() == _idhad ) hadron = inc.first;
    else throw Error("DISKinematics projector could not find the correct beam.");

    if ( &(dislep.in().getHepMCParticle()) == &(hadron.getHepMCParticle()) ) {
      throw Error("DISKinematics projector could not find the correct beam.");
    }

    const FourMomentum pLepIn = dislep.in().getMomentum();
    const FourMomentum pLepOut = dislep.out().getMomentum();
    const FourMomentum pHad = hadron.getMomentum();
    const FourMomentum pGamma = pLepIn - pLepOut;
    const FourMomentum tothad = pGamma + pHad;
    _theQ2 = -pGamma.mass2();
    _theW2 = tothad.mass2();
    _theX = Q2()/(2.0 * pGamma * pHad);
    _theY = (pGamma * pHad) / (pLepIn * pHad);
    _theS = invariant(pLepIn + pHad);

    // Calculate boost vector for boost into HCM-system
    _hcm.setBoost(-boostVector(tothad));
    // Boost the gamma and lepton
    FourMomentum pGammaHCM =  _hcm.transform(pGamma);
    // rotate so photon in x-z plane
    _hcm.preMult(Matrix3(Vector3::mkZ(), -pGammaHCM.azimuthalAngle()));
    // rotate so photon along z axis
    pGammaHCM = _hcm.transform(pGamma);
    _hcm.preMult(Matrix3(Vector3::mkY(), -pGammaHCM.polarAngle()));
    // rotate by 180 if along -z 
    pGammaHCM = _hcm.transform(pGamma);
    if(pGammaHCM.z() > 0.0) _hcm.preMult(Matrix3(Vector3::mkY(), pi));
    // finally rotate so outgoing lepton at phi=0
    FourMomentum pLepOutHCM =  _hcm.transform(pLepOut);
    _hcm.preMult(Matrix3(Vector3::mkZ(), -pLepOutHCM.azimuthalAngle()));
    // Boost to Breit frame
    const double bz = 1 - 2*x();
    _breit = LorentzTransform(Vector3::mkZ() * bz).combine(_hcm);
  }


  int DISKinematics::compare(const Projection & p) const {
    const DISKinematics& other = pcast<DISKinematics>(p);
    return \
      mkNamedPCmp(other, "Lepton") || 
      cmp(_idhad, other._idhad);
  }


}
