// -*- C++ -*-

#include "Rivet/Rivet.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {

  void DISKinematics::project(const Event& e) {
    Log& log = getLog();

    const DISLepton& dislep = e.applyProjection(_lepton);
    const ParticlePair& inc = e.applyProjection(_beams).getBeams();
    Particle hadron;

    if ( inc.second.getPdgId() == _idhad ) hadron = inc.second;
    else if ( inc.first.getPdgId() == _idhad ) hadron = inc.first;
    else throw runtime_error("DISKinematics projector could not find the correct beam.");

    if ( &(dislep.in().getHepMCParticle()) == &(hadron.getHepMCParticle()) ) {
      throw runtime_error("DISKinematics projector could not find the correct beam.");
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
    const FourMomentum pGammaHCM =  _hcm.transform(pGamma);
    const FourMomentum pLepOutHCM =  _hcm.transform(pLepOut);
    /// @todo Does this commute?
    // Rotate DIS lepton on to \phi = 0
    _hcm.rotate(Vector3::Z(), -pLepOutHCM.azimuthalAngle());
    // Rotate photon on to z-axis
    /// @todo Is this putting the gamma on the -ve z-axis?
    _hcm.rotate(Vector3::Y(), -pGammaHCM.polarAngle());
    log << Log::DEBUG << "pGammaHCM = " << pGammaHCM << ", " << _hcm.transform(pGamma) << endl;
    log << Log::DEBUG << "pLepOutHCM = " << pLepOutHCM << ", " << _hcm.transform(pLepOut) << endl;
    // Boost to Breit frame
    _breit = _hcm;
    const double bz = 1 - 2*x();
    _breit.setBoost(Vector3::Z() * -bz);
  }


  int DISKinematics::compare(const Projection & p) const {
    const DISKinematics& other = dynamic_cast<const DISKinematics&>(p);
    return \
      pcmp(_lepton, other._lepton) || 
      pcmp(_beams, other._beams) || 
      cmp(_idhad, other._idhad);
  }


}
