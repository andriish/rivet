// $Id: $
#ifndef RIVET_RIVETCLHEP_HH 
#define RIVET_RIVETCLHEP_HH 1

#define CLHEP_MAX_MIN_DEFINED
#define CLHEP_SQR_DEFINED
#define CLHEP_ABS_DEFINED

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Vector/LorentzRotation.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace CLHEP {
  typedef HepLorentzRotation LorentzRotation;
  typedef HepLorentzVector LorentzVector;
  typedef HepRotation Rotation;
  typedef Hep3Vector Vector3;
}

namespace Rivet {
  using CLHEP::LorentzRotation;
  using CLHEP::LorentzVector;
  using CLHEP::Rotation;
  using CLHEP::Vector3;
}

#endif // RIVET_RIVETCLHEP_HH
