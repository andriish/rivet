#include "Rivet/Particle.hh"
#include "Rivet/Projections/Beam.hh"

int main() {
  using namespace Rivet;
  using namespace std;

  const FourMomentum pbeam1a(1*GeV, 0, 0, 1*TeV), pbeam1b(1*GeV, 0, 0, -1*TeV);
  const Particle beam1a(PID::PROTON, pbeam1a), beam1b(PID::PROTON, pbeam1b);
  //
  const Vector3 beta1a = cmsBoostBetaVec(pbeam1a, pbeam1b);
  const Vector3 beta1b = cmsBoostBetaVec(beam1a, beam1b);
  cout << "Beta_symm = " << beta1a << " or " << beta1b << endl;
  //
  const Vector3 gamma1a = cmsBoostGammaVec(pbeam1a, pbeam1b);
  const Vector3 gamma1b = cmsBoostGammaVec(beam1a, beam1b);
  cout << "Gamma_symm = " << gamma1a << " or " << gamma1b << endl;

  cout << endl;

  const FourMomentum pbeam2a(1*GeV, 0, 0, 10*GeV), pbeam2b(2*GeV, 0, 0, 0);
  // const Particle beam2a(, pbeam1), beam2b(, pbeam1);
  //
  const Vector3 beta2a = cmsBoostBetaVec(pbeam2a, pbeam2b);
  //const Vector3 beta2b = cmsBoostBetaVec(beam2a, beam2b);
  cout << "Beta_asymm = " << beta2a << endl; // << " or " << beta2b << endl;
  //
  const Vector3 gamma2a = cmsBoostGammaVec(pbeam2a, pbeam2b);
  //const Vector3 gamma2b = cmsBoostGammaVec(beam2a, beam2b);
  cout << "Gamma_asymm = " << gamma2a << endl; // << " or " << gamma2b << endl;

  const LorentzTransform trfb = LorentzTransform::mkFrameTransformFromBeta(cmsBoostBetaVec(pbeam1a, pbeam1b));
  const FourMomentum pbeam2c = trfb.transform(pbeam2a);
  const FourMomentum pbeam2d = trfb.transform(pbeam2b);
  cout << "Beta-boosted_asymm = " << pbeam2c << " + " << pbeam2d << " = " << (pbeam2c + pbeam2d) << endl;

  const LorentzTransform trfy = LorentzTransform::mkFrameTransformFromGamma(cmsBoostGammaVec(pbeam1a, pbeam1b));
  const FourMomentum pbeam2e = trfy.transform(pbeam2a);
  const FourMomentum pbeam2f = trfy.transform(pbeam2b);
  cout << "Gamma-boosted_asymm = " << pbeam2e << " + " << pbeam2f << " = " << (pbeam2e + pbeam2f) << endl;

  return 0;
}
