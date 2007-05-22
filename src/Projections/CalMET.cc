// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CentralEtHCM class.
//

#include "Rivet/Projections/CalMET.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;
using namespace std;


int CalMET::compare(const Projection & p) const {
  const CalMET & other =
    dynamic_cast<const CalMET &>(p);
  return cmp(_met, other._met);
}


void CalMET::initialize(double etaMax, bool addMuons) { 
  _etaMax = etaMax;
  _addMuons = addMuons;
}


void CalMET::project(const Event & e) {
  const FinalState & fs = e.applyProjection(*_fs);
  _met = 0.;
  _metx = 0.;
  _mety = 0.;
  for ( int i = 0, N = fs.particles().size(); i < N; ++i ) {
    if (abs(fs.particles()[i].getPdgId()) != 12 && //no nu_e
	abs(fs.particles()[i].getPdgId()) != 14 && //no nu_mu
	abs(fs.particles()[i].getPdgId()) != 16 && //no nu_tau
	(abs(fs.particles()[i].getPdgId()) != 13 || _addMuons) //no mu
	&& fabs(fs.particles()[i].getMomentum().eta()) < _etaMax )
      _metx -= fs.particles()[i].getMomentum().px();
      _mety -= fs.particles()[i].getMomentum().py();
  }

  _met = sqrt(_metx * _metx  + _mety * _mety);

  return;
}

// RivetInfo CalMET::getInfo() const {
//   return Projection::getInfo() + _fs->getInfo();
// }

