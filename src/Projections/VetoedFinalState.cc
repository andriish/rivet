// -*- C++ -*-

#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/Cmp.hh"
#include <algorithm>

using namespace std;

namespace Rivet {


  int VetoedFinalState::compare(const Projection& p) const {
    const VetoedFinalState& other = dynamic_cast<const VetoedFinalState&>(p);
    const int fscmp = pcmp(*_fsproj, *other._fsproj);
    if (fscmp != 0) return fscmp;
    if (_vetoCodes == other._vetoCodes) return 0;
    if (_vetoCodes < other._vetoCodes) return -1; else return 1;
  }


  void VetoedFinalState::project(const Event& e) {
    Log log = getLog();
    const FinalState& fs = e.applyProjection(*_fsproj);
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());
    const ParticleVector& fsps = fs.particles();
    for (ParticleVector::const_iterator p = fsps.begin(); p != fsps.end(); ++p) {
      if (log.isActive(Log::DEBUG)) {
        stringstream codes; 
        codes << "{";
        for (map<long,vector<double> >::const_iterator code = _vetoCodes.begin(); 
             code != _vetoCodes.end(); ++code) {
	  codes << " " << code->first;
	}
        codes << " }";
        log << Log::DEBUG << p->getPdgId() << " vs. veto codes = " 
            << codes.str() << " (" << _vetoCodes.size() << ")" << endl;
      }
      const long pdgid = p->getPdgId();
      map<long,vector<double> >::iterator iter = _vetoCodes.find(pdgid);
      if ( (iter == _vetoCodes.end())) {
	log << Log::DEBUG << "Storing with PDG code " << pdgid << " pt " << p->getMomentum().perp() << endl;
	_theParticles.push_back(*p);
      } else {
	vector<double> range = iter->second;
	double pt = p->getMomentum().perp();
	if (range.size()>1) log << Log::DEBUG << "ID " << pdgid << " pt range " << range.at(0) << "," << range.at(1) << endl;
	if ( range.size()>1 &&
	     (pt < range.at(0) || (pt>range.at(1) && range.at(1)>0.0))) {
	  log << Log::DEBUG << "Storing with PDG code " << pdgid << " pt " << p->getMomentum().perp() << endl;
	  _theParticles.push_back(*p);
	} else {
	  log << Log::DEBUG << "Vetoed with PDG code " << pdgid << " pt " << p->getMomentum().perp()<< endl;
	}
      }
    }
  } 

}


