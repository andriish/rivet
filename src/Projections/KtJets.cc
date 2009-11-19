// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/KtJets.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {


  int KtJets::compare(const Projection& p) const {
    const KtJets& other = dynamic_cast<const KtJets&>(p);
    return \
      mkNamedPCmp(other "FS") ||
      cmp(_type, other._type) ||
      cmp(_angle, other._angle) ||
      cmp(_recom, other._recom) ||
      cmp(_rparameter, other._rparameter);
  }


  void KtJets::project(const Event& e) {
    // Project into final state
    const FinalState& fs = applyProjection<FinalState>(e, "FS");

    // Store 4 vector data about each particle into vecs
    vector<KtJet::KtLorentzVector> vecs;
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      HepMC::FourVector fv = p->momentum();
      // Store the FourVector in the KtLorentzVector form
      KtJet::KtLorentzVector ktlv(fv.px(), fv.py(), fv.pz(), fv.e());
      vecs.push_back(ktlv);
    }
    if (_pktev) delete _pktev;
    _pktev = new KtJet::KtEvent(vecs, _type, _angle, _recom, _rparameter);
  }


  vector<double> KtJets::getYSubJet(const KtJet::KtLorentzVector& jet) const {
    map<int,vector<double> >::iterator iter = _yscales.find(jet.getID());
    if (iter == _yscales.end()) {
      KtJet::KtEvent subj = KtJet::KtEvent(jet, _angle, _recom);
      vector<double> yMergeVals;
      for (int i=1; i<5; ++i) {
        if (subj.getNConstituents() > i){
          yMergeVals.push_back(subj.getYMerge(i));
        }
      }
      _yscales.insert(make_pair( jet.getID(), yMergeVals ));
      return yMergeVals;
    } else {
      // This was cached.
      return iter->second;
    }

  }


}
