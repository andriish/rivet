// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Cmp.hh"

namespace Rivet {

  int FastJets::compare(const Projection& p) const {
    const FastJets& other = dynamic_cast<const FastJets&>(p);
    return \
      pcmp(*_fsproj, *other._fsproj) || 
      cmp(_type, other._type) ||
      cmp(_recom, other._recom) ||
      cmp(_rparameter, other._rparameter);
  }


  void FastJets::project(const Event& e) {
    // Project into final state

    const FinalState& fs = e.applyProjection(*_fsproj);

    // Store 4 vector data about each particle into vecs
    vector<fastjet::PseudoJet> vecs;
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      HepMC::FourVector fv = p->getMomentum();
      // Store the FourVector in the KtLorentzVector form
      fastjet::PseudoJet psj(fv.px(),fv.py(),fv.pz(),fv.e());
      vecs.push_back(psj);
    }
    if (_cseq) delete _cseq;
    _cseq = new fastjet::ClusterSequence(vecs, _jdef);
  }


  vector<double> FastJets::getYSubJet(const fastjet::PseudoJet& jet) const {
    map<int,vector<double> >::iterator iter = _yscales.find(jet.cluster_hist_index());

    if (iter == _yscales.end()) {
      fastjet::ClusterSequence cseq(_cseq->constituents(jet),_jdef);
      vector<double> yMergeVals;
      for (int i=1; i<4; ++i) {
	// multiple the dmerge value by R^2 so that it corresponds to a
	// relative k_t (fastjet has 1/R^2 in the dij distance by default)
        yMergeVals.push_back(cseq.exclusive_dmerge(i)*_jdef.R()*_jdef.R()/jet.perp2());
        //yMergeVals.push_back(cseq.exclusive_dmerge(i)*_jdef.R()*_jdef.R());
      }
      _yscales.insert(make_pair( jet.cluster_hist_index(), yMergeVals ));
      return yMergeVals;
    } else {
      // This was cached.
      return iter->second;
    }

  }


}
