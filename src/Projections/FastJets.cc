// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Cmp.hh"


namespace Rivet {


  int FastJets::compare(const Projection& p) const {
    const FastJets& other = dynamic_cast<const FastJets&>(p);
    return \
      mkNamedPCmp(other, "FS") || 
      cmp(_jdef.jet_algorithm(), other._jdef.jet_algorithm()) ||
      cmp(_jdef.recombination_scheme(), other._jdef.recombination_scheme()) ||
      cmp(_jdef.plugin(), other._jdef.plugin()) ||
      cmp(_jdef.R(), other._jdef.R());
  }


  void FastJets::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    // Store 4 vector data about each particle into vecs
    PseudoJets vecs;
    for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
      HepMC::FourVector fv = p->getMomentum();
      fastjet::PseudoJet psj(fv.px(),fv.py(),fv.pz(),fv.e());
      vecs.push_back(psj);
    }
    fastjet::ClusterSequence cs(vecs, _jdef);
    _cseq = cs;
  }


  Jets FastJets::getJets() const {
    Jets rtn;
    const PseudoJets pjets = _cseq.inclusive_jets();
    for (PseudoJets::const_iterator pj = pjets.begin(); pj != pjets.end(); ++pj) {
      Jet j;
      const PseudoJets parts = getClusterSeq().constituents(*pj);
      for (PseudoJets::const_iterator p = parts.begin(); p != parts.end(); ++p) {
        const FourMomentum particle(p->E(), p->px(), p->py(), p->pz());
        j.addParticle(particle);
      }
      rtn.push_back(j);
    }
    return rtn;
  }


  vector<double> FastJets::getYSubJet(const fastjet::PseudoJet& jet) const {
    map<int,vector<double> >::iterator iter = _yscales.find(jet.cluster_hist_index());
    if (iter == _yscales.end()) {
      fastjet::ClusterSequence subjet_cseq(_cseq.constituents(jet), _jdef);
      vector<double> yMergeVals;
      for (int i = 1; i < 4; ++i) {
        // Multiply the dmerge value by R^2 so that it corresponds to a
        // relative k_T (fastjet has 1/R^2 in the d_ij distance by default)
        const double ktmerge = subjet_cseq.exclusive_dmerge(i) * _jdef.R()*_jdef.R();
        yMergeVals.push_back(ktmerge/jet.perp2());
      }
      _yscales.insert(make_pair( jet.cluster_hist_index(), yMergeVals ));
      return yMergeVals;
    } else {
      // This was cached.
      return iter->second;
    }
    
  }
  
  
}
