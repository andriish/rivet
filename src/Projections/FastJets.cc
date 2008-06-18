// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"

using namespace fastjet;

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


  PseudoJet particleToPseudojet(const Particle& p) {
    const FourMomentum& fv = p.getMomentum();
    return PseudoJet(fv.px(), fv.py(), fv.pz(), fv.E());
  }


  void FastJets::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const ParticleVector particles = fs.particles();
    if (!particles.empty()) {
      // Store 4 vector data about each particle into vecs
      vector<PseudoJet> vecs(particles.size());
      transform(particles.begin(), particles.end(), vecs.begin(), particleToPseudojet);
      getLog() << Log::DEBUG << "Running FastJet ClusterSequence construction" << endl;
      ClusterSequence cs(vecs, _jdef);
      _cseq = cs;
    }
  }


  Jets FastJets::_pseudojetsToJets(const PseudoJets& pjets) const {
    Jets rtn;
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
      ClusterSequence subjet_cseq(_cseq.constituents(jet), _jdef);
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
  
  fastjet::PseudoJet FastJets::splitJet(fastjet::PseudoJet jet, double& last_R) 
    const { 

    cout << "Jet:" << jet.m() << "," << jet.e() << endl;
    
    fastjet::PseudoJet parent1, parent2;
    fastjet::PseudoJet split(0.0, 0.0, 0.0, 0.0);

    // Build a new cluster sequence just using the consituents of this jet.
    ClusterSequence cs(_cseq.constituents(jet), _jdef);

    // Get the jet back again
    fastjet::PseudoJet remadeJet = cs.inclusive_jets()[0];
    cout << "Jet2:" << remadeJet.m() << "," << remadeJet.e() << endl;

    while (cs.has_parents(remadeJet, parent1, parent2)) {
      cout << "Parents:" << parent1.m() << "," << parent2.m() << endl; 
      if (parent1.m2() < parent2.m2()) {
	fastjet::PseudoJet tmp;
	tmp = parent1; parent1 = parent2; parent2 = tmp;
      }
      
      double ktdist = parent1.kt_distance(parent2);
      double rtycut2 = 0.3*0.3;
      
      if (parent1.m() < ((2.0*remadeJet.m())/3.0) && ktdist > rtycut2*remadeJet.m2()) {
	break;
      } else {
	remadeJet = parent1;
      }
    }

    last_R = 0.5 * sqrt(parent1.squared_distance(parent2));    

    split.reset(parent1.px(), parent1.py(), parent1.pz(), parent1.E());

    return split;
  }


  fastjet::PseudoJet FastJets::filterJet(fastjet::PseudoJet jet, double& stingy_R) const { 

    stingy_R = 0.3 < stingy_R ? 0.3 : stingy_R;
    fastjet::JetDefinition stingy_jet_def(fastjet::cambridge_algorithm,
					  stingy_R);
    //FlavourRecombiner recom;
    //stingy_jet_def.set_recombiner(&recom);
    fastjet::ClusterSequence scs(_cseq.constituents(jet), stingy_jet_def);
    std::vector<fastjet::PseudoJet> stingy_jets = sorted_by_pt(scs.inclusive_jets());
    
    fastjet::PseudoJet reconst_jet(0.0, 0.0, 0.0, 0.0);
    
    for (unsigned isj = 0; isj < std::min(3U, (unsigned int)stingy_jets.size()); isj++) {
      reconst_jet += stingy_jets[isj];
    }	
    return reconst_jet;
  }
  
}
