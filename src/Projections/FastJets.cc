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



  void FastJets::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const ParticleVector particles = fs.particles();
    _particles.clear();
    vector<PseudoJet> vecs;
  
    if (!particles.empty()) {
      // Store 4 vector data about each particle into vecs

      int counter = 1;
      foreach (const Particle& p, particles) {
        const FourMomentum fv = p.momentum();
        PseudoJet pJet(fv.px(), fv.py(), fv.pz(), fv.E());
        pJet.set_user_index(counter);
        vecs.push_back(pJet);
        _particles[counter] = p;
        ++counter;
      }
      
      getLog() << Log::DEBUG << "Running FastJet ClusterSequence construction" << endl;
      ClusterSequence cs(vecs, _jdef);
      _cseq = cs;
    } else {
      //ClusterSequence emptycs();
      //_cseq = emptycs;
      getLog() << Log::DEBUG << "No tracks... aborting" << endl;
      throw Error("Can't handle this situation, since FastJet's API doesn't allow us to pass null track vectors or set up a null ClusterSequence :(");
    }
  }



  Jets FastJets::_pseudojetsToJets(const PseudoJets& pjets) const {
    Jets rtn;
    foreach (const PseudoJet& pj, pjets) {
      Jet j;
      const PseudoJets parts = getClusterSeq().constituents(pj);
      foreach (const PseudoJet& p, parts) {
        map<int, Particle>::const_iterator found = _particles.find(p.user_index());
        if (found != _particles.end()) {
          // New way keeping full particle info
          j.addParticle(found->second);
        } else {
          // Old way storing just the momentum
          const FourMomentum particle(p.E(), p.px(), p.py(), p.pz());
          j.addParticle(particle);
        }
      }
      rtn.push_back(j);
    }
    return rtn;
  }



  vector<double> FastJets::getYSubJet(const fastjet::PseudoJet& jet) const {
    map<int,vector<double> >::iterator iter = _yscales.find(jet.cluster_hist_index());
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
  }
  


  fastjet::PseudoJet FastJets::splitJet(fastjet::PseudoJet jet, double& last_R) const { 
    // Sanity cuts
    if (jet.E() <= 0 || _cseq.constituents(jet).size() <= 1) {
      return jet; 
    }

    // Build a new cluster sequence just using the consituents of this jet.
    ClusterSequence cs(_cseq.constituents(jet), _jdef);

    // Get the jet back again
    fastjet::PseudoJet remadeJet = cs.inclusive_jets()[0];
    getLog() << Log::DEBUG << "Jet2:" << remadeJet.m() << "," << remadeJet.e() << endl;

    fastjet::PseudoJet parent1, parent2;
    fastjet::PseudoJet split(0.0, 0.0, 0.0, 0.0);
    while (cs.has_parents(remadeJet, parent1, parent2)) {
      getLog() << Log::DEBUG << "Parents:" << parent1.m() << "," << parent2.m() << endl; 
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
    split.reset(remadeJet.px(), remadeJet.py(), remadeJet.pz(), remadeJet.E());
    return split;
  }



  fastjet::PseudoJet FastJets::filterJet(fastjet::PseudoJet jet, double& stingy_R, const double def_R) const { 
    if (jet.E()<=0 || _cseq.constituents(jet).size()==0) { return jet; }
    if (stingy_R==0.) { stingy_R=def_R; }

    stingy_R = def_R < stingy_R ? def_R : stingy_R;
    fastjet::JetDefinition stingy_jet_def(fastjet::cambridge_algorithm, stingy_R);

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
