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
    const FourMomentum& fv = p.momentum();
    return PseudoJet(fv.px(), fv.py(), fv.pz(), fv.E());
  }

  void FastJets::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    const ParticleVector particles = fs.particles();
  
    if (!particles.empty()) {
      // Store 4 vector data about each particle into vecs
      vector<PseudoJet> vecs(particles.size());
      transform(particles.begin(), particles.end(), vecs.begin(), particleToPseudojet);
      _particles.clear();
      _particles.insert(particles.begin(), particles.end());
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

        //older way adding just the particle momentum to the jet
        //j.addParticle(particle);
        
        // Look up the particle by its PT to get ID information

        Particle toLookup;
        toLookup.setMomentum(particle);
        set<Particle>::const_iterator high = _particles.upper_bound(toLookup);
                
        if(high != _particles.begin()){
          if(high==_particles.end()){
            //check this particle has a PT within 1% of the hardest particle
            set<Particle>::const_iterator end=_particles.end();
            if(!_particles.empty()){
              --end;
              double maxPT2 = end->momentum().pT2();
              if(fuzzyEquals(maxPT2, particle.pT2(), 0.01 * maxPT2)){
                --high;
              }else{
                // this is bad
                if(!_particles.empty()) throw Error("particle to lookup has PT above the max input to the jet finder!");
              }
            }
          }else{
            --high;
          }
        }else{
          if(!_particles.empty()) throw Error("FastJets failed to lookup particle information!");
        }
                
        getLog() << Log::TRACE << " Particle in jet Vs. full particle record = " 
        <<particle.px()<<", "<<particle.py()<<", "<<particle.pz()<<", "<<particle.E()<<" Vs. "
        <<high->momentum().px()<<", "<< high->momentum().py()<<", "<< high->momentum().pz()
        <<", "<<high->momentum().E()<<endl;
        
        
        //add the complete particle information to the jet
        j.addParticle(*high);
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
  
  fastjet::PseudoJet FastJets::splitJet(fastjet::PseudoJet jet, double& last_R)
    const { 

    if (jet.E()<=0 || _cseq.constituents(jet).size()<=1) { return jet; }

    fastjet::PseudoJet parent1, parent2;
    fastjet::PseudoJet split(0.0, 0.0, 0.0, 0.0);

    // Build a new cluster sequence just using the consituents of this jet.
    
    ClusterSequence cs(_cseq.constituents(jet), _jdef);

    // Get the jet back again
    fastjet::PseudoJet remadeJet = cs.inclusive_jets()[0];
    //cout << "Jet2:" << remadeJet.m() << "," << remadeJet.e() << endl;

    while (cs.has_parents(remadeJet, parent1, parent2)) {
      //cout << "Parents:" << parent1.m() << "," << parent2.m() << endl; 
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
