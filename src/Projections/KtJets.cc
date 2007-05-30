// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/KtJets.hh"
#include "Rivet/Projections/Cmp.hh"

using namespace Rivet;

int KtJets::compare(const Projection & p) const {
  const KtJets & other = dynamic_cast<const KtJets &>(p);
  return pcmp(*fsproj, *other.fsproj) || cmp(type_, other.type_) ||
    cmp(angle_, other.angle_) || cmp(recom_, other.recom_) ||
    cmp(rparameter_, other.rparameter_);
}


void KtJets::project(const Event & e) {
  vector<KtJet::KtLorentzVector> vecs;

  // Project into final state
  const FinalState& fs = e.applyProjection(*fsproj);
  
  // Store 4 vector data about each particle into vecs
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    HepMC::FourVector fv = p->getMomentum();
    // store the FourVector in the KtLorentzVector form
    KtJet::KtLorentzVector ktlv(fv.px(), fv.py(), fv.pz(), fv.e());
    vecs.push_back(ktlv);
  }
  if ( pktev_ ) delete pktev_;

  pktev_ = new KtJet::KtEvent(vecs, type_, angle_, recom_, rparameter_);

}

vector<double> KtJets::getYSubJet(const KtJet::KtLorentzVector &jet) const {

  map<int,vector<double> >::iterator iter = yscales_.find(jet.getID());

  if( iter == yscales_.end() ) {
    
    KtJet::KtEvent subj = KtJet::KtEvent(jet, angle_, recom_);
    vector<double> yMergeVals = vector<double>();
    for(int i=1; i<4; ++i) {
      yMergeVals.push_back(subj.getYMerge(i));
    }
    yscales_.insert(make_pair( jet.getID(), yMergeVals ));
    return yMergeVals;

  } else {
    // This was cached.
    return iter->second;

  }
}
