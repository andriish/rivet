// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/KtJets.hh"

using namespace Rivet;

KtJets::~KtJets() {}

/* The compare method is incomplete? */

int KtJets::compare(const Projection & p) const {
  //const KtJets & other = dynamic_cast<const KtJets &>(p);
  return 0;
}

void KtJets::project(const Event & e) {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);

  vector<KtJet::KtLorentzVector> vecs;

  type_ = (*this).type_;
      
  // Project into final state
  const FinalState& fs = e.applyProjection(fsproj);
  
  // Store 4 vector data about each particle into vecs
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    HepMC::FourVector fv = p->momentum;
    // store the FourVector in the KtLorentzVector form
    KtJet::KtLorentzVector ktlv(fv.px(), fv.py(), fv.pz(), fv.e());
    vecs.push_back(ktlv);
  }
  vecs_ = vecs;
  type_ = this->type_;
  angle_ = this->angle_;
  recom_ = this-> recom_;
  rparameter_ = this-> rparameter_;
}
