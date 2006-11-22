// -*- C++ -*-

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedMultiplicity.hh"
#include "Rivet/Projections/Cmp.hh"
#include "HepPDT/ParticleID.hh"

using namespace Rivet;

ChargedMultiplicity::~ChargedMultiplicity() {}

int ChargedMultiplicity::compare(const Projection & p) const {
  const ChargedMultiplicity & other = dynamic_cast<const ChargedMultiplicity &>(p);
  return cmp(etamin, other.etamin) || cmp(etamax, other.etamax);
}

void ChargedMultiplicity::project(const Event & e) {
  Logger& log = getLogger();
  log.setPriority(LogPriority::INFO);

  /// @todo Clear something
  const FinalStateProjection& fs = event.addProjection(fsproj);
  unsigned int chmult(0), unchmult(0);
  unsigned int particleNum(0);
  for (ParticleVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    ++particleNum;
    HepPDT::ParticleID pInfo(p->id);
    bool isHadron = pInfo.isHadron();
    if (pInfo.threeCharge() != 0) {
      ++chmult;
      if (isHadron) ++hchmult;
      log << LogPriority::DEBUG << "Incrementing charged multiplicity   = " << chmult 
          << " (" << hchmult << " hadrons)" << endlog;
    } else {
      ++unchmult;
      if (isHadron) ++hunchmult;
      log << LogPriority::DEBUG << "Incrementing uncharged multiplicity = " << unchmult 
          << " (" << hunchmult << " hadrons)" << endlog;
    }
  }
  log << LogPriority::INFO << "Event charged multiplicity   = " << chmult << endlog;
  log << LogPriority::INFO << "Event uncharged multiplicity = " << unchmult << endlog;



  for ( GenEvent::particle_const_iterator pi = e.genEvent().particles_begin();
        pi != e.genEvent().particles_end(); ++pi ) {
    if ( (*pi)->status() == 1 &&
         (*pi)->momentum().eta() > etamin &&
         (*pi)->momentum().eta() < etamax)
      theParticles.push_back(Particle(**pi));
  }
}

