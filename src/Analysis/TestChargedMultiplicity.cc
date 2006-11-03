// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TestMultiplicity class.
//

#include "Rivet/Analysis/TestChargedMultiplicity.h"
#include "HepMC/ParticleDataTable.h"
#include "HepMC/ParticleData.h"

using namespace Rivet;
using namespace std;
using namespace HepMC;

TestChargedMultiplicity::~TestChargedMultiplicity() {}

void TestChargedMultiplicity::init() {
  // Book histogram here.
}

void TestChargedMultiplicity::analyze(const Event & event) {
  const FinalStateProjection& fs = event.addProjection(fsproj);
  unsigned int chmult(0), unchmult(0);
  
  unsigned int particleNum(0);
  for (PVector::const_iterator p = fs.particles().begin(); p != fs.particles().end(); ++p) {
    particleNum++;
    int pdgCode = p->id;
    ParticleDataTable pdgTable;
    ParticleData* pdgData = pdgTable.find(pdgCode);
    if (pdgData) {
      if (pdgData->is_hadron()) {
        if (pdgData->charge() != 0) {
          ++chmult;
          cout << "Incrementing charged multiplicity = " << chmult << endl;
        } else {
          ++unchmult;
          cout << "Incrementing uncharged multiplicity = " << unchmult << endl;
        }
      }
    } else {
      cerr << "Null pointer to PDG data for particle " << particleNum << endl;
      /// @todo Make a default string rep for a Particle
    }
  }
  cout << "Event charged multiplicity = " << chmult << endl;
  cout << "Event uncharged multiplicity = " << unchmult << endl;
  // Fill histogram here.
}

void TestChargedMultiplicity::finalize() {}

RivetInfo TestChargedMultiplicity::getInfo() const {
  return AnalysisBase::getInfo() + fsproj.getInfo();
}
