#ifndef RIVET_PARTICLENAME_HH
#define RIVET_PARTICLENAME_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Particle.fhh"

namespace Rivet {


  /// Enumeration of available beam particles (using PDG IDs where available)
  enum ParticleName {
    ELECTRON = 11,
    POSITRON = -11,
    PROTON = 2212,
    ANTIPROTON = -2212,
    PHOTON = 22,
    NEUTRON = 2112,
    ANTINEUTRON = -2112,
    MUON = 13,
    ANTIMUON = -13,
    NU_E = 12,
    NU_EBAR = -12,
    NU_MU = 14,
    NU_MUBAR = -14,
    NU_TAU = 16,
    NU_TAUBAR = -16,
    PIPLUS = 211,
    PIMINUS = -211,
    K0L = 130,
    K0S = 310,
    KPLUS = 321,
    KMINUS = -321,
    LAMBDA = 3122,
    LAMBDABAR = -3122,
    XIMINUS = 3312,
    XIPLUS = -3312,
    OMEGAMINUS = 3334,
    OMEGAPLUS = -3334,
    TAU = 15,
    ANTITAU = -15,
    EMINUS = 11,
    EPLUS = -11,
    P = 2212,
    PBAR = -2212,
    GLUON = 21,
    GAMMA = 22,
    WPLUSBOSON = 24,
    WMINUSBOSON = -24,
    ZBOSON = 23,
    HIGGS = 25,
    DQUARK = 1,
    UQUARK = 2,
    SQUARK = 3,
    CQUARK = 4,
    BQUARK = 5,
    TQUARK = 6,
    ANY = 10000,
    PHOTOELECTRON,
    PHOTOPOSITRON,
    PHOTOMUON,
    PHOTOANTIMUON,
    PHOTOTAU,
    PHOTOANTITAU
  };


  /// Convenience maker of particle ID pairs.
  inline std::pair<PdgId,PdgId> make_pdgid_pair(PdgId a, PdgId b) {
    return make_pair<PdgId,PdgId>(a, b);
  }

  /// Convenience maker of particle ID pairs.
  // inline std::pair<PdgId,PdgId> make_pdgid_pair(ParticleName aname, ParticleName bname) {
  //   return make_pdgid_pair(aname, bname);
  // }

  /// Convenience maker of particle ID pairs.
  inline std::pair<PdgId,PdgId> make_pdgid_pair(const std::pair<ParticleName,ParticleName>& pnamepair) {
    return make_pdgid_pair(pnamepair.first, pnamepair.second);
  }

  /// Typedef for a map of beam particle name enums to strings.
  typedef std::map<PdgId, std::string> ParticleNameMap;


  /// Typedef for a map of beam particle name strings to enums.
  typedef std::map<std::string, PdgId> ParticleNameMapR;


  static ParticleNameMap s_pnames;


  /// Function which returns a map from beam particle name strings to the corresponding enums.
  inline ParticleNameMapR getParticleNamesRMap() {
    ParticleNameMapR bpmapr;
    for (ParticleNameMap::const_iterator bp = s_pnames.begin(); bp != s_pnames.end(); ++bp) {
      bpmapr[bp->second] = bp->first;
    }
    return bpmapr;
  }


  /// Typedef for a collection of beam particle name enums.
  typedef std::vector<PdgId> ParticleNameList;


  /// Function which returns a vector of all the beam particle values in
  /// the ParticleName enum.
  inline ParticleNameList getParticleNameEnums() {
    ParticleNameList names;
    for (ParticleNameMap::const_iterator bp = s_pnames.begin(); bp != s_pnames.end(); ++bp) {
      names.push_back(bp->first);
    }
    return names;
  }


  /// Function which converts a particle name string to a ParticleName enum
  inline ParticleName getParticleNameEnum(const std::string& pname) {
    return (ParticleName) Rivet::getParticleNamesRMap()[pname];
  }



  /// Function which returns a vector of all the beam particle name strings.
  inline std::vector<std::string> getParticleNames() {
    vector<string> names;
    for (ParticleNameMap::const_iterator bp = s_pnames.begin(); bp != s_pnames.end(); ++bp) {
      names.push_back(bp->second);
    }
    return names;
  }


  /// Print a ParticleName as a string.
  inline std::string toString(const ParticleName& p) {
    return s_pnames[p];
  }


  /// Print a PdgId as a named string.
  inline std::string toParticleName(PdgId p) {
    if (s_pnames.find(p) != s_pnames.end()) {
      return s_pnames[p];
    }
    ostringstream ss;
    ss << p;
    return ss.str();
  }


  /// Allow ParticleName to be passed to an iostream.
  inline std::ostream& operator<<(std::ostream& os, const ParticleName& p) {
    os << toString(p);
    return os;
  }


  /////////////////////////////////////////////////
  // Beams

  /// Print a BeamPair as a string.
  inline std::string toString(const BeamPair& pair) {
    string out = "[" +
      toParticleName(pair.first) + ", " +
      toParticleName(pair.second) + "]";
    return out;
  }

  /// Allow BeamPair to be passed to an ostream.
  inline std::ostream& operator<<(std::ostream& os, const BeamPair& bp) {
    os << toString(bp);
    return os;
  }


}


#endif
