// $Id: $
#ifndef RIVET_BEAMPARTICLE_H 
#define RIVET_BEAMPARTICLE_H 1

#include "Rivet/Rivet.hh"

namespace Rivet {

  /// Enumeration of available beam particles (using PDG IDs where available)
  enum BeamParticle { 
    ELECTRON = 11, 
    POSITRON = -11, 
    PROTON = 2212, 
    ANTIPROTON = -2212, 
    PHOTON = 22, 
    NEUTRON = 2112, 
    ANTINEUTRON = 2112, 
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
    TAU = 15, 
    ANTITAU = -15,
    EMINUS = 11, 
    EPLUS = -11, 
    P = 2212, 
    PBAR = -2212,
    GAMMA = 22,
    ANY = 10000,
    PHOTOELECTRON,
    PHOTOPOSITRON,
    PHOTOMUON,     
    PHOTOANTIMUON,
    PHOTOTAU,      
    PHOTOANTITAU
  };


  /// Typedef for a map of beam particle name enums to strings.
  typedef std::map<BeamParticle, std::string> BeamParticleMap;


  /// Typedef for a map of beam particle name strings to enums.
  typedef std::map<std::string, BeamParticle> BeamParticleMapR;


  /// Function which returns a map from beam particle enums to the corresponding name strings.
  inline BeamParticleMap getKnownBeamParticles() {
    BeamParticleMap bpmap;
    bpmap[ELECTRON] = "ELECTRON";
    bpmap[POSITRON] = "POSITRON";
    bpmap[PROTON] = "PROTON";
    bpmap[ANTIPROTON] = "ANTIPROTON";
    bpmap[PHOTON] = "PHOTON";
    bpmap[NEUTRON] = "NEUTRON";
    bpmap[ANTINEUTRON] = "ANTINEUTRON";
    bpmap[MUON] = "MUON";
    bpmap[ANTIMUON] = "ANTIMUON";
    bpmap[NU_E] = "NU_E";
    bpmap[NU_EBAR] = "NU_EBAR";
    bpmap[NU_MU] = "NU_MU";
    bpmap[NU_MUBAR] = "NU_MUBAR";
    bpmap[NU_TAU] = "NU_TAU";
    bpmap[NU_TAUBAR] = "NU_TAUBAR";
    bpmap[PIPLUS] = "PIPLUS"; 
    bpmap[PIMINUS] = "PIMINUS";
    bpmap[TAU] = "TAU"; 
    bpmap[ANTITAU] = "ANTITAU";
    bpmap[PHOTOELECTRON] = "PHOTOELECTRON";
    bpmap[PHOTOPOSITRON] = "PHOTOPOSITRON";
    bpmap[PHOTOMUON] = "PHOTOMUON";
    bpmap[PHOTOANTIMUON] = "PHOTOANTIMUON";
    bpmap[PHOTOTAU] = "PHOTOTAU"; 
    bpmap[PHOTOANTITAU] = "PHOTOANTITAU";
    return bpmap;
  }

  /// Function which returns a map from beam particle name strings to the corresponding enums.
  inline BeamParticleMapR getKnownBeamParticlesR() {
    BeamParticleMap bpmap = getKnownBeamParticles();
    BeamParticleMapR bpmapr;
    for (BeamParticleMap::const_iterator bp = bpmap.begin(); bp != bpmap.end(); ++bp) {
      bpmapr[bp->second] = bp->first;
    }
    return bpmapr;
  }


  /// Typedef for a collection of beam particle name enums.
  typedef std::vector<BeamParticle> BeamParticleList;


  /// Function which returns a vector of all the beam particle values in 
  /// the BeamParticle enum.
  inline BeamParticleList getKnownBeamParticleEnums() {
    BeamParticleList names;
    BeamParticleMap bpmap = getKnownBeamParticles();
    for (BeamParticleMap::const_iterator bp = bpmap.begin(); bp != bpmap.end(); ++bp) {
      names.push_back(bp->first);
    }
    return names;
  }


  /// Function which returns a vector of all the beam particle name strings.
  inline std::vector<std::string> getKnownBeamParticleNames() {
    std::vector<std::string> names;
    BeamParticleMap bpmap = getKnownBeamParticles();
    for (BeamParticleMap::const_iterator bp = bpmap.begin(); bp != bpmap.end(); ++bp) {
      names.push_back(bp->second);
    }
    return names;
  }
  
  
}


#endif // RIVET_BEAMPARTICLE_H
