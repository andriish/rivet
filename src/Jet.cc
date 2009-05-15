#include "Rivet/Jet.hh"

namespace Rivet {

  Jet::Jet() 
    : ParticleBase()
  {
    clear();
  }


  Jet& Jet::setParticles(vector<FourMomentum> particles) {
    _particles = particles;
    _resetCaches();
    return *this;
  }


  Jet& Jet::addParticle(FourMomentum particle) {
    _particles.push_back(particle);
    _resetCaches();
    return *this;
  }


  Jet& Jet::addParticle(const Particle& particle) {
    _fullParticles.push_back(particle);
    _particles.push_back(particle.momentum());
    _resetCaches();
    return *this;
  }
    

  bool Jet::containsParticle(const Particle& particle) const {
    const int barcode = particle.genParticle().barcode();
    foreach (const Particle& p, _fullParticles) {
      const GenParticle& part = p.genParticle();
      /// @todo Prefer to compare GenEvent pointers?
      if (part.barcode() == barcode) return true;
    }
    return false;
  }


  bool Jet::containsParticleId(PdgId pid) const {
    foreach (const Particle& p, _fullParticles) {
      if (p.pdgId() == pid) return true;
    }
    return false;
  }


  Jet& Jet::clear() {
    _particles.clear();
    _fullParticles.clear();
    _resetCaches();
    return *this;
  }


  double Jet::ptWeightedEta() const {
    _calcPtAvgs();
    assert(_okPtWeightedEta);
    return _ptWeightedEta;
  }


  double Jet::ptWeightedPhi() const {
    _calcPtAvgs();
    assert(_okPtWeightedPhi);
    return _ptWeightedPhi;
  }


  double Jet::eta() const {
    _calcAvgs();
    assert(_okEta);
    return _eta;
  }
  

  double Jet::phi() const {
    _calcAvgs();
    assert(_okPhi);
    return _phi;
  }


  const FourMomentum& Jet::momentum() const {
    _calcMomVector();
    return _momentum;
  }

    
  FourMomentum& Jet::momentum() { 
    _calcMomVector();
    return _momentum;
  }

    
  double Jet::ptSum() const {
    if (!_okTotalPt) {
      double ptsum(0.0);
      for (const_iterator p = this->begin(); p != this->end(); ++p) {
        ptsum += p->pT();
      }
      _totalPt = ptsum;
      _okTotalPt = true;
    }
    return _totalPt;
  }


  double Jet::EtSum() const {
    if (!_okTotalEt) {
      double Etsum(0.0);
      for (const_iterator p = this->begin(); p != this->end(); ++p) {
        Etsum += p->Et();
      }
      _totalEt = Etsum;
      _okTotalEt = true;
    }
    return _totalEt;
  }


  void Jet::_resetCaches() const {
    _okPhi = false;
    _okEta = false;
    _okPtWeightedPhi = false;
    _okPtWeightedEta = false;
    _okTotalPt = false;
    _okTotalEt = false;
    _okMomentum = false;
  }
  
  
  void Jet::_calcMomVector() const {
    if (!_okMomentum) {
      _momentum = accumulate(begin(), end(), FourMomentum());
      _okMomentum = true;
    }
  }


  void Jet::_calcPtAvgs() const {
    if (!_okPtWeightedEta || !_okPtWeightedPhi) {
      double ptwetasum(0.0), ptwphisum(0.0), ptsum(0.0);
      double phibegin = 0.0;
      for (const_iterator p = this->begin(); p != this->end(); ++p) {
        double pt = p->pT();
        ptsum += pt;
        ptwetasum += pt * p->pseudorapidity();
        
        if (p == this->begin()) {
          phibegin = p->azimuthalAngle();
        } else {
          const double dphi = p->azimuthalAngle() - phibegin;
          ptwphisum += pt * mapAngleMPiToPi(dphi);
        }
      }
      _totalPt = ptsum;
      _okTotalPt = true;
      _ptWeightedEta = ptwetasum / ptSum();
      _okPtWeightedEta = true;
      _ptWeightedPhi = phibegin + ptwphisum / ptSum();
      _ptWeightedPhi = mapAngleMPiToPi(_ptWeightedPhi);
      _okPtWeightedPhi = true;
    }
  }
  

  void Jet::_calcAvgs() const {
    if (!_okEta || !_okPhi) {
      double etasum(0.0), phisum(0.0);
      double phibegin = 0.0;
      for (const_iterator p = this->begin(); p != this->end(); ++p) {
        etasum += p->pseudorapidity();
        if (p == this->begin()) {
          phibegin = p->azimuthalAngle();
        } else {
          const double dphi = p->azimuthalAngle() - phibegin;
          phisum += mapAngleMPiToPi(dphi);
        }
      }
      const double dnum = _particles.size();
      _eta = etasum / dnum;
      _okEta = true;
      _phi = phibegin + phisum / dnum;
      _phi = mapAngleMPiToPi(_phi);
      _okPhi = true;
    }  
  }
  
}
