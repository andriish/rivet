// -*- C++ -*-
#ifndef RIVET_RivetHepMC_HH
#define RIVET_RivetHepMC_HH

#ifdef ENABLE_HEPMC_3
#include "HepMC3/HepMC3.h"
#include "HepMC3/Relatives.h"
#include "HepMC3/Reader.h"

namespace Rivet{
  namespace RivetHepMC = HepMC3;
  using RivetHepMC::ConstGenParticlePtr;
  using RivetHepMC::ConstGenVertexPtr;
  using RivetHepMC::Relatives;
  using RivetHepMC::ConstGenHeavyIonPtr;
  
  using HepMC_IO_type = RivetHepMC::Reader;

  using PdfInfo = RivetHepMC::GenPdfInfo;
}

#else
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/Version.h"
#include "HepMC/GenRanges.h"
#include "HepMC/IO_GenEvent.h"


namespace Rivet{
  namespace RivetHepMC = HepMC;
  
  // HepMC 2.07 provides its own #defines
  #define ConstGenParticlePtr const HepMC::GenParticle*
  #define ConstGenVertexPtr const HepMC::GenVertex*
  #define ConstGenHeavyIonPtr const HepMC::HeavyIon*
  
  /// @brief Replicated the HepMC3 Relatives syntax using HepMC2 IteratorRanges
  /// This is necessary mainly because of capitalisation differences
  class Relatives{
    
    public:
    
    constexpr Relatives(HepMC::IteratorRange relo): _internal(relo){}
    
    constexpr HepMC::IteratorRange operator()() const {return _internal;}
    operator HepMC::IteratorRange() const {return _internal;}
    
    const static Relatives PARENTS;
    const static Relatives CHILDREN;
    const static Relatives ANCESTORS;
    const static Relatives DESCENDANTS;
    
    private:
    const HepMC::IteratorRange _internal;
    
  };
  
  using HepMC_IO_type = HepMC::IO_GenEvent;
  using PdfInfo = RivetHepMC::PdfInfo;

}
  
#endif

#include "Rivet/Tools/RivetSTL.hh"
#include "Rivet/Tools/Exceptions.hh"


namespace Rivet {

  using RivetHepMC::GenEvent;
  using ConstGenEventPtr = std::shared_ptr<const GenEvent>;
  /// @todo Use mcutils?

  namespace HepMCUtils{
  
    std::vector<ConstGenParticlePtr> particles(ConstGenEventPtr ge);
    std::vector<ConstGenParticlePtr> particles(const GenEvent *ge);
    std::vector<ConstGenVertexPtr>   vertices(ConstGenEventPtr ge);
    std::vector<ConstGenVertexPtr>   vertices(const GenEvent *ge);
    std::vector<ConstGenParticlePtr> particles(ConstGenVertexPtr gv, const Relatives &relo);
    std::vector<ConstGenParticlePtr> particles(ConstGenParticlePtr gp, const Relatives &relo);
    int uniqueId(ConstGenParticlePtr gp);
    int particles_size(ConstGenEventPtr ge);
    int particles_size(const GenEvent *ge);
    std::vector<ConstGenParticlePtr> beams(const GenEvent *ge);
    std::shared_ptr<HepMC_IO_type> makeReader(std::istream &istr);
    bool readEvent(std::shared_ptr<HepMC_IO_type> io, std::shared_ptr<GenEvent> evt);
  };
}

#endif
