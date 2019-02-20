// -*- C++ -*-

#include "Rivet/Tools/RivetHepMC.hh"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/ReaderAsciiHepMC2.h"

namespace Rivet{
  
  namespace HepMCUtils{
    
    std::vector<ConstGenParticlePtr> particles(ConstGenEventPtr ge){
      return ge->particles();
    }
    
    std::vector<ConstGenParticlePtr> particles(const GenEvent *ge){
      assert(ge);
      return ge->particles();
    }
    
    std::vector<ConstGenVertexPtr> vertices(ConstGenEventPtr ge){
      return ge->vertices();
    }
    
    std::vector<ConstGenVertexPtr> vertices(const GenEvent *ge){
      assert(ge);
      return ge->vertices();
    }
    
    std::vector<ConstGenParticlePtr> particles(ConstGenVertexPtr gv, const Relatives &relo){
      return relo(gv);
    }
    
    std::vector<ConstGenParticlePtr> particles(ConstGenParticlePtr gp, const Relatives &relo){
      return relo(gp);
    }
    
    int uniqueId(ConstGenParticlePtr gp){
      return gp->id();
    }
    
    std::vector<ConstGenParticlePtr> beams(const GenEvent *ge){
      return ge->beams();
    }
    
    bool readEvent(std::shared_ptr<HepMC_IO_type> io, std::shared_ptr<GenEvent> evt){
      io->read_event(*evt);
      return !io->failed();
    }
    
    shared_ptr<HepMC_IO_type> makeReader(std::istream &istr){
      if(&istr == &std::cin) return make_shared<RivetHepMC::ReaderAsciiHepMC2>(istr);
      
      istr.seekg(istr.beg);
      std::string line1, line2;
      
      while(line1.empty()){
        std::getline(istr, line1);
      }
      
      while(line2.empty()){
        std::getline(istr, line2);
      }
      
      istr.seekg(istr.beg);
      
      // if this is absent it doesn't appear to be a HepMC file :(
      if(line1.find("HepMC::Version") == std::string::npos) return nullptr;
      
      shared_ptr<HepMC_IO_type> result;
      
      // Looks like the new HepMC 3 format!
      if(line2.find("HepMC::Asciiv3") != std::string::npos){
        result = make_shared<RivetHepMC::ReaderAscii>(istr);
      }else{
        // assume old HepMC 2 format from here
        result = make_shared<RivetHepMC::ReaderAsciiHepMC2>(istr);
      }
      
      if(result->failed()) result.reset((RivetHepMC::Reader*)nullptr);
      return result;
    }
  }
}
