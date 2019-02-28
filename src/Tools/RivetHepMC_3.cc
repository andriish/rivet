// -*- C++ -*-

#include "Rivet/Tools/RivetHepMC.hh"
#include "Rivet/Tools/ReaderCompressedAscii.hh"
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
    
    int particles_size(ConstGenEventPtr ge){
      return particles(ge).size();
    }

    int particles_size(const GenEvent *ge){
      return particles(ge).size();
    }
    
    int uniqueId(ConstGenParticlePtr gp){
      return gp->id();
    }
    
    std::pair<ConstGenParticlePtr,ConstGenParticlePtr> beams(const GenEvent *ge) {
      std::vector<ConstGenParticlePtr> beamlist = ge->beams();
      if ( beamlist.size() < 2 ) return std::pair<ConstGenParticlePtr,ConstGenParticlePtr>();
      return std::make_pair(beamlist[0], beamlist[1]);
    }
    
    bool readEvent(std::shared_ptr<HepMC_IO_type> io, std::shared_ptr<GenEvent> evt){
      return io->read_event(*evt) && !io->failed();
    }
    
    shared_ptr<HepMC_IO_type> makeReader(std::istream & istr,
                                         std::string * errm) {
      shared_ptr<HepMC_IO_type> ret;
      
      // First scan forward and check if there is some hint as to what
      // kind of file we are looking att.
      int ntry = 10;
      std::string header;
      int filetype = -1;
      while  ( ntry ) {
        std::getline(istr, header);
        if ( header.empty() ) continue;
        if ( header.substr(0, 34) == "HepMC::Asciiv3-START_EVENT_LISTING" ) {
          filetype = 3;
          break;
        }
        if ( header.substr(0, 44) == "HepMC::CompressedAsciiv3-START_EVENT_LISTING" ) {
          filetype = 4;
          break;
        }
        if ( header.substr(0, 38) == "HepMC::IO_GenEvent-START_EVENT_LISTING" ) {
          filetype = 2;
          break;
        }
        --ntry;
      }

      
      if ( filetype == 3 )
        ret = make_shared<RivetHepMC::ReaderAscii>(istr);
      else if ( filetype == 4 )
        ret = make_shared<Rivet::ReaderCompressedAscii>(istr);
      else
        ret = make_shared<RivetHepMC::ReaderAsciiHepMC2>(istr);
      if ( filetype == 0 && errm )
        *errm += "Could not determine file type. Assuming HepMC2 file. ";
      // Check that everything was ok.
      if ( ret->failed() ) {
        if ( errm ) *errm = "Problems reading from HepMC file.";
        ret = shared_ptr<HepMC_IO_type>();
      }

      return ret;
    }

  }
}
