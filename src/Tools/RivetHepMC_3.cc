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
      if ( beamlist.size() < 2 ) {
        std::cerr << "CANNOT FIND ANY BEAMS!" << std::endl;
        return std::pair<ConstGenParticlePtr,ConstGenParticlePtr>();
      }
      return std::make_pair(beamlist[0], beamlist[1]);
    }
    
    bool readEvent(std::shared_ptr<HepMC_IO_type> io, std::shared_ptr<GenEvent> evt){
      io->read_event(*evt);
      return !io->failed();
    }
    
    shared_ptr<HepMC_IO_type> makeReader(std::istream & istr,
                                         std::string * errm) {
      shared_ptr<HepMC_IO_type> ret;
      
      // First scan forward and check if there is some hint as to what
      // kind of file we are looking att.
      int nchar = 256;
      std::string header;
      while ( nchar-- && !istr.eof() )
        header += char(istr.get());
      // If this stream was too short to contain any reasonable number
      // of events, just give up.
      if ( !istr ) {
        if ( errm ) *errm = "Could not find HepMC header information";
        return shared_ptr<HepMC_IO_type>();
      }

      // Now reset the stream to its original state ...
      for ( int i = header.length() - 1; i >= 0; --i )
        istr.putback(header[i]);

      // ... and check which kind of format it was and create the
      // corresponding reader. First try the HepM3 ascii format.
      if ( header.find("HepMC::Asciiv3-START_EVENT_LISTING") !=
           std::string::npos )
        ret = make_shared<RivetHepMC::ReaderAscii>(istr);

      // Check if the file is written by WriterRoot or WriterRootTree.
      else if ( header.substr(0, 4) == "root" ) {
        if ( errm ) *errm = "Rivet cancurrently not read HepMC root files.";
        return ret;
      }

      // The default is to assume it is a good old HepMC2 ascii file.
      else {
        ret = make_shared<RivetHepMC::ReaderAsciiHepMC2>(istr);
      }

      // Check that everything was ok.
      if ( ret->failed() ) {
        if ( errm ) *errm = "Problems reading from HepMC file.";
        ret = shared_ptr<HepMC_IO_type>();
      }

      return ret;
    }

  }
}
