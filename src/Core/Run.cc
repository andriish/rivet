// -*- C++ -*-
#include "Rivet/Run.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "zstr/zstr.hpp"
#include <limits>
#include <iostream>

#ifdef ENABLE_HEPMC_3
#include "HepMC3/ReaderFactory.h"
#endif

using std::cout;
using std::endl;

namespace Rivet {
union magic_t {
    uint8_t bytes[4];
    uint32_t number;
};

  Run::Run(AnalysisHandler& ah)
    : _ah(ah), _fileweight(1.0), _xs(NAN)
  { }


  Run::~Run() { }


  Run& Run::setCrossSection(const double xs) {
    _xs = xs;
    return *this;
  }


  Run& Run::setListAnalyses(const bool dolist) {
    _listAnalyses = dolist;
    return *this;
  }


  // Fill event and check for a bad read state
  bool Run::readEvent() {
    /// @todo Clear rather than new the GenEvent object per-event?
    _evt.reset(new GenEvent());
    if (!HepMCUtils::readEvent(_hepmcReader, _evt)){
      Log::getLog("Rivet.Run") << Log::DEBUG << "Read failed. End of file?" << endl;
      return false;
    }
    // Rescale event weights by file-level weight, if scaling is non-trivial
    if (_fileweight != 1.0) {
      for (size_t i = 0; i < (size_t) _evt->weights().size(); ++i) {
        _evt->weights()[i] *= _fileweight;
      }
    }
    return true;
  }



  bool Run::openFile(const std::string& evtfile, double weight) {
    // Set current weight-scaling member
    _fileweight = weight;

    // In case makeReader fails.
    std::string errormessage;

#ifdef ENABLE_HEPMC_3
if (evtfile == "-") {
       printf("Reading events from standard input assuming IO_GenEvent format\n");
#ifdef HAVE_LIBZ
      _istr = make_shared<zstr::istream>(std::cin);
      _hepmcReader = std::shared_ptr<RivetHepMC::Reader>((RivetHepMC::Reader*) ( new RivetHepMC::ReaderAsciiHepMC2(*_istr)));
#else
      _hepmcReader = std::shared_ptr<RivetHepMC::Reader>((RivetHepMC::Reader*) ( new RivetHepMC::ReaderAsciiHepMC2(std::cin)));
#endif
} 
else 
{
_hepmcReader = RivetHepMC::deduce_reader(evtfile);
#ifdef HAVE_LIBZ
if (!_hepmcReader)
{
        printf("No succes with deduction of file type. Test if the file is compressed.\n");
        std::ifstream file_test(evtfile);
        magic_t my_magic = {0x1f, 0x8b, 0x08, 0x08};
        magic_t file_magic;
        file_test.read((char *) file_magic.bytes, sizeof(file_magic));
        if ( file_magic.number == my_magic.number )
        {
        printf("File is compressed\n");
        printf("Reading events from compressed file assuming IO_GenEvent format\n");
        _istr = make_shared<zstr::ifstream>(evtfile);
        _hepmcReader = std::shared_ptr<RivetHepMC::Reader>((RivetHepMC::Reader*) ( new RivetHepMC::ReaderAsciiHepMC2(*_istr)));
        }
}
#endif
}

#else    
    // Set up HepMC input reader objects
    if (evtfile == "-") {
      #ifdef HAVE_LIBZ
      _istr = make_shared<zstr::istream>(std::cin);
      _hepmcReader = HepMCUtils::makeReader(*_istr, &errormessage);
      #else
      _hepmcReader = HepMCUtils::makeReader(std::cin, &errormessage);
      #endif
    } else {
      if ( !fileexists(evtfile) )
        throw Error("Event file '" + evtfile + "' not found");
      #ifdef HAVE_LIBZ
      // NB. zstr auto-detects if file is deflated or plain-text
      _istr = make_shared<zstr::ifstream>(evtfile.c_str());
      #else
      _istr = make_shared<std::ifstream>(evtfile.c_str());
      #endif
      _hepmcReader = HepMCUtils::makeReader(*_istr, &errormessage);

    }
#endif

    if (_hepmcReader == nullptr) {
      Log::getLog("Rivet.Run")
        << Log::ERROR << "Read error on file '" << evtfile << "' "
        << errormessage << endl;
      return false;
    }
    return true;
  }


  bool Run::init(const std::string& evtfile, double weight) {
    if (!openFile(evtfile, weight)) return false;

    // Read first event to define run conditions
    bool ok = readEvent();
    if (!ok) return false;
    if(HepMCUtils::particles(_evt).size() == 0){
      Log::getLog("Rivet.Run") << Log::ERROR << "Empty first event." << endl;
      return false;
    }

    // Initialise AnalysisHandler with beam information from first event
    _ah.init(*_evt);

    // Set cross-section from command line
    if (!std::isnan(_xs)) {
      Log::getLog("Rivet.Run")
        << Log::DEBUG << "Setting user cross-section = " << _xs << " pb" << endl;

      _ah.setCrossSection(make_pair(_xs, 0.0));
    }

    // List the chosen & compatible analyses if requested
    if (_listAnalyses) {
      for (const std::string& ana : _ah.analysisNames()) {
        cout << ana << endl;
      }
    }

    return true;
  }


  bool Run::processEvent() {
    // Analyze event
    _ah.analyze(*_evt);

    return true;
  }


  bool Run::finalize() {
    _evt.reset();

    _ah.finalize();

    return true;
  }


}
