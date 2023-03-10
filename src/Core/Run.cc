// -*- C++ -*-
#include "Rivet/Run.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "Rivet/Tools/RivetHepMC.hh"
#include "zstr/zstr.hpp"
#include <limits>
#include <iostream>

using std::cout;
using std::endl;

namespace Rivet {


  /// Byte/number conversion via union, for HepMC file inspection
  union magic_t {
    uint8_t bytes[2];
    uint16_t number;
  };


  Run::Run(AnalysisHandler& ah)
    : _ah(ah)
  { }


  Run::~Run() { }


  Run& Run::setCrossSection(double xs) {
    _xs = xs;
    return *this;
  }


  Run& Run::setListAnalyses(bool dolist) {
    _listAnalyses = dolist;
    return *this;
  }


  Log& Run::getLog() const {
    return Log::getLog("Rivet.Run");
  }


  // Fill event and check for a bad read state
  bool Run::readEvent() {
    /// @todo Clear rather than new the GenEvent object per-event?
    _evt.reset(new GenEvent());
    if (!HepMCUtils::readEvent(_hepmcReader, _evt)) {
      MSG_DEBUG("Read failed. End of file?");
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

    #ifdef RIVET_ENABLE_HEPMC_3
    if (evtfile == "-") {
      // Turn off the buffering to make IO faster and make ungetc work on cin
      std::basic_ios<char>::sync_with_stdio(false);
      #ifdef HAVE_LIBZ
      _istr = make_shared<zstr::istream>(std::cin);
      #else
      _istr = make_shared<std::istream>(std::cin);
      #endif
      // Use standard HepMC3 deduction on stream. For HepMC3 < 3.2.0 the function is implemented in Rivet
      _hepmcReader = RivetHepMC::deduce_reader(*_istr);
    } else {
      // Use standard HepMC3 deduction on file
      _hepmcReader = RivetHepMC::deduce_reader(evtfile);
      // Check if the file is compressed, if the deduction fails
      /// @todo Can we move this into the RivetHepMC.hh header? This is a *lot* of HepMC-specific noise for the Run manager class
      if (!_hepmcReader) {
        MSG_INFO("No success with deduction of file type. Test if the file is compressed");
        std::ifstream file_test(evtfile);
        magic_t my_magic = {{0x1f, 0x8b}};
        magic_t file_magic;
        file_test.read((char *) file_magic.bytes, sizeof(file_magic));
        if (file_magic.number == my_magic.number) {
          MSG_INFO("File is compressed");
          #ifdef HAVE_LIBZ
          _istr = make_shared<zstr::ifstream>(evtfile);
          _hepmcReader = RivetHepMC::deduce_reader(*_istr);
          #else
          MSG_INFO("No zlib support.");
          #endif
        } else {
          // File is not compressed. Open stream and let the code below to handle it
          MSG_INFO("File is not compressed. No succes with deduction of file type.");
          _istr = make_shared<std::ifstream>(evtfile);
        }
      }
    }

    // If the deduction has failed, check custom formats
    /// @todo Move this into the RivetHepMC.hh header: this is a *lot* of HepMC-specific noise for the Run manager class
    if (!_hepmcReader) {
      std::vector<std::string> head;
      head.push_back("");
      size_t back=0;
      size_t backnonempty=0;
      while (back < 200 && backnonempty < 100 && _istr) { ///< @todo Document this
        char c = _istr->get();
        back++;
        if (c == '\n') {
          if (head.back().length() != 0)
            head.push_back("");
        } else {
          head.back() += c;
          backnonempty++;
        }
      }
      if (!_istr) MSG_INFO("Info in deduce_reader: input stream is too short or invalid.");
      for (size_t i = 0; i < back; ++i) _istr->unget();
      if (strncmp(head.at(0).c_str(), "HepMC::Version", 14) == 0 &&
          strncmp(head.at(1).c_str(), "HepMC::CompressedAsciiv3-START_EVENT_LISTING", 44) == 0) {
        MSG_INFO("Info in deduce_reader: Attempt CompressedAsciiv3");
        //_hepmcReader= make_shared<Rivet::RivetHepMC::ReaderCompressedAscii>(_istr);
      }
    }
    #endif


    #ifndef RIVET_ENABLE_HEPMC_3
    // Use Rivet's own file format deduction (which uses the one in
    // HepMC3 if needed).
    _hepmcReader = HepMCUtils::makeReader(evtfile, _istr, &errormessage);

    // Check that it worked.
    #endif
    if (_hepmcReader == nullptr) {
      MSG_ERROR("Read error in file '" << evtfile << "' " << errormessage);
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
      MSG_ERROR("Empty first event.");
      return false;
    }

    // Initialise AnalysisHandler with beam information from first event
    _ah.init(*_evt);

    // Set cross-section from command line
    if (notNaN(_xs)) {
      MSG_DEBUG("Setting user cross-section = " << _xs << " pb");
      _ah.setCrossSection(make_pair(_xs, 0.0), true);
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
    _ah.analyze(*_evt);
    return true;
  }


  bool Run::finalize() {
    _evt.reset();
    _ah.finalize();
    return true;
  }


}
