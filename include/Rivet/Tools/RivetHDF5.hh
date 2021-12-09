// -*- C++ -*-
#ifndef RIVET_RivetHDF5_HH
#define RIVET_RivetHDF5_HH
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetPaths.hh"
#include "highfive/H5File.hpp"

namespace Rivet {
  namespace H5 {


    using namespace HighFive;


    /// Read HDF5 file @a filename
    HighFive::File readFile(const string& filename) {
      const string filepath = findAnalysisDataFile(filename);
      if (filepath.empty()) throw IOError("Failed to load HDF5 file " + filename);
      try {
        return HighFive::File(filepath, HighFive::File::ReadOnly);
      } catch (...) {
        throw IOError("Failed to load HDF5 file " + filename);
      }
    }


    /// Read HDF5 data from dataset @a dsname in file @a filename, into the provided @a rtndata container
    template <typename T>
    bool readData(const string& filename, const string& dsname, T& rtndata) {
      try {
        HighFive::File h5file = readFile(filename);
        DataSet dataset = h5file.getDataSet(dsname);
        dataset.read(rtndata);
      } catch (...) {
        return false;
      }
      return true;
    }


    /// Read HDF5 data from dataset @a dsname in file @a filename, into a newly constructed container
    template <typename T>
    T readData(const string& filename, const string& dsname) {
      T rtn;
      readData(filename, dsname, rtn);
      return rtn;
    }


  }
}

#endif
