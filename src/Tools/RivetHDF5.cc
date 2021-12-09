// -*- C++ -*-
#include "Rivet/Tools/RivetHDF5.hh"

namespace Rivet {
  namespace H5 {


    // /// Read HDF5 data from datatset @a dsname in file @a filename, into the provided @a rtndata container
    // template <typename T>
    // bool readData(const string& filename, const string& dsname, T& rtndata) {
    //   const string filepath = findAnalysisDataFile(filename);
    //   if (filepath.empty()) return false;

    //   try {
    //     HighFive::File h5file(filepath, HighFive::File::ReadOnly);

    //     DataSet dataset = h5file.getDataSet(dsname);
    //     dataset.read(rtndata);
    //   } catch (...) {
    //     return false;
    //   }
    //   return true;
    // }


    // /// Read HDF5 data from datatset @a dsname in file @a filename, into a newly constructed container
    // template <typename T>
    // T readData(const string& filename, const string& dsname) {
    //   T rtn;
    //   readData(filename, dsname, rtn);
    //   return rtn;
    // }


    void _h5_compile_test() {
      auto x = readData<vector<int>>("foobar.h5", "/foo/bar/baz");
      auto y = readData<vector<vector<double>>>("foobar.h5", "/foo/bar/baz");
      auto z = readData<vector<int>>("foobar.h5", "/foo/bar/baz");
    }


  }
}
