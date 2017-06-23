// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include <random>
#if defined(_OPENMP)
#include "omp.h"
#endif

namespace Rivet {


  // Return a uniformly sampled random number between 0 and 1
  double rand01() {
    // return rand() / (double)RAND_MAX;
    //
    // static random_device rd;
    // static mt19937 gen(rd());
    // return generate_canonical<double, 10>(gen);
    //
    // static mt19937 gen(12345);
    // return generate_canonical<double, 10>(gen);
    //
    #if defined(_OPENMP)
    static map<int,mt19937> gens;
    const int nthread = omp_get_thread_num();
    if (gens.find(nthread) == gens.end()) {
      seed_seq seq{1,2,3,4,5};
      vector<uint32_t> seeds(nthread+1);
      seq.generate(seeds.begin(), seeds.end());
      gens[nthread] = mt19937(seeds[nthread]);
      // cout << "Thread " << nthread+1 << ", seed=" << seeds[nthread] << " (" << gens.size() << " RNGs)" << endl;
    }
    mt19937& g = gens[nthread];
    #else
    static mt19937 g(12345);
    #endif
    const double r = generate_canonical<double, 10>(g);
    return r;
  }


}
