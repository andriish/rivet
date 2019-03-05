// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include <random>
#if defined(_OPENMP)
#include "omp.h"
#endif

namespace Rivet {


  // Return a thread-safe random number generator
  std::mt19937& rng() {
    #if defined(_OPENMP)
    static map<int,std::mt19937> gens;
    const int nthread = omp_get_thread_num();
    if (gens.find(nthread) == gens.end()) {
      // Make seeds for each thread, either via the standard seed generator or based on a fixed seed from the environment
      vector<uint32_t> seeds(nthread+1);
      const uint32_t envseed = getEnvParam<uint32_t>("RIVET_RANDOM_SEED", 0);
      if (envseed > 0) {
        std::iota(seeds.begin(), seeds.end(), envseed);
      } else {
        std::seed_seq seq{1,2,3,4,5};
        seq.generate(seeds.begin(), seeds.end());
      }
      gens[nthread] = std::mt19937(seeds[nthread]);
      // cout << "Thread " << nthread+1 << ", seed=" << seeds[nthread] << " (" << gens.size() << " RNGs)" << '\n';
    }
    std::mt19937& g = gens[nthread];
    #else
    static std::mt19937 g(12345);
    #endif
    return g;
  }


  // Return a uniformly sampled random number between 0 and 1
  double rand01() {
    // return rand() / (double)RAND_MAX;
    return std::generate_canonical<double, 32>(rng()); ///< @todo What's the "correct" number of bits of randomness?
  }


  // Return a Gaussian/normal sampled random number with the given mean and width
  double randnorm(double loc, double scale) {
    std::normal_distribution<> d(loc, scale);
    return d(rng());
  }


  // Return a log-normal sampled random number
  double randlognorm(double loc, double scale) {
    std::lognormal_distribution<> d(loc, scale);
    return d(rng());
  }


}
