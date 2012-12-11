// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class ATLAS_2011_I894867 : public Analysis {
  public:

    ATLAS_2011_I894867()
      : Analysis("ATLAS_2011_I894867")
    {    }

  public:

    void init() {
      addProjection(FinalState(), "FS");
      _h_sigma = bookHistogram1D(1, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      const FinalState& fs = applyProjection<FinalState>(event, "FS");
      if (fs.size() < 2) vetoEvent; // need at least two particles to calculate gaps

      // Calculate gap sizes and midpoints
      const ParticleVector particlesByEta = fs.particlesByEta(); // sorted from minus to plus
      vector<double> gaps, midpoints;
      for (size_t ip = 1; ip < particlesByEta.size(); ++ip) {
        const Particle& p1 = particlesByEta[ip-1];
        const Particle& p2 = particlesByEta[ip];
        const double gap = p2.momentum().eta() - p1.momentum().eta();
        const double mid = (p2.momentum().eta() + p1.momentum().eta()) / 2.;
        assert(gap >= 0);
        gaps.push_back(gap);
        midpoints.push_back(mid);
      }

      // Find the midpoint of the largest gap to separate X and Y systems
      int imid = std::distance(gaps.begin(), max_element(gaps.begin(), gaps.end()));
      double gapcenter = midpoints[imid];

      FourMomentum mxFourVector, myFourVector;
      foreach(const Particle& p, fs.particlesByEta()) {
        if (p.momentum().eta() > gapcenter) {
          mxFourVector += p.momentum();
        } else {
          myFourVector += p.momentum();
        }
      }
      const double M2 = max(mxFourVector.mass2(), myFourVector.mass2());
      const double xi = M2/sqr(sqrtS()); // sqrt(s)=7000 GeV, note that units cancel
      if (xi < 5*10e-6) vetoEvent;

      _h_sigma->fill(sqrtS()/GeV, weight);
    }


    void finalize() {
      scale(_h_sigma, crossSection()/millibarn/sumOfWeights());
    }


  private:

    AIDA::IHistogram1D* _h_sigma;

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2011_I894867);

}
