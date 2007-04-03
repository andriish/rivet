// -*- C++ -*-
#ifndef RIVET_HZ95108_H
#define RIVET_HZ95108_H

#include "Rivet/Analysis/Analysis.hh"
#include "Rivet/Projections/FinalStateHCM.hh"
#include "Rivet/Projections/CentralEtHCM.hh"
#include "Rivet/RivetAIDA.fhh"


namespace Rivet {

  /// This analysis simply measures the total multiplicity.
  class HZ95108 : public Analysis {

  public:

    /// The default constructor.
    inline HZ95108()
      : lepton(beams, 11, -11), diskin(beams, lepton, 2212),
        fsproj(lepton, diskin, fsp), y1hcm(fsproj) 
    { }

  public:

    /// The name of this analysis is "HZ95108"
    inline std::string name() const {
      return "HZ95108";
    }


    /**
     * Initialize this analysis object. A concrete class should here
     * book all necessary histograms. An overridden function must make
     * sure it first calls the base class function.
     */
    virtual void init();

    /**
     * Analyze one event. A concrete class should here apply the
     * necessary projections on the \a event and fill the relevant
     * histograms. An overridden function must make sure it first calls
     * the base class function.
     */
    virtual void analyze(const Event & event);

    /**
     * Finalize this analysis object. A concrete class should here make
     * all necessary operations on the histograms. Writing the
     * histograms to a file is, however, done by the Rivet class. An
     * overridden function must make sure it first calls the base class
     * function.
     */
    virtual void finalize();

    /**
     * Return the RivetInfo object of this analysis object. Derived
     * classes should re-implement this function to return the combined
     * RivetInfo object of this object and of any Projection objects
     * upon which this depends.
     */
    virtual RivetInfo getInfo() const;

  protected:

    /**
     * Calculate the bin number from the DISKinematics projection.
     */
    int getbin(const DISKinematics & dk);


  private:

    /**
     * The Beam projector used.
     */
    Beam beams;

     /**
     * The DISLepton projector used.
     */
    DISLepton lepton;

     /**
     * The DISKinematics projector used.
     */
    DISKinematics diskin;

   /**
     * The FinalState projector used.
     */
    FinalState fsp;

   /**
     * The FinalStateHCM projector used.
     */
    FinalStateHCM fsproj;

    /**
     * The CentralEtHCM projector used.
     */
    CentralEtHCM y1hcm;

    /**
     * Some integer constants used.
     */
    static const int nb = 24, nbin = 9;

    /**
     * Some double constants used.
     */
    static const double xmin, xmax;

    /**
     * Histograms for the Et flows
     */
    vector<AIDA::IHistogram1D *> hEtFlow, hEtFlowStat;

    /**
     * Histograms for averages in different kinematical bins.
     */
    AIDA::IHistogram1D * hAvEt, * hAvX, * hAvQ2, * hN;

    /**
     * Helper vector;
     */
    vector<double> nev;

  private:

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    HZ95108 & operator=(const HZ95108 &);

  };

}

#endif /* RIVET_HZ95108_H */
