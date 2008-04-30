// -*- C++ -*-
#ifndef RIVET_DELPHI_1996_S3430090_HH
#define RIVET_DELPHI_1996_S3430090_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Implementation of DELPHI event shape paper.
  /// @author Andy Buckley
  class DELPHI_1996_S3430090 : public Analysis {

  public:

    /// Default constructor.
    DELPHI_1996_S3430090() 
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(*new Beam(), "Beams");
      const ChargedFinalState& cfs = addProjection(*new ChargedFinalState(), "FS");
      #ifdef HAVE_JADE
      addProjection(*new FastJets(cfs, FastJets::JADE, 0.7), "JadeJets");
      addProjection(*new FastJets(cfs, FastJets::DURHAM, 0.7), "DurhamJets");
      #endif
      addProjection(*new Sphericity(cfs), "Sphericity");
      addProjection(*new ParisiTensor(cfs), "Parisi");
      const Thrust& thrust = addProjection(*new Thrust(cfs), "Thrust");
      addProjection(*new Hemispheres(thrust), "Hemispheres");
    }


    /// Factory method.
    static Analysis* create() { 
      return new DELPHI_1996_S3430090(); 
    }


    /// @name Publication metadata
    //@{
    /// Get a description of the analysis.
    string getSpiresId() const {
      return "3430090";
    }
    /// Get a description of the analysis.
    string getDescription() const {
      return "Delphi MC tuning on event shapes and identified particles.";
    }
    /// Experiment which performed and published this analysis.
    string getExpt() const {
      return "DELPHI";
    }
    /// When published (preprint year according to SPIRES).
    string getYear() const {
      return "1996";
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /// Hide the assignment operator
    DELPHI_1996_S3430090& operator=(const DELPHI_1996_S3430090&);

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the 
    /// inclusive single particle distributions' normalisations.
    double _weightedTotalPartNum;

    /// @name Histograms
    //@{
    AIDA::IHistogram1D *_histPtTIn;
    AIDA::IHistogram1D *_histPtTOut;
    AIDA::IHistogram1D *_histPtSIn;
    AIDA::IHistogram1D *_histPtSOut;

    AIDA::IHistogram1D *_histRapidityT;
    AIDA::IHistogram1D *_histRapidityS;

    AIDA::IHistogram1D *_histScaledMom, *_histLogScaledMom;

    AIDA::IProfile1D   *_histPtTOutVsXp, *_histPtVsXp;

    AIDA::IHistogram1D *_hist1MinusT; 
    AIDA::IHistogram1D *_histTMajor; 
    AIDA::IHistogram1D *_histTMinor; 
    AIDA::IHistogram1D *_histOblateness; 

    AIDA::IHistogram1D *_histSphericity;
    AIDA::IHistogram1D *_histAplanarity;
    AIDA::IHistogram1D *_histPlanarity;

    AIDA::IHistogram1D *_histCParam;
    AIDA::IHistogram1D *_histDParam;

    AIDA::IHistogram1D *_histHemiMassD;
    AIDA::IHistogram1D *_histHemiMassH;
    AIDA::IHistogram1D *_histHemiMassL;
               
    AIDA::IHistogram1D *_histHemiBroadW;
    AIDA::IHistogram1D *_histHemiBroadN;
    AIDA::IHistogram1D *_histHemiBroadT;
    AIDA::IHistogram1D *_histHemiBroadD;

    AIDA::IHistogram1D *_histDiffRate2Durham;
    AIDA::IHistogram1D *_histDiffRate2Jade; 
    AIDA::IHistogram1D *_histDiffRate3Durham;
    AIDA::IHistogram1D *_histDiffRate3Jade;
    AIDA::IHistogram1D *_histDiffRate4Durham;
    AIDA::IHistogram1D *_histDiffRate4Jade;

    AIDA::IHistogram1D *_histEEC, *_histAEEC;

    AIDA::IHistogram1D *_histMultiCharged;
    //@}
  
  };

}

#endif
