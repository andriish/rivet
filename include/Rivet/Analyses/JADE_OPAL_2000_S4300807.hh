// -*- C++ -*-
#ifndef RIVET_JADE_OPAL_2000_S4300807_HH
#define RIVET_JADE_OPAL_2000_S4300807_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {

  /**
   * @brief Jet rates in e+e- at OPAL and JADE
   * @author Frank Siegert
   *
   * @par Run conditions
   *
   * @arg LEP1 beam energy: \f$ \sqrt{s} = \$f 91.2 GeV
   * @arg Run with generic QCD events.
   */
  class JADE_OPAL_2000_S4300807 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    JADE_OPAL_2000_S4300807(double sqrts, int nr_R_Jade,
                            int nr_R_Durham, int nr_y_Durham) :
      _sqrts(sqrts), _nr_R_Jade(nr_R_Jade),
      _nr_R_Durham(nr_R_Durham), _nr_y_Durham(nr_y_Durham)
    {
      setBeams(ELECTRON, POSITRON); 
      addProjection(Beam(), "Beams");
      const FinalState fs;
      addProjection(fs, "FS");
      #ifdef HAVE_JADE
      addProjection(FastJets(fs, FastJets::JADE, 0.7), "JadeJets");
      addProjection(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
      #endif
    }

    //@}


    /// @name Publication metadata
    //@{
    /// A short description of the analysis.
    string spiresId() const {
      return "4300807";
    }
    /// A short description of the analysis.
    virtual string summary() const {
      return "Jet rates in e+e- at JADE [35-43 GeV] and OPAL [91.2-189 GeV].";
    }
    /// Experiment which performed and published this analysis.
    string experiment() const {
      return "JADE_OPAL";
    }
    /// Collider on which the experiment ran.
    string collider() const {
      if (_sqrts<90.0) {
        return "DESY PETRA";
      }
      else if (_sqrts>90.0 && _sqrts<92.0) {
        return "LEP Run I";
      }
      else {
        return "LEP Run 2";
      }
    }
    /// When published (preprint year according to SPIRES).
    string year() const {
      return "2000";
    }
    /// Names & emails of analysis authors.
    vector<string> authors() const {
      vector<string> ret;
      ret += "Frank Siegert <frank.siegert@durham.ac.uk>";
      return ret;
    }
    /// A full description of the analysis.
    string description() const {
      ostringstream os;
      os << "Differential and integrated jet rates for Durham and JADE " 
         << "jet algorithms at sqrt(s) = " << _sqrts << ".";
      return os.str();
    }
    /// Information about the events needed as input for this analysis.
    string runInfo() const {
      ostringstream os;
      os << "e+ e- collisions: " << endl << endl
         << "* e+ e- -> jet jet (+ jets) at " << _sqrts << " GeV. "
         << "* no cuts needed" << endl;
      return os.str();
    }
    string status() const {
      return "VALIDATED";
    }
    /// Journal, and preprint references.
    vector<string> references() const {
      vector<string> ret;
      ret.push_back("Eur.Phys.J.C17:19-51,2000");
      ret.push_back("arXiv:hep-ex/0001055");
      return ret;
    }
    //@}


    /// @name Analysis methods
    //@{
    virtual void init();
    virtual void analyze(const Event& event);
    virtual void finalize();
    //@}


  private:

    /// @name Histograms
    //@{
    AIDA::IDataPointSet *_h_R_Jade[5];
    AIDA::IDataPointSet *_h_R_Durham[5];

    AIDA::IHistogram1D *_h_y_Durham[4];
    //@}

    double _sqrts;
    int _nr_R_Jade, _nr_R_Durham, _nr_y_Durham;

  };


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class JADE_OPAL_2000_S4300807_35GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_35GEV() :
      JADE_OPAL_2000_S4300807(35.0, 7, 16, 24) {}

    static Analysis* create() { return new JADE_OPAL_2000_S4300807_35GEV(); }

    string summary() const { return "Jet rates in e+e- at JADE [35 GeV]."; }
  };


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class JADE_OPAL_2000_S4300807_44GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_44GEV() :
      JADE_OPAL_2000_S4300807(44.0, 8, 17, 25) {}

    static Analysis* create() { return new JADE_OPAL_2000_S4300807_44GEV(); }

    string summary() const { return "Jet rates in e+e- at JADE [44 GeV]."; }
  };


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class JADE_OPAL_2000_S4300807_91GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_91GEV() :
      JADE_OPAL_2000_S4300807(91.2, 9, 18, 26) {}

    static Analysis* create() { return new JADE_OPAL_2000_S4300807_91GEV(); }

    string summary() const { return "Jet rates in e+e- at OPAL [91 GeV]."; }
  };


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class JADE_OPAL_2000_S4300807_133GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_133GEV() :
      JADE_OPAL_2000_S4300807(133.0, 10, 19, 27) {}

    static Analysis* create() { return new JADE_OPAL_2000_S4300807_133GEV(); }

    string summary() const { return "Jet rates in e+e- at OPAL [133 GeV]."; }
  };


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class JADE_OPAL_2000_S4300807_161GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_161GEV() :
      JADE_OPAL_2000_S4300807(161.0, 11, 20, 28) {}

    static Analysis* create() { return new JADE_OPAL_2000_S4300807_161GEV(); }

    string summary() const { return "Jet rates in e+e- at OPAL [161 GeV]."; }
  };


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class JADE_OPAL_2000_S4300807_172GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_172GEV() :
      JADE_OPAL_2000_S4300807(172.0, 12, 21, 29) {}

    static Analysis* create() { return new JADE_OPAL_2000_S4300807_172GEV(); }

    string summary() const { return "Jet rates in e+e- at OPAL [172 GeV]."; }
  };


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class JADE_OPAL_2000_S4300807_183GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_183GEV() :
      JADE_OPAL_2000_S4300807(183.0, 13, 22, 30) {}

    static Analysis* create() { return new JADE_OPAL_2000_S4300807_183GEV(); }

    string summary() const { return "Jet rates in e+e- at OPAL [183 GeV]."; }
  };


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  class JADE_OPAL_2000_S4300807_189GEV : public JADE_OPAL_2000_S4300807 {
  public:
    JADE_OPAL_2000_S4300807_189GEV() :
      JADE_OPAL_2000_S4300807(189.0, 14, 23, 31) {}

    static Analysis* create() { return new JADE_OPAL_2000_S4300807_189GEV(); }

    string summary() const { return "Jet rates in e+e- at OPAL [189 GeV]."; }
  };

}

#endif
