%module rivet

%{
  #define SWIG_FILE_WITH_INIT
  #include "Rivet/Analysis.hh"
  #include "Rivet/AnalysisHandler.hh"
  #include "Rivet/AnalysisLoader.hh"
  #include "Rivet/Run.hh"
  #include "Rivet/Tools/Logging.hh"
  #include "Rivet/Event.hh"
  #include "Rivet/Particle.hh"
  #include "Rivet/ParticleName.hh"
  #include "Rivet/Projections/Beam.hh"
%}

// STL stuff
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_map.i"
%template(StrList) std::vector<std::string>;
%template(DblPair) std::pair<double, double>;
%template(DblPairList) std::vector< std::pair<double, double> >;

// Histo format enum
%include "Rivet/HistoFormat.hh"

// Particle ID stuff
%include "Rivet/Particle.fhh"
%include "Rivet/ParticleName.hh"
%template(PdgIdPair) std::pair<Rivet::PdgId,Rivet::PdgId>;
%template(PdgIdPairList) std::vector<Rivet::BeamPair>;
//%template(PdgIdPairList) std::vector< std::pair<PdgId,PdgId> >;

// Logging interface
%template(LogLevelMap) std::map<std::string, int>;
%ignore operator<<;
namespace Rivet {
  %rename(setLogLevel) Log::setLevel(const std::string&, int);
}
%include "Rivet/Tools/Logging.hh"


// Rivet class mappings
namespace Rivet {

  std::string version();

  class Event {
    Event();
    Event(const HepMC::GenEvent&);
    const HepMC::GenEvent& genEvent() const;
    double weight() const;
  };

  class Particle {
    Particle();
    bool hasGenParticle() const;
    const HepMC::GenParticle& genParticle() const;
    const long pdgId() const;
  };

  ParticlePair beams(const Event& e);

  BeamPair beamIds(const HepMC::GenEvent& e) {
    return beamIds(Event(e));
  }

  double sqrtS(const Event& e);


  // Mapping of just the metadata parts of the Analysis API
  class Analysis {
  public:
    virtual std::string name() const;
    virtual std::string spiresId() const;
    virtual std::string summary() const;
    virtual std::string description() const;
    virtual std::string runInfo() const;
    virtual std::string experiment() const;
    virtual std::string collider() const;
    virtual std::string year() const;
    virtual const std::vector<BeamPair>& requiredBeams() const;
    virtual const std::vector<std::pair<double,double> >& energies() const;
    virtual std::vector<std::string> authors() const;
    virtual std::vector<std::string> references() const;
    virtual std::vector<std::string> todos() const;
    virtual std::string status() const;
    virtual std::string bibKey() const;
    virtual std::string bibTeX() const;
    virtual const bool isCompatible(const ParticleName& beam1, 
                                    const ParticleName& beam2) const;
    virtual const bool isCompatible(const BeamPair& beams) const;
    //AnalysisHandler& handler() const;
    bool needsCrossSection() const;
  private:
    Analysis();
  };


  class AnalysisHandler {
  public:
    AnalysisHandler(std::string basefilename="Rivet", 
                    std::string runname="", 
                    HistoFormat storetype=AIDAML);
    std::string runName() const;
    size_t numEvents() const;
    double sumOfWeights() const;
    double sqrtS() const;
    const ParticlePair& beams() const;
    const BeamPair& beamIds() const;
    std::vector<std::string> analysisNames();
    AnalysisHandler& addAnalysis(const std::string& analysisname);
    AnalysisHandler& addAnalyses(const std::vector<std::string>& analysisnames);
    AnalysisHandler& removeAnalysis(const std::string& analysisname);
    AnalysisHandler& removeAnalyses(const std::vector<std::string>& analysisnames);
    AnalysisHandler& removeIncompatibleAnalyses(const BeamPair& beams);
    void init();
    void init(const HepMC::GenEvent& event);
    void analyze(const HepMC::GenEvent& event);
    void finalize();
    bool needCrossSection();
    AnalysisHandler& setCrossSection(double xs);
    void commitData();
  };


  class AnalysisLoader {
  public:
    static std::vector<std::string> analysisNames();
    static Analysis* getAnalysis(const std::string& analysisname);
  };


}

%include "Rivet/Run.hh"
