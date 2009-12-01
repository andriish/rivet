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


// Histo format enum
%include "Rivet/HistoFormat.hh"


// Particle ID stuff
%include "Rivet/ParticleName.hh"
%template(BeamPair) std::pair<long,long>;


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

  typedef std::vector<Particle> ParticleVector;
  typedef std::pair<Particle, Particle> ParticlePair;

  ParticlePair beams(const Event& e);

  BeamPair beamIds(const HepMC::GenEvent& e) {
    return beamIds(Event(e));
  }

  double sqrtS(const Event& e);


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
    virtual std::string status() const;
    virtual const std::vector<std::pair<double,double> >& energies() const;
    virtual std::vector<std::string> authors() const;
    virtual std::vector<std::string> references() const;
    virtual const BeamPair& beams() const;
    virtual const BeamPair& requiredBeams() const;
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
    AnalysisHandler(std::string basefilename="Rivet", std::string runname="", 
                    HistoFormat storetype=AIDAML);
    std::string runName() const;
    size_t numEvents() const;
    double sumOfWeights() const;
    std::vector<std::string> analysisNames();
    AnalysisHandler& addAnalysis(const std::string& analysisname);
    AnalysisHandler& addAnalyses(const std::vector<std::string>& analysisnames);
    AnalysisHandler& removeAnalysis(const std::string& analysisname);
    AnalysisHandler& removeAnalyses(const std::vector<std::string>& analysisnames);
    AnalysisHandler& removeIncompatibleAnalyses(const BeamPair& beams);
    void init(int i=0, int N=0);
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
