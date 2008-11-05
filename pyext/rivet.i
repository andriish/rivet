%module rivet

%{
  #define SWIG_FILE_WITH_INIT
  #include "Rivet/Analysis.hh"
  #include "Rivet/AnalysisHandler.hh"
  #include "Rivet/AnalysisLoader.hh"
  #include "Rivet/Tools/Logging.hh"
%}


%include "std_string.i"
%include "std_vector.i"
%include "std_set.i"
%include "std_map.i"
%template(StrList) std::vector<std::string>;
%template(StrSet) std::set<std::string>;
%template(LogLevelMap) std::map<std::string, int>;


%include "Rivet/ParticleName.hh"
%include "Rivet/HistoFormat.hh"


%ignore operator<<;
namespace Rivet {
  %rename(setLogLevel) Log::setLevel(const std::string&, int);
}
%include "Rivet/Tools/Logging.hh"


namespace Rivet {

  class Analysis {
  public:
    virtual std::string getName() const;
    virtual std::string getSpiresId() const;
    virtual std::string getDescription() const;
    virtual std::string getExpt() const;
    virtual std::string getYear() const;
    virtual std::vector<std::string> getReferences() const;
    //virtual const Cuts getCuts() const;
    virtual const BeamPair& getBeams() const;
    //virtual const bool isCompatible(const std::string& quantity, const double value) const;
    virtual const bool isCompatible(const ParticleName& beam1, const ParticleName& beam2) const;
    virtual const bool isCompatible(const BeamPair& beams) const;
    //AnalysisHandler& getHandler() const;
    bool needsCrossSection() const;
  private:
    Analysis();
  };


  class AnalysisHandler {
  public:
    //AnalysisHandler(AIDA::IAnalysisFactory& afac, string basefilename="Rivet", HistoFormat storetype=AIDAML);
    AnalysisHandler(std::string basefilename="Rivet", HistoFormat storetype=AIDAML);
    size_t numEvents() const;
    double sumOfWeights() const;
    AnalysisHandler& addAnalysis(const std::string& analysisname);
    AnalysisHandler& addAnalyses(const std::vector<std::string>& analysisnames);
    void init(int i=0, int N=0);
    void analyze(const HepMC::GenEvent& event);
    void finalize();
    bool needCrossSection();
    AnalysisHandler& setCrossSection(double xs);
    void commitData();
  };


  class AnalysisLoader {
  public:
    static std::set<std::string> getAllAnalysisNames();
    static Analysis* getAnalysis(const std::string& analysisname);
    static void closeAnalysisBuilders();    
  };


}
