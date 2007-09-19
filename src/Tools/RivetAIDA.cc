#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Tools/Utils.hh"
#include "LWH/AnalysisFactory.h"
#include "TinyXML/tinyxml.h"
#include <sstream>
using namespace std;



/// "Plugin" function to return an AIDA system (LWH impl.)
extern "C" AIDA::IAnalysisFactory* AIDA_createAnalysisFactory() {
  return new LWH::AnalysisFactory();
}



namespace Rivet {

  // Forward declaration of generated function.
  const string getInstalledDataPath();


  const string getDataPath(string papername) {
    return getInstalledDataPath() + "/" + papername + ".aida";
  }


  const map<string, BinEdges> getBinEdges(string papername) {
    // Get filename
    const string xmlfile = getDataPath(papername);

    // Open AIDA XML file
    TiXmlDocument doc(xmlfile);
    doc.LoadFile();
    if (doc.Error()) {
      string err = "Error in " + string(doc.Value());
      err += ": " + string(doc.ErrorDesc());
      throw runtime_error(err);
    }

    // Return value, to be populated
    map<string, BinEdges> plotbinedges;

    try {      
      // Walk down tree to get to the <paper> element
      const TiXmlNode* aidaN = doc.FirstChild("aida");
      if (!aidaN) throw runtime_error("Couldn't get <aida> root element");
      for (const TiXmlNode* dpsN = aidaN->FirstChild("dataPointSet"); dpsN; dpsN = dpsN->NextSibling()) {
        /// @todo Check AIDA path for "^/HepData" to make sure that this is a reference histogram.
        const TiXmlElement* dpsE = dpsN->ToElement();
        const string plotname = dpsE->Attribute("name");
        list<double> edges;
        for (const TiXmlNode* dpN = dpsN->FirstChild("dataPoint"); dpN; dpN = dpN->NextSibling()) {
          const TiXmlNode* xMeasN = dpN->FirstChild("measurement");
          if (xMeasN) {
            const TiXmlElement* xMeasE = xMeasN->ToElement();
            const string centreStr = xMeasE->Attribute("value");
            const string errplusStr = xMeasE->Attribute("errorPlus"); 
            const string errminusStr = xMeasE->Attribute("errorMinus"); 
            //if (!centreStr) throw runtime_error("Couldn't get a valid bin centre");
            //if (!errplusStr) throw runtime_error("Couldn't get a valid bin err+");
            //if (!errminusStr) throw runtime_error("Couldn't get a valid bin err-");
            istringstream ssC(centreStr);
            istringstream ssP(errplusStr);
            istringstream ssM(errminusStr);
            double centre, errplus, errminus;
            ssC >> centre; ssP >> errplus; ssM >> errminus;
            cout << "  " << centre << " + " << errplus << " - " << errminus << endl;
            const double lowedge = centre - errminus;
            const double highedge = centre + errplus;
            edges.push_back(lowedge);
            edges.push_back(highedge);
          } else {
            cerr << "Couldn't get <measurement> tag" << endl;
            /// @todo Throw an exception here?
          }
        }

        //cout << edges.size() << " edges -> " << edges.size()/2 << " bins" << endl;

        // Remove duplicates (the careful testing is why we haven't used a set)
        for (list<double>::iterator e = edges.begin(); e != edges.end(); ++e) {
          list<double>::iterator e2 = e;
          while (e2 != edges.end()) {
            if (e != e2) {
              if (fuzzyEquals(*e, *e2)) {
                edges.erase(e2++);
              }
            }
            ++e2;
          }
        }

        //cout << edges.size() << " edges after dups removal (should be #bins+1)" << endl;

        // Add to the map
        //BinEdges edgesvec = 
        plotbinedges[plotname] = BinEdges(edges.begin(), edges.end()); 
      }

    }
    // Write out the error
    /// @todo Rethrow as a general XML failure. 
    catch (exception& e) {
      cerr << e.what() << endl;
      throw;
    }

    // Return
    return plotbinedges;
  }

}
