#include <iostream>
#include <string>
#include <list>
#include <sstream>
#include <stdexcept>

#include "TinyXML/tinyxml.h"

using namespace std;


int main(int argc, char* argv[]) {
  string xmlpath = "test.xml";
  if (argc > 1) xmlpath = argv[1];

  //popen("ls -l");

  TiXmlDocument doc(xmlpath);
  doc.LoadFile();
  if ( doc.Error() ) {
    cerr << "Error in " << doc.Value() << ": " << doc.ErrorDesc() << endl;
    return EXIT_FAILURE;
  }


  cout << "Getting bin edges from document..." << endl;
  try {

    // Walk down tree to get to the <paper> element
    const TiXmlNode* aidaN = doc.FirstChild("aida");
    if (!aidaN) throw runtime_error("Couldn't get <aida> root element");
    for (const TiXmlNode* dpsN = aidaN->FirstChild("dataPointSet"); dpsN; dpsN = dpsN->NextSibling()) {
      const TiXmlElement* dpsE = dpsN->ToElement();
      cout << "DataPointSet: " << dpsE->Attribute("name") << endl;
      list<double> edges;
      for (const TiXmlNode* dpN = dpsN->FirstChild("dataPoint"); dpN; dpN = dpN->NextSibling()) {
        const TiXmlNode* xMeasN = dpN->FirstChild("measurement");
        if (xMeasN) {
          const TiXmlElement* xMeasE = xMeasN->ToElement();
          const string centreStr = xMeasE->Attribute("value");
          const string errplusStr = xMeasE->Attribute("errorPlus"); 
          const string errminusStr = xMeasE->Attribute("errorMinus"); 
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
        }
      }

      // Print out the edges
      cout << "Raw: [";
      for (list<double>::const_iterator e = edges.begin(); e != edges.end(); ++e) {
        cout << *e << " ";
      }
      cout << "]" << endl;

      // Remove duplicates(within fractional 10^-3)
      for (list<double>::iterator e = edges.begin(); e != edges.end(); ++e) {
        list<double>::iterator e2 = e;
        while (e2 != edges.end()) {
          if (e != e2 && *e == *e2) {
            edges.erase(e2++);
          } else {
            ++e2;
          }
        }
      }

      // Print out the edges
      cout << "Stripped: [";
      for (list<double>::const_iterator e = edges.begin(); e != edges.end(); ++e) {
        cout << *e << " ";
      }
      cout << "]" << endl;
      cout << endl;
    }

  } 
  catch (exception& e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
