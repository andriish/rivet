
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>

#include "TinyXML/tinyxml.h"

using namespace std;


int main(int argc, char* argv[]) {
  string xmlpath = "test.xml";
  if (argc > 1) xmlpath = argv[1];


  TiXmlDocument doc(xmlpath);
  doc.LoadFile();
  if ( doc.Error() ) {
    cerr << "Error in " << doc.Value() << ": " << doc.ErrorDesc() << endl;
    return EXIT_FAILURE;
  }


  cout << "Getting bin edges from document..." << endl;
  try {
    TiXmlNode* hepml = doc.FirstChild("hepml");
    if (!hepml) throw runtime_error("Couldn't get <hepml> root element");
    TiXmlNode* data = hepml->FirstChild("data");
    if (!data) throw runtime_error("Couldn't get <data> element");
    TiXmlNode* paper = data->FirstChild("paper")->ToElement();
    if (!paper) throw runtime_error("Couldn't get <paper> element");
    TiXmlElement* p = paper->ToElement();
    cout << "SPIRES ID: " << p->Attribute("irn") << endl;

    TiXmlNode* ds = paper->FirstChild("dataset");
    if (!ds) throw runtime_error("Couldn't get <dataset> element");
    TiXmlNode* xaxis = ds->FirstChild("xaxis");
    if (!xaxis) throw runtime_error("Couldn't get <xaxis> element");
    TiXmlNode* bins = xaxis->FirstChild("bins");
    if (!bins) throw runtime_error("Couldn't get <bins> element");

    for( TiXmlNode* bin = bins->FirstChild("bin"); bin; bin = bin->NextSibling()) {
      TiXmlElement* b = bin->ToElement();
      double low, high;
      istringstream ss1(b->Attribute("lowedge"));
      istringstream ss2(b->Attribute("highedge"));
      ss1 >> low; ss2 >> high;
      cout << low << "-" << high << endl;
    }
  } catch (exception& e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  

  return EXIT_SUCCESS;
}
