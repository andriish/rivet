
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

    // Walk down tree to get to the <paper> element
    const TiXmlNode* hepml = doc.FirstChild("hepml");
    if (!hepml) throw runtime_error("Couldn't get <hepml> root element");
    const TiXmlNode* data = hepml->FirstChild("data");
    if (!data) throw runtime_error("Couldn't get <data> element");
    const TiXmlNode* paper = data->FirstChild("paper")->ToElement();
    if (!paper) throw runtime_error("Couldn't get <paper> element");
    const TiXmlElement* p = paper->ToElement();
    cout << "SPIRES ID: " << p->Attribute("irn") << endl;

    // Walk down tree to get to the <bins> element of the first dataset's first x-axis
    const TiXmlNode* ds = paper->FirstChild("dataset");
    if (!ds) throw runtime_error("Couldn't get <dataset> element");
    const TiXmlNode* xaxis = ds->FirstChild("xaxis");
    if (!xaxis) throw runtime_error("Couldn't get <xaxis> element");
    const TiXmlNode* bins = xaxis->FirstChild("bins");
    if (!bins) throw runtime_error("Couldn't get <bins> element");

    // Get the bin edges for each <bin>
    for( const TiXmlNode* bin = bins->FirstChild("bin"); bin; bin = bin->NextSibling()) {
      const TiXmlElement* b = bin->ToElement();
      double low, high;
      const char* lowstr = b->Attribute("low");
      const char* highstr = b->Attribute("high");
      if (!lowstr) throw runtime_error("Couldn't get a valid 'low' attribute on <bin>");
      if (!highstr) throw runtime_error("Couldn't get a valid 'high' attribute on <bin>");
      istringstream ss1(lowstr);
      istringstream ss2(highstr);
      ss1 >> low; ss2 >> high;
      cout << low << "-" << high << endl;
    }

  } 
  catch (exception& e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
