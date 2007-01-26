#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <string>

#if defined(XERCES_NEW_IOSTREAMS)
#include <iostream>
#else
#include <iostream.h>
#endif

XERCES_CPP_NAMESPACE_USE
using namespace std;


int main (int argc, char* args[]) {
    try {
        XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Error during initialization! :\n" << message << endl;
        XMLString::release(&message);
        return 1;
    }

    XercesDOMParser* parser = new XercesDOMParser();
    parser->setValidationScheme(XercesDOMParser::Val_Always); // optional
    parser->setDoNamespaces(true); // optional

    ErrorHandler* errHandler = (ErrorHandler*) new HandlerBase();
    parser->setErrorHandler(errHandler);

    string xmlFile = "test.xml";
    if (argc == 2) xmlFile = args[1];
    try {
        parser->parse(xmlFile.c_str());
    }
    catch (const XMLException& toCatch) {
        char* message = XMLString::transcode(toCatch.getMessage());
        cout << "Exception message is: \n" << message << endl;
        XMLString::release(&message);
        return -1;
    }
    catch (const DOMException& toCatch) {
        char* message = XMLString::transcode(toCatch.msg);
        cout << "Exception message is: \n" << message << endl;
        XMLString::release(&message);
        return -1;
    }
    catch (...) {
      cout << "Unexpected Exception - does " << xmlFile << " exist?" << endl;
        return -1;
    }

    DOMNode* doc = parser->getDocument();
    if (doc) {
      DOMNodeList* nodes = doc->getChildNodes();
      if (nodes) {
        cout << "Num nodes = " << nodes->getLength() << endl;
        for (XMLSize_t i = 0; i < nodes->getLength(); ++i) {
          DOMNode* node = nodes->item(i);
          if (node) {
            cout << i << ": " << *( node->getNodeName() ) << endl;
            cout << "Children: " << node->getChildNodes()->getLength() << endl;
            for (XMLSize_t j = 0; i < node->getChildNodes()->getLength(); ++i) {
              DOMNode* node2 = node->getChildNodes()->item(j);
              if (node2) cout << j << ": " << *( node2->getNodeValue() ) << endl;
            }
          }
        }
      }
    }

    delete parser;
    delete errHandler;
    return EXIT_SUCCESS;
}

// #include <xercesc/util/PlatformUtils.hpp>
// #include <xercesc/dom/DOM.hpp>

// XERCES_CPP_NAMESPACE_USE 
  
// int main(int argc, char* argv[])
// {
//   try {
//     XMLPlatformUtils::Initialize();
//   }
//   catch (const XMLException& toCatch) {
//     // Do your failure processing here
//     return 1;
//   }

//   // Do your actual work with Xerces-C++ here.

//   XMLPlatformUtils::Terminate();

//   // Other terminations and cleanup.
//   return EXIT_SUCCESS;
// }

