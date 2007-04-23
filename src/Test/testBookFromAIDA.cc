#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"

using namespace Rivet;

int main(int argc, char* argv[]) {
  string paper = "HepEx0112029";
  if (argc > 1) paper = argv[1];

  cout << "Path to data file: " << getDataPath(paper) << endl;
  const map<string, BinEdges> datasets = getBinEdges(paper);

  // Print out the edges
  for (map<string, BinEdges>::const_iterator d = datasets.begin(); d != datasets.end(); ++d) {
    cout << "Dataset " << d->first << ": ";
    cout << "[";
    for (vector<double>::const_iterator e = d->second.begin(); e != d->second.end(); ++e) {
      cout << *e << " ";
    }
    cout << "]" << endl;
  }

  return EXIT_SUCCESS;
}
