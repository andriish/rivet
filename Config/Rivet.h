#ifndef RIVET_Rivet_H
#define RIVET_Rivet_H

#include <typeinfo>
#include <set>
#include <string>
#include <vector>
#include <stdexcept>

/**
 * This is the main namespace in which all Rivet classes are defined.
 */
namespace Rivet {

using std::set;
using std::type_info;
using std::string;
using std::less;
using std::vector;
using std::pair;
using std::runtime_error;

const double MaxRapidity = 100000.0;

}

#endif
