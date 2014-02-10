#ifndef RIVET_Rivet_HH
#define RIVET_Rivet_HH

#include "Rivet/Config/RivetConfig.hh"
#include "Rivet/Config/BuildOptions.hh"

// Bring selected STL classes into the Rivet namespace
#include "Rivet/Tools/RivetSTL.hh"

// Macro to help with overzealous compiler warnings
#ifdef UNUSED
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#else
# define UNUSED(x) x
#endif


namespace Rivet {

  /// A sensible default maximum value of rapidity for Rivet analyses to use.
  static const double MAXRAPIDITY = 100000.0;

  /// A function to get the Rivet version string
  string version();

}

// HepMC headers and helper functions
#include "Rivet/Tools/RivetHepMC.hh"

// Now import some Rivet classes
#include "Rivet/Exceptions.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Math/Matrices.hh"
#include "Rivet/Math/Units.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/RivetPaths.hh"

#include "Rivet/ParticleName.hh"
#include "Rivet/Particle.hh"

#endif
