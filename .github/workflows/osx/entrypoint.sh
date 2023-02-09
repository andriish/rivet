#!/bin/sh -l
set -x
brew install wget coreutils gcc autoconf automake
brew tap davidchall/hep
brew install hepmc3 yoda fastjet fjcontrib 
TOP=$(pwd)

autoreconf --force --install --verbose .
automake -a --force
configure  --prefix=$TOP/INSTALL --disable-doxygen --with-yoda=$(yoda-config --prefix ) --with-hepmc3=$(HepMC3-config --prefix) --with-fjcontrib=$(fastjet-config --prefix) --with-fastjet=$(fastjet-config --prefix)
make -j 4
make install
find $TOP/INSTALL
