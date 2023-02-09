#!/bin/sh -l
set -x
brew install wget coreutils gcc autoconf automake cython libtool
brew tap davidchall/hep
brew install hepmc3 yoda fastjet 
TOP=$(pwd)

mkdir FJ
cd FJ
wget https://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.048.tar.gz
tar zxfv fjcontrib-1.048.tar.gz
./configure --prefix=$(fastjet-config --prefix)
make fragile-shared-install
make install


cd $TOP
autoreconf --force --install --verbose .
automake -a --force
./configure  --prefix=$TOP/INSTALL --disable-doxygen --with-yoda=$(yoda-config --prefix ) --with-hepmc3=$(HepMC3-config --prefix) --with-fjcontrib=$(fastjet-config --prefix) --with-fastjet=$(fastjet-config --prefix)
make -j 4
make install
find $TOP/INSTALL


