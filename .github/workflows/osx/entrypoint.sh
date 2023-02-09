#!/bin/sh -l
set -x
brew install wget coreutils gcc autoconf automake cython libtool gnu-which
export PATH=/usr/local/opt/cython/bin:$PATH
cp /usr/local/opt/cython/bin/cython /usr/local/opt/cython/bin/cython2
cp /usr/local/opt/cython/bin/cython /usr/local/opt/cython/bin/cython-2
cp /usr/local/opt/cython/bin/cython /usr/local/opt/cython/bin/Cython-2
cp /usr/local/opt/cython/bin/cython /usr/local/opt/cython/bin/Cython2
ls -lah /usr/local/opt/cython/bin
which -a cython
which -a cython2
which -a cython-2
which -a Cython-2
which -a Cython2

brew tap davidchall/hep
brew install hepmc3 yoda fastjet 
TOP=$(pwd)

mkdir FJ
cd FJ
wget https://fastjet.hepforge.org/contrib/downloads/fjcontrib-1.048.tar.gz
tar zxfv fjcontrib-1.048.tar.gz
cd fjcontrib-1.048
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


