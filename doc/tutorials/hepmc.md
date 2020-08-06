## Linking against HepMC

Monte Carlo event generators deliver generated events to Rivet in the HepMC (http://hepmc.web.cern.ch/hepmc/) format. While HepMC is currently in release series 3, not all generators and experiments have caught up with these developments. Rivet therefore supports both HepMC2 and HepMC3, but encourage users to integrate HepMC3 in applications.

### Bootstrap installation (HepMC2)

The Rivet [bootstrap script](installation.md) currently installs the latest version of HepMC2. In order to use HepMC3, one must therefore perform a manual installation of Rivet dependencies.

### Manual installation and linking of HepMC3

Rivet is tested against HepMC-3.1.1, later versions might not be supported by the current Rivet release. HepMC-3.1.1 can be downloaded and unpacked:

```
wget http://hepmc.web.cern.ch/hepmc/releases/HepMC3-3.1.1.tar.gz
tar xvf HepMC3-3.1.1.tar.gz
```
Make a build directory for HepMC and configure the build with `cmake`. The installation directory should be set to a user defined location, instead of the dummy location shown here. Note that HepMC *must* have with static libraries built, and should be built with the search function enabled (ie. the default settings should not be changed).

```
mkdir hepmc3-build
cd hepmc3-build
cmake -DHEPMC3_ENABLE_ROOTIO=OFF -DCMAKE_INSTALL_PREFIX=<hepmc installation path> ../HepMC3-3.1.1

```
HepMC3 can then be compiled and installed to the selected location:

```
make -jN && make install
```
We are now ready to configure Rivet (assuming that the other dependencies (YODA, Fastjet and Fastjet contrib) are available in the default path), using the static HepMC library. In the Rivet directory:
```
./configure --with-hepmc3=<hepmc installation path> --with-hepmc3-libname=HepMC3_static
```
And Rivet can then be compiled and installed as usual:
```
make -jN && make install
```
