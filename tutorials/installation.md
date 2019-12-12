## Installation

The easiest way to start using Rivet is via the Docker system (see https://www.docker.com/), which is like a lightweight Linux virtual machine that's easily installed and run on Mac, Linux and Windows machines. See the dedicated [Docker](docker.md) tutorial page for more detail.

The rest of these instructions are mainly aimed at users who want to natively install and run a release of Rivet on their machine.

GettingStartedForDevelopers has some additional or replacement steps for people wishing to check out the development version from the repository and build from there.

# Native installation

The simplest Rivet installation uses a "bootstrap" script to install Rivet and all its dependencies from release tarballs.

If you are installing Rivet on an Apple Mac, please take a look at the MacInstallationIssues first. If you want to use packages from CERN's CVMFS network file system to assist your installation, or to directly use Rivet from there, see the bottom of this page.

# Installation of Rivet and all dependencies

'''Prerequisite''' Python header files are required. On Ubuntu,
you can use this command to install the necessary files system-wide:

```
sudo apt-get install python-dev
```


1. Download the bootstrap script''' into a temporary working directory, and make it executable:
```
  cd /scratch/rivet
  wget https://phab.hepforge.org/source/rivetbootstraphg/browse/3.0.2/rivet-bootstrap?view=raw -O rivet-bootstrap
  chmod +x rivet-bootstrap
```
(Replace the version string as appropriate if you want to install other versions of Rivet.)

2. Check the options.''' Look at the header of the script to see all variables which you can set, e.g. to skip installation of certain dependencies if they are available in your system:
```
  less rivet-bootstrap ## and read...
```


3. Run the script.''' By default it will install to `$PWD/local`, where `$PWD` is the current directory. If you need to change that, specify the corresponding values on the command line. Examples:
```
./rivet-bootstrap
  # or, e.g.
INSTALL_PREFIX=$HOME/software/rivet MAKE="make -j8" ./rivet-bootstrap
```
We will refer to the installation root path as `$PREFIX`.

# Setting up the environment

After the script grinds away for a while, it will tell you that it is finished and how to set up a runtime environment (similar to that used inside the build script) for running Rivet. A sourceable rivetenv.(c)sh script is provided for (c)sh shell users to help set up this environment. Here's how to set up the environment and then test the `rivet` program's help feature and analysis listing:

```
  source $PREFIX/rivetenv.sh
  rivet --help
  rivet --list-analyses
```

If that works, everything is installed correctly. If you are using the `bash` shell in your terminal, then Rivet will offer you programmable tab completion: try typing {{{rivet}}} and pressing the Tab key!

'''You may wish to add the environment variable settings to your `~/.bashrc` shell config file, so that Rivet will work without needing any special session setup.'''

You can now check out the FirstRivetRun guide.

