FROM hepstore/rivet:3.1.2-ubuntu-gcc-hepmc2-py2
LABEL maintainer="rivet-developers@cern.ch"

WORKDIR /work

# TODO: switch to Py3 base image at the same time as 2.8.x or 3.0.x
#RUN wget --no-verbose https://launchpad.net/mg5amcnlo/2.0/2.8.x/+download/MG5_aMC_v2.8.1.tar.gz -O- | tar xz
RUN wget --no-verbose https://launchpad.net/mg5amcnlo/2.0/2.7.x/+download/MG5_aMC_v2.7.3.tar.gz -O- | tar xz

WORKDIR /work
