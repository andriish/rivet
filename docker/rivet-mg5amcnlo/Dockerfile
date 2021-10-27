ARG RIVET_VERSION
FROM hepstore/rivet:${RIVET_VERSION}-ubuntu-gcc-hepmc2-py3
LABEL maintainer="rivet-developers@cern.ch"

WORKDIR /work

ARG MG5_URL
RUN wget --no-verbose $MG5_URL -O- | tar xz
#RUN wget --no-verbose https://launchpad.net/mg5amcnlo/3.0/3.1.x/+download/MG5_aMC_v3.1.0.tar.gz -O- | tar xz
#RUN wget --no-verbose https://launchpad.net/mg5amcnlo/2.0/2.8.x/+download/MG5_aMC_v2.8.1.tar.gz -O- | tar xz
#RUN wget --no-verbose https://launchpad.net/mg5amcnlo/2.0/2.7.x/+download/MG5_aMC_v2.7.3.tar.gz -O- | tar xz

RUN true \
  && apt-get update -y && apt-get install -y gnuplot-nox && apt-get autoclean -y \
  && cd MG5_aMC_v* \
  && echo -e "set auto_update 0\ninstall pythia8\n" > mg5setup.txt \
  && bin/mg5_aMC mg5setup.txt \
  && rm mg5setup.txt

RUN lhapdf install NNPDF23_lo_as_0130_qed

WORKDIR /work
