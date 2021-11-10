ARG RIVET_VERSION
FROM hepstore/rivet:${RIVET_VERSION}
LABEL maintainer="rivet-developers@cern.ch"

WORKDIR /work

ARG MG5_URL
RUN wget --no-verbose $MG5_URL -O- | tar xz

RUN true \
  && apt-get update -y && apt-get install -y gnuplot-nox && apt-get autoclean -y \
  && cd MG5_aMC_v* \
  && echo -e "set auto_update 0\ninstall pythia8\n" > mg5setup.txt \
  && bin/mg5_aMC mg5setup.txt \
  && rm mg5setup.txt

RUN lhapdf install NNPDF23_lo_as_0130_qed

WORKDIR /work
