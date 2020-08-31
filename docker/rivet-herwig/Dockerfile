FROM ubuntu:20.04
LABEL maintainer="rivet-developers@cern.ch"
SHELL ["/bin/bash", "-c"]

ARG RIVET_VERSION=3.1.2
ARG YODA_VERSION=1.8.3
ARG THEPEG_VERSION=2.2.1
ARG HERWIG_VERSION=7.2.1

ARG CXX_CMD=g++
ARG CC_CMD=gcc
ARG FC_CMD=gfortran

ENV CXX=${CXX_CMD}
ENV CC=${CC_CMD}
ENV FC=${FC_CMD}

RUN export DEBIAN_FRONTEND=noninteractive \
    && if test "$CXX_CMD" = "g++"; then CXX_PKG=g++; else CXX_PKG=clang; fi \
    && if test "$CC_CMD" = "gcc"; then CC_PKG=gcc; else CC_PKG=clang; fi \
    && if test "$FC_CMD" = "gfortran"; then FC_PKG=gfortran; else FC_PKG=flang; fi \
    && apt-get update -y \
    && apt-get install -y ${CXX_PKG} ${CC_PKG} ${FC_PKG} \
    && if test "$CXX_CMD" = "clang++"; then update-alternatives --install /usr/bin/g++ g++ /usr/bin/clang++ 2; fi \
    && if test "$CXX_CMD" = "clang++"; then update-alternatives --install /usr/bin/c++ c++ /usr/bin/clang++ 2; fi \
    && if test "$CC_CMD" = "clang"; then update-alternatives --install /usr/bin/gcc gcc /usr/bin/clang 2; fi \
    && if test "$CC_CMD" = "clang"; then update-alternatives --install /usr/bin/cc cc /usr/bin/clang 2; fi \
    && if test "$FC_CMD" = "flang"; then update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/flang 2; fi \
    && apt-get install -y apt-utils tzdata \
    && apt-get install -y \
      make automake autoconf libtool cmake rsync \
      git wget tar less bzip2 findutils nano file \
      zlib1g-dev \
    && apt-get install -y python3 python3-dev python3-venv python3-wheel \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3 2 \
    && update-alternatives --install /usr/bin/python-config python-config /usr/bin/python3-config 2 \
    && wget --no-verbose https://bootstrap.pypa.io/get-pip.py -O get-pip.py \
    && python get-pip.py \
    && apt-get -y autoclean

RUN export DEBIAN_FRONTEND=noninteractive \
    && python -m venv /herwig \
    && . /herwig/bin/activate \
    && pip install Cython six \
    && wget https://herwig.hepforge.org/downloads/herwig-bootstrap \
    && chmod +x herwig-bootstrap \
    && ./herwig-bootstrap /herwig -j$(nproc --ignore=1) \
      --thepeg-version=${THEPEG_VERSION} --herwig-version=${HERWIG_VERSION} \
      --rivet-version=${RIVET_VERSION} --yoda-version=${YODA_VERSION} \
      --without-hjets \
    && mv herwig-bootstrap /herwig/src

RUN echo ". /herwig/bin/herwig-activate" >> /etc/bash.bashrc

WORKDIR /work
CMD /bin/bash
