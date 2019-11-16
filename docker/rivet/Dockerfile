FROM fedora:30
LABEL maintainer="rivet@projects.hepforge.org"

RUN dnf update -y \
    && dnf install -y \
      make gcc-c++ gcc-gfortran redhat-rpm-config \
      wget tar less bzip2 findutils which nano zlib-devel \
      python python-devel python-requests file python-matplotlib gsl-devel \
      texlive-latex-bin texlive-texconfig-bin texlive-pst-tools texlive-pst-arrow \
      ghostscript ImageMagick texlive-dvips texlive-relsize \
      texlive-cm texlive-hyphen-base texlive-collection-fontsrecommended \
    && dnf clean all

RUN mkdir /code && cd /code \
    && wget https://phab.hepforge.org/source/rivetbootstraphg/browse/3.0.2/rivet-bootstrap?view=raw -O rivet-bootstrap \
    && chmod +x rivet-bootstrap \
    && INSTALL_PREFIX=/usr/local INSTALL_GSL=0 INSTALL_RIVETDEV=0 MAKE="make -j7" ./rivet-bootstrap \
    && echo "source /usr/local/share/Rivet/rivet-completion" > /etc/profile.d/rivet-completion.sh \
    && echo "source /usr/local/share/YODA/yoda-completion" > /etc/profile.d/yoda-completion.sh \
    && texconfig rehash \
    && rm -rf /code

ENV LD_LIBRARY_PATH /usr/local/lib
ENV PYTHONPATH /usr/local/lib64/python2.7/site-packages

WORKDIR /work
