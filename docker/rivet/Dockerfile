FROM hepstore/hepbase-ubuntu-gcc-hepmc2-py3
LABEL maintainer="rivet@projects.hepforge.org"

ENV RIVET_VERSION=3.1.1

RUN apt-get update -y && apt-get install -y \
      texlive-latex-recommended texlive-fonts-recommended \
      graphicsmagick

RUN mkdir /code && cd /code \
    && wget https://gitlab.com/hepcedar/rivetbootstrap/raw/${RIVET_VERSION}/rivet-bootstrap \
    && chmod +x rivet-bootstrap \
    && INSTALL_PREFIX=/usr/local INSTALL_CYTHON=0 \
       INSTALL_RIVET=1 INSTALL_YODA=0 \
       INSTALL_HEPMC=0 HEPMC_VERSION=3.2.1 \
       INSTALL_FASTJET=0 INSTALL_FJCONTRIB=0 \
       MAKE="make -j6" ./rivet-bootstrap \
    && echo "source /usr/local/share/Rivet/rivet-completion" > /etc/profile.d/rivet-completion.sh \
    && echo "source /usr/local/share/YODA/yoda-completion" > /etc/profile.d/yoda-completion.sh \
    && texconfig rehash \
    && rm -rf /code

# ENV LD_LIBRARY_PATH /usr/local/lib
# ENV PYTHONPATH /usr/local/lib64/python2.7/site-packages

WORKDIR /work
