FROM fedora:27
LABEL maintainer="rivet@projects.hepforge.org"
#RUN dnf update -y
RUN dnf install -y make gcc-c++ gcc-gfortran redhat-rpm-config wget tar less bzip2 findutils python python-devel file python-matplotlib gsl-devel
RUN dnf install -y texlive-latex-bin texlive-texconfig-bin texlive-pst-tools which ghostscript ImageMagick texlive-dvips texlive-relsize
RUN dnf install -y nano zlib-devel
RUN dnf install -y texlive-cm texlive-hyphen-base
RUN dnf install -y texlive-collection-fontsrecommended
RUN dnf clean all
RUN mkdir /code && cd /code && \
    wget https://phab.hepforge.org/source/rivetbootstraphg/browse/2.7.0/rivet-bootstrap?view=raw -O rivet-bootstrap && \
    chmod +x rivet-bootstrap && \
    INSTALL_PREFIX=/usr/local INSTALL_GSL=0 MAKE="make -j7" ./rivet-bootstrap
RUN rm -rf /code
ENV LD_LIBRARY_PATH /usr/local/lib
ENV PYTHONPATH /usr/local/lib64/python2.7/site-packages
RUN texconfig rehash
#COPY Rivet.yoda /tmp/Rivet.yoda
#RUN rivet-mkhtml /tmp/Rivet.yoda
#RUN ln -s /usr/local/share/Rivet/rivet-completion /etc/bash_completion.d/rivet-completion.sh
#RUN ln -s /usr/local/share/YODA/yoda-completion /etc/bash_completion.d/yoda-completion.sh
RUN echo "source /usr/local/share/Rivet/rivet-completion" > /etc/profile.d/rivet-completion.sh
RUN echo "source /usr/local/share/YODA/yoda-completion" > /etc/profile.d/yoda-completion.sh
