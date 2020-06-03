FROM hepstore/rivet:3.1.1 AS build
LABEL maintainer="rivet-developers@cern.ch"

FROM registry.fedoraproject.org/fedora-minimal:30

COPY --from=build /usr/local/bin/* /usr/local/bin/
COPY --from=build /usr/local/lib /usr/local/lib/
COPY --from=build /usr/local/lib64/python2.7/site-packages /usr/local/lib64/python2.7/site-packages/
COPY --from=build /usr/local/share /usr/local/share/
COPY --from=build /usr/local/*env*.sh /usr/local/
COPY --from=build /etc/profile.d/*-completion.sh /etc/profile.d/

RUN microdnf install python libgomp && \
    microdnf clean all && \
    rm /usr/local/bin/rivet-build* /usr/local/bin/rivet-mkhtml /usr/local/bin/make-plots* /usr/local/bin/yodaplot

ENV LD_LIBRARY_PATH /usr/local/lib
ENV PYTHONPATH /usr/local/lib64/python2.7/site-packages
WORKDIR /work
