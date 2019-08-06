FROM centos:centos7

RUN yum install -y gcc gcc-c++ \
                   libtool libtool-ltdl \
                   make cmake \
                   git \
                   pkgconfig \
                   sudo \
                   automake autoconf \
                   yum-utils rpm-build && \
    yum clean all

RUN yum install -y zlib-devel glibc-static libstdc++-static

ENV FLAVOR=rpmbuild OS=centos DIST=el7

RUN mkdir -p /libs

COPY libs/ /libs

COPY Makefile_assembly /

COPY Makefile_scaffold /

RUN sed -i 's/libs/\/libs/g' Makefile_assembly
RUN sed -i 's/libs/\/libs/g' Makefile_scaffold

WORKDIR /in 

ENV LD_LIBRARY_PATH='/libs:${LD_LIBRARY_PATH}'

#COPY "docker_make.sh" /usr/local/bin/

#ENTRYPOINT ["/usr/local/bin/docker_make.sh"]
