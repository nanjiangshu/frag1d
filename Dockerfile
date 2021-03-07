FROM nanjiang/common-ubuntu
LABEL maintainer "Nanjiang Shu (nanjiang.shu@scilifelab.se)"
LABEL version "1.4"

#================================
# Install basics
#===============================
RUN apt-get update -y &&\
    apt-get install -y apt-utils  \
                       curl wget bc \
                       python python-dev python-pip \
                       build-essential  \
                       make  \
                       locales-all \
                       lib32ncurses5 lib32z1 \
                       blast2

#================================
#  Add frag1d source code
#===============================
WORKDIR /app
# add the source code to WORKDIR /app
ADD frag1d ./frag1d

RUN mkdir -p /scratch/ /static/
# building frag1d
RUN cd /app/frag1d/&& \
    make && make install

#================================
# Setting ENVs
#===============================
ENV USER_DIRS "/app"
ENV BLASTDB /data/blastdb
ENV BLASTBIN /usr/bin

CMD ["/bin/bash" ]
