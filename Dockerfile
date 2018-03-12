FROM ubuntu:latest

ENV GFDL_WORK /tmp
ENV GFDL_BASE /isca
ENV GFDL_DATA /data
ENV GFDL_ENV docker

COPY . /isca



RUN useradd -ms /bin/bash isca
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y \
    git python3 python3-pip \
    libnetcdf-dev libpnetcdf-dev libnetcdff-dev \
    libhdf5-openmpi-dev wget tcl tcl-dev \
    build-essential gcc gfortran

RUN wget https://sourceforge.net/projects/modules/files/Modules/modules-4.1.1/modules-4.1.1.tar.gz/download
RUN tar xvf download

RUN cd modules-4.1.1 && ./configure && make && make install
RUN /bin/bash -c "source /usr/local/Modules/init/profile.sh"
RUN pip3 install -r /isca/src/extra/python/requirements.txt
RUN pip3 install -e /isca/src/extra/python


RUN chown -R isca /isca
RUN mkdir -p /data && chown -R isca /data

WORKDIR /isca
VOLUME /data

USER isca
RUN python3 /isca/exp/test_cases/held_suarez/held_suarez_test_case.py 