FROM ubuntu:latest

ENV GFDL_WORK /tmp
ENV GFDL_BASE /isca
ENV GFDL_DATA /data
ENV GFDL_ENV docker

RUN useradd -ms /bin/bash isca
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y \
    git python3 python3-pip \
    libnetcdf-dev libpnetcdf-dev libnetcdff-dev \
    libhdf5-openmpi-dev wget tcl tcl-dev \
    build-essential gcc gfortran

# RUN wget https://sourceforge.net/projects/modules/files/Modules/modules-4.1.1/modules-4.1.1.tar.gz/download
# RUN tar xvf download

# RUN cd modules-4.1.1 && ./configure && make && make install
# RUN /bin/bash -c "source /usr/local/Modules/init/profile.sh"

COPY . /isca
RUN chown -R isca /isca

RUN pip3 install -r /isca/src/extra/python/requirements.txt
RUN pip3 install -e /isca/src/extra/python


RUN mkdir -p /data && chown -R isca /data
VOLUME /data

WORKDIR /isca
USER isca

RUN python3 -c "import isca; cb = isca.IscaCodeBase.from_directory('/isca'); cb.compile()"
#CMD python3 /isca/exp/held_suarez.py --compile --up-to -i 3 -n 2