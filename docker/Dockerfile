FROM continuumio/miniconda3:25.3.1-1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    flex \
    git \
    libfreetype6-dev \
    libxml2-dev \
    libtool \
    m4 \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    zlib1g-dev \
    curl \
 && rm -rf /var/cache/apt/*

## Python packages
COPY docker/environment.yml /opt/environment.yml
RUN conda env create -f /opt/environment.yml

ENV CONDA_ENV_NAME=multi-maps
ENV CONDA_PREFIX=/opt/conda/envs/${CONDA_ENV_NAME}
ENV PATH=${CONDA_PREFIX}/bin:$PATH

RUN pip3 install git+https://github.com/satijalab/panhumanpy.git@d4c69410e7dc23213b7a848c8c1d8274e333272a
RUN curl -L -o /opt/pan-human-azimuth-crosswalk.csv https://cdn.humanatlas.io/digital-objects/ctann/pan-human-azimuth/v1.0/assets/pan-human-azimuth-crosswalk.csv

COPY bin/ /opt
COPY data /opt/data

WORKDIR /opt

CMD ["/bin/bash"]