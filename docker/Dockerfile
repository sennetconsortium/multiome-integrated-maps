FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    flex \
    git \
    libfreetype6-dev \
    libigraph-dev \
    libxml2-dev \
    libtool \
    m4 \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    zlib1g-dev \
    curl \
    unzip \
 && rm -rf /var/cache/apt/*

# Install AWS CLI using the official AWS method (via curl and unzip)
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf aws awscliv2.zip

RUN python3 -m venv /opt/venv \
 && /opt/venv/bin/pip install --upgrade pip
COPY docker/requirements.txt /opt/requirements.txt
RUN /opt/venv/bin/pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

COPY bin /opt
COPY data /opt/data
RUN chmod -R +x /opt

WORKDIR /opt

ENV PATH="/opt/venv/bin:$PATH"
