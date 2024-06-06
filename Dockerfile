FROM ubuntu:22.04
#  Install miniconda
RUN  apt-get update && apt-get install -y wget
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py39_23.9.0-0-Linux-x86_64.sh -O ~/miniconda.sh && \
  /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=/opt/conda/bin:${PATH}

COPY environment.yml .
RUN \
   conda install -c conda-forge mamba \ 
   && mamba env update -n root -f environment.yml \
   && conda clean -a

# Install software
RUN apt-get update && \
    apt-get install -y gfortran \
    python3 \
    zlib1g-dev \
    libgomp1 \
    procps \
    libx11-6 \
    bc
    
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Install mutserve (not as conda package available)

ENV MUTSERVE_VERSION=2.0.1
RUN mkdir /opt/mutserve
WORKDIR "/opt/mutserve"
RUN wget https://github.com/seppinho/mutserve/releases/download/v${MUTSERVE_VERSION}/mutserve.zip
#COPY files/mutserve.zip .
RUN unzip mutserve.zip && \
    rm mutserve.zip
ENV PATH="/opt/mutserve:${PATH}"


# Install haplocheck 1.3.3

ENV HAPLOCHECK_VERSION=1.3.3
RUN mkdir /opt/haplocheck
WORKDIR "/opt/haplocheck"
RUN wget https://github.com/genepi/haplocheck/releases/download/v${HAPLOCHECK_VERSION}/haplocheck.zip
RUN unzip haplocheck.zip && \
    rm haplocheck.zip
ENV PATH="/opt/haplocheck:${PATH}"

ENV HAPLOGREP_VERSION=3.2.2
RUN mkdir /opt/haplogrep
WORKDIR "/opt/haplogrep"
RUN wget https://github.com/genepi/haplogrep3/releases/download/v${HAPLOGREP_VERSION}/haplogrep3-${HAPLOGREP_VERSION}-linux.zip && \
    unzip haplogrep3-${HAPLOGREP_VERSION}-linux.zip && \
    rm haplogrep3-${HAPLOGREP_VERSION}-linux.zip && \
    ./haplogrep3 install-tree phylotree-fu-rcrs@1.2
ENV PATH="/opt/haplogrep:${PATH}"

WORKDIR "/opt"
RUN wget https://github.com/jbangdev/jbang/releases/download/v0.91.0/jbang-0.91.0.zip && \
    unzip -q jbang-*.zip && \
    mv jbang-0.91.0 jbang  && \
    rm jbang*.zip

ENV PATH="/opt/jbang/bin:${PATH}"
WORKDIR "/opt"
COPY ./bin/VariantMerger.java ./
RUN jbang export portable -O=VariantMerger.jar VariantMerger.java

COPY ./bin/CoverageEstimation.java ./
RUN jbang export portable -O=CoverageEstimation.jar CoverageEstimation.java


