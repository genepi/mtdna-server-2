FROM nfcore/base:1.14
LABEL Hansi Weissensteiner <hansi.weissensteiner@i-med.ac.at>

COPY environment.yml .
RUN \
   conda env update -n root -f environment.yml \
&& conda clean -a

# RUN apt-get --allow-releaseinfo-change update

# Install mutserve (not as conda package available)

RUN mkdir /opt/mutserve
WORKDIR "/opt/mutserve"
#RUN wget https://github.com/seppinho/mutserve/releases/download/v2.0.0-rc12/mutserve.zip
COPY files/mutserve.zip .
RUN unzip mutserve.zip && \
    rm mutserve.zip
ENV PATH="/opt/mutserve:${PATH}"


# Install haplocheck 1.3.3

RUN mkdir /opt/haplocheck
WORKDIR "/opt/haplocheck"
RUN wget https://github.com/genepi/haplocheck/releases/download/v1.3.3/haplocheck.zip
RUN unzip haplocheck.zip && \
    rm haplocheck.zip
ENV PATH="/opt/haplocheck:${PATH}"

RUN mkdir /opt/haplogrep
WORKDIR "/opt/haplogrep"
RUN wget https://github.com/genepi/haplogrep3/releases/download/v3.2.1/haplogrep3-3.2.1-linux.zip && \
    unzip haplogrep3-3.2.1-linux.zip && \
    rm haplogrep3-3.2.1-linux.zip
ENV PATH="/opt/haplogrep:${PATH}"


