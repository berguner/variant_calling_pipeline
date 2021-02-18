FROM ubuntu:bionic as preprocess

RUN apt-get update && apt-get install -y \
	g++ \
	gcc \
	make \
	autoconf \
	wget \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	parallel \
	pigz \
	alien \
	zip

# Compile BWA
WORKDIR /root/source
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
	&& tar -xf bwa-0.7.17.tar.bz2
WORKDIR /root/source/bwa-0.7.17
RUN make

# Compile bcftools
WORKDIR /root/source
RUN wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 \
	&& tar -xf bcftools-1.11.tar.bz2 
WORKDIR /root/source/bcftools-1.11
RUN ./configure \
	&& make

# Compile samtools
WORKDIR /root/source
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 \
	&& tar -xf samtools-1.11.tar.bz2
WORKDIR /root/source/samtools-1.11
RUN ./configure  --without-curses \
	&& make

# Compile samblaster
WORKDIR /root/source
RUN wget https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.24/samblaster-v.0.1.24.tar.gz \
	&& tar -xf samblaster-v.0.1.24.tar.gz
WORKDIR /root/source/samblaster-v.0.1.24
RUN make

# Compile fastp
WORKDIR /root/source
RUN wget https://github.com/OpenGene/fastp/archive/v0.20.1.tar.gz \
	&& tar -xf v0.20.1.tar.gz
WORKDIR /root/source/fastp-0.20.1
RUN make


# Start building the final container

FROM broadinstitute/gatk:4.1.9.0

################## METADATA ######################
LABEL about.summary="Collection of tools for the BSF variant_calling_pipeline"
LABEL about.home="https://www.biomedical-sequencing.org"
LABEL about.tags="Genomics"
LABEL about.maintainer="Bekir Erguener @ CeMM/BSF"
################## MAINTAINER ######################

COPY --from=preprocess /root/source/bwa-0.7.17/bwa /usr/local/bin/bwa
COPY --from=preprocess /root/source/bcftools-1.11/bcftools /usr/local/bin/bcftools
COPY --from=preprocess /root/source/samtools-1.11/samtools /usr/local/bin/samtools
COPY --from=preprocess /root/source/samblaster-v.0.1.24/samblaster /usr/local/bin/samblaster
COPY --from=preprocess /root/source/fastp-0.20.1/fastp /usr/local/bin/fastp
COPY --from=preprocess /usr/bin/parallel /usr/bin/parallel

CMD ["/bin/bash"]

