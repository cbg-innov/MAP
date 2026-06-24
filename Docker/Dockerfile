FROM mambaorg/micromamba:2.5-debian13-slim

ARG TARGETARCH
ENV TARGETARCH=${TARGETARCH}
ENV MAMBA_ROOT_PREFIX=/home/mamba
ENV ENV_NAME=map

#set working directory
WORKDIR /MAP

# Copy files from this directory into image
COPY . .

#mamba environment
RUN micromamba create -y -n map -c conda-forge -c bioconda -c anaconda \
    python=3.12 \
    snakemake-minimal=9.16.3 \
    r-base=4.5.2 \
    pear=0.9.6 \
    seqkit=2.12.0 \
    libgcc=15.2.0 \
    cutadapt=5.2 \
    seqtk=1.5 \
    vsearch=2.30.4 \
    parallel=20260122 \
    rename=1.601 \
    bc=1.07.1 \
    zip=3.0 \
    chopper=0.13.0 \
# and python libs.
    pandas=3.0.0 \
    openpyxl=3.1.5 \
    biopython=1.86 \
# and all the Rpkgs
    r-tidyr=1.3.2 \
    r-data.table=1.17.8 \
    bioconductor-biostrings=2.78 \
    r-scales=1.4.0 \
    r-ggplot2=4.0.2 \
    r-patchwork=1.3.2 \
    r-openxlsx=4.2.8.1 \
    r-plyr=1.8.9 \
    r-readxl=1.4.5 \
    r-plotly=4.12.0 \
    r-forcats=1.0.1 \
    r-shiny=1.13.0 \
    r-kableextra=1.4.0 \
    r-treemapify=2.6.0 \
    r-paletteer=1.7.0 \
    r-htmltools=0.5.9 \
    r-dt=0.34.0 \
    r-leaflet=2.2.3 \
    && micromamba clean --all --yes

USER root

RUN apt update && apt install -y curl && rm -rf /var/lib/apt/lists/* && \
    if [ "$TARGETARCH" = "amd64" ]; then \
        echo "Installing Quarto AMD64"; \
        curl -L https://github.com/quarto-dev/quarto-cli/releases/download/v1.8.27/quarto-1.8.27-linux-amd64.deb \
            -o /tmp/quarto.deb && \
        dpkg -i /tmp/quarto.deb && rm /tmp/quarto.deb; \
    elif [ "$TARGETARCH" = "arm64" ]; then \
        echo "Installing Quarto ARM64"; \
        curl -L https://github.com/quarto-dev/quarto-cli/releases/download/v1.8.27/quarto-1.8.27-linux-arm64.deb \
            -o /tmp/quarto.deb && \
        dpkg -i /tmp/quarto.deb && rm /tmp/quarto.deb; \
    else \
        echo "Unsupported architecture: $TARGETARCH" && exit 1; \
    fi

#Activate Map environment.
RUN echo "micromamba activate map" >> ~/.bashrc

#Make directory structure
RUN mkdir SCRIPTS
RUN mkdir Metabarcoding

#Move files to correct place
RUN mv *.sh SCRIPTS
RUN mv *.py SCRIPTS
RUN mv *.R SCRIPTS
RUN mv *.qmd SCRIPTS
RUN mv *.fastq.gz Metabarcoding
RUN mv *.xlsx Metabarcoding

# Ensure Map environment when image is run.
SHELL ["micromamba", "run", "-n", "map", "/bin/bash", "-c"]

# iNEXT installation doesn't work well with mamba. Install it from CRAN.
RUN R -e "install.packages('iNEXT', repos='https://cloud.r-project.org', type='source')"

#Download and unzip latest BOLDDistilled lib (curl needs map environment running, hence placement here at the end)
RUN curl https://us-sea-1.linodeobjects.com/boldistilled/sintax.zip --output reflib_latest.zip
RUN python -m zipfile -e reflib_latest.zip .

#change dir name to REFS, remove zip file.
RUN mv sintax/ REFS/
RUN python -m zipfile -e reference_seqs_327K.fasta.zip ./REFS/
RUN rm reflib_latest.zip reference_seqs_327K.fasta.zip