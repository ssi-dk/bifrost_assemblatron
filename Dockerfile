FROM \
    ssidk/bifrost-base:2.0.5

LABEL \
    name="bifrost-assemblatron" \
    description="Docker environment for assemblatron in bifrost" \
    version="2.0.5" \
    DBversion="31/07/2019" \
    maintainer="kimn@ssi.dk;"

RUN \
    conda install -yq -c conda-forge -c bioconda -c defaults bbmap==38.58; \
    conda install -yq -c conda-forge -c bioconda -c defaults skesa==2.3.0; \
    conda install -yq -c conda-forge -c bioconda -c defaults minimap2==2.17; \
    conda install -yq -c conda-forge -c bioconda -c defaults samtools==1.9; \
    conda install -yq -c conda-forge -c bioconda -c defaults cyvcf2==0.11.4; \
    conda install -yq -c conda-forge -c bioconda -c defaults prokka==1.14.0; \
    conda install -yq -c conda-forge -c bioconda -c defaults quast==5.0.2; \
    cd /bifrost; \
    git clone https://github.com/ssi-dk/bifrost-assemblatron.git assemblatron; \

ADD \
    https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/adapters.fasta /bifrost_resources/
RUN \
    chmod +r /bifrost_resources/adapters.fasta

ENTRYPOINT [ "/bifrost/min_read_check/launcher.py"]
CMD [ "/bifrost/min_read_check/launcher.py", "--help"]