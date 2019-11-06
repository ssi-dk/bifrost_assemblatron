FROM ssidk/bifrost-base:2.0.5

ARG version="2.0.5"
ARG last_updated="31/07/2019"
ARG name="assemblatron"
ARG full_name="bifrost-${name}"

LABEL \
    name=${name} \
    description="Docker environment for ${full_name}" \
    version=${version} \
    resource_version=${last_updated} \
    maintainer="kimn@ssi.dk;"

#- Tools to install:start---------------------------------------------------------------------------
RUN \
    conda install -yq -c conda-forge -c bioconda -c defaults bbmap==38.58; \
    conda install -yq -c conda-forge -c bioconda -c defaults skesa==2.3.0; \
    conda install -yq -c conda-forge -c bioconda -c defaults minimap2==2.17; \
    conda install -yq -c conda-forge -c bioconda -c defaults samtools==1.9; \
    conda install -yq -c conda-forge -c bioconda -c defaults cyvcf2==0.11.4; \
    conda install -yq -c conda-forge -c bioconda -c defaults prokka==1.14.0; \
    # Don't use conda for Quast they cap the python version which causes issues with install
    pip install -q quast==5.0.2; 
#- Tools to install:end ----------------------------------------------------------------------------

#- Additional resources (files/DBs): start ---------------------------------------------------------
RUN cd /bifrost_resources && \
    wget -q https://raw.githubusercontent.com/ssi-dk/bifrost/master/setup/adapters.fasta && \
    chmod +r adapters.fasta
#- Additional resources (files/DBs): end -----------------------------------------------------------

#- Source code:start -------------------------------------------------------------------------------
RUN cd /bifrost && \
    git clone --branch ${version} https://github.com/ssi-dk/${full_name}.git ${name};
#- Source code:end ---------------------------------------------------------------------------------

#- Set up entry point:start ------------------------------------------------------------------------
ENV PATH /bifrost/${name}/:$PATH
ENTRYPOINT ["launcher.py"]
CMD ["launcher.py", "--help"]
#- Set up entry point:end --------------------------------------------------------------------------