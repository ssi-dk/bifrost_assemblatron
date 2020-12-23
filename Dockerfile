# This is intended to run in Local Development (dev) and Github Actions (test/prod)
# BUILD_ENV options (dev, test, prod) dev for local testing and test for github actions testing on prod ready code
ARG BUILD_ENV="prod"
ARG MAINTAINER="kimn@ssi.dk;"
ARG BIFROST_COMPONENT_NAME="bifrost_assemblatron"
ARG FORCE_DOWNLOAD=true


#---------------------------------------------------------------------------------------------------
# Programs for all environments
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.8.2 as build_base
ARG BIFROST_COMPONENT_NAME
ARG FORCE_DOWNLOAD
LABEL \
    BIFROST_COMPONENT_NAME=${BIFROST_COMPONENT_NAME} \
    description="Docker environment for ${BIFROST_COMPONENT_NAME}" \
    environment="${BUILD_ENV}" \
    maintainer="${MAINTAINER}"
RUN \
    conda install -yq -c conda-forge -c bioconda -c default snakemake-minimal==5.7.1; \
    conda install -yq -c conda-forge -c bioconda -c defaults bbmap==38.58; \
    conda install -yq -c conda-forge -c bioconda -c defaults skesa==2.3.0; \
    conda install -yq -c conda-forge -c bioconda -c defaults minimap2==2.17; \
    conda install -yq -c conda-forge -c bioconda -c defaults samtools==1.9; \
    conda install -yq -c conda-forge -c bioconda -c defaults cyvcf2==0.11.4; \
    # Note prokka has a 1 year deadline due to tbl2asn. 1.14.6 was made available Feb 20th
    conda install -yq -c conda-forge -c bioconda -c defaults prokka==1.14.6; \
    # Don't use conda for Quast they cap the python version which causes issues with install
    pip install -q quast==5.0.2;


#---------------------------------------------------------------------------------------------------
# Base for dev environement
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.8.2 as build_dev
ARG BIFROST_COMPONENT_NAME
COPY --from=build_base / /
COPY /components/${BIFROST_COMPONENT_NAME} /bifrost/components/${BIFROST_COMPONENT_NAME}
COPY /lib/bifrostlib /bifrost/lib/bifrostlib
WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}/
RUN \
    pip install -r requirements.txt; \
    pip install --no-cache -e file:///bifrost/lib/bifrostlib; \
    pip install --no-cache -e file:///bifrost/components/${BIFROST_COMPONENT_NAME}/

#---------------------------------------------------------------------------------------------------
# Base for production environment
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.8.2 as build_prod
ARG BIFROST_COMPONENT_NAME
COPY --from=build_base / /
WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}
COPY ./ ./
RUN \
    pip install file:///bifrost/components/${BIFROST_COMPONENT_NAME}/

#---------------------------------------------------------------------------------------------------
# Base for test environment (prod with tests)
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.8.2 as build_test
ARG BIFROST_COMPONENT_NAME
COPY --from=build_base / /
WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}
COPY ./ ./
RUN \
    pip install -r requirements.txt \
    pip install file:///bifrost/components/${BIFROST_COMPONENT_NAME}/


#---------------------------------------------------------------------------------------------------
# Additional resources
#---------------------------------------------------------------------------------------------------
FROM build_${BUILD_ENV}
# NA


#- Set up entry point:start ------------------------------------------------------------------------
WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}
ENTRYPOINT ["python3", "-m", "bifrost_assemblatron"]
CMD ["python3", "-m", "bifrost_assemblatron", "--help"]
#- Set up entry point:end --------------------------------------------------------------------------
