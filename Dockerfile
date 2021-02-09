# This is intended to run in Local Development (dev) and Github Actions (test/prod)
# BUILD_ENV options (dev, test, prod) dev for local testing and test for github actions testing on prod ready code
ARG BUILD_ENV="prod"
ARG MAINTAINER="kimn@ssi.dk;"
ARG BIFROST_COMPONENT_NAME="bifrost_assemblatron"

#---------------------------------------------------------------------------------------------------
# Programs for all environments
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.8.2 as build_base
ONBUILD ARG BIFROST_COMPONENT_NAME
ONBUILD ARG BUILD_ENV
ONBUILD ARG MAINTAINER
LABEL \
    BIFROST_COMPONENT_NAME=${BIFROST_COMPONENT_NAME} \
    description="Docker environment for ${BIFROST_COMPONENT_NAME}" \
    environment="${BUILD_ENV}" \
    maintainer="${MAINTAINER}"
RUN \
    conda install -yq -c conda-forge -c bioconda -c default snakemake-minimal==5.31.1; \
    conda install -yq -c conda-forge -c bioconda -c default bbmap==38.58; \
    conda install -yq -c conda-forge -c bioconda -c default skesa==2.4.0; \
    conda install -yq -c conda-forge -c bioconda -c default minimap2==2.17; \
    conda install -yq -c conda-forge -c bioconda -c default samtools==1.11; \
    conda install -yq -c conda-forge -c bioconda -c default cyvcf2==0.30.1; \
    # Don't use conda for Quast they cap the python version which causes issues with install
    pip install -q quast==5.0.2;


#---------------------------------------------------------------------------------------------------
# Base for dev environement
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.8.2 as build_dev
ONBUILD ARG BIFROST_COMPONENT_NAME
ONBUILD COPY --from=build_base / /
ONBUILD COPY /components/${BIFROST_COMPONENT_NAME} /bifrost/components/${BIFROST_COMPONENT_NAME}
ONBUILD COPY /lib/bifrostlib /bifrost/lib/bifrostlib
ONBUILD WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}/
ONBUILD RUN \
    pip install -r requirements.txt; \
    pip install --no-cache -e file:///bifrost/lib/bifrostlib; \
    pip install --no-cache -e file:///bifrost/components/${BIFROST_COMPONENT_NAME}/

#---------------------------------------------------------------------------------------------------
# Base for production environment
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.8.2 as build_prod
ONBUILD ARG BIFROST_COMPONENT_NAME
ONBUILD COPY --from=build_base / /
ONBUILD WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}
ONBUILD COPY ./ ./
ONBUILD RUN \
    pip install file:///bifrost/components/${BIFROST_COMPONENT_NAME}/

#---------------------------------------------------------------------------------------------------
# Base for test environment (prod with tests)
#---------------------------------------------------------------------------------------------------
FROM continuumio/miniconda3:4.8.2 as build_test
ONBUILD ARG BIFROST_COMPONENT_NAME
ONBUILD COPY --from=build_base / /
ONBUILD WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}
ONBUILD COPY ./ ./
ONBUILD RUN \
    pip install -r requirements.txt \
    pip install file:///bifrost/components/${BIFROST_COMPONENT_NAME}/


#---------------------------------------------------------------------------------------------------
# Additional resources
#---------------------------------------------------------------------------------------------------
FROM build_${BUILD_ENV}
ONBUILD ARG BIFROST_COMPONENT_NAME
# NA


#- Set up entry point:start ------------------------------------------------------------------------
WORKDIR /bifrost/components/${BIFROST_COMPONENT_NAME}
ENTRYPOINT ["python3", "-m", "bifrost_assemblatron"]
CMD ["python3", "-m", "bifrost_assemblatron", "--help"]
#- Set up entry point:end --------------------------------------------------------------------------
