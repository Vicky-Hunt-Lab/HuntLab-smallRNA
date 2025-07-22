FROM mambaorg/micromamba:latest

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER docker_env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python -m pip install /tmp

RUN rm -rf /tmp/*

ENV PATH=${PATH}:/opt/conda/bin