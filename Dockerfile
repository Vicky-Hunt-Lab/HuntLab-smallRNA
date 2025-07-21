FROM mambaorg/micromamba:latest

COPY --chown=$MAMBA_USER:$MAMBA_USER docker_env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python -m pip install /tmp

RUN rm -rf /tmp/*