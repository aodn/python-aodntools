FROM ubuntu:20.04

ARG BUILDER_UID=9999
ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Australia/Hobart
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV PATH /home/builder/.local/bin:$PATH
ENV PYTHON_VERSION 3.8.13

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    rm -rf /var/lib/apt/lists/*

RUN add-apt-repository ppa:rael-gc/rvm && apt-get update

RUN apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
    libmagic1 \
    libudunits2-dev \
    python3-dev \
    wget \
    libffi-dev \
    # Pyenv pre-requisites (from https://github.com/pyenv/pyenv/wiki#suggested-build-environment)
    make build-essential libssl-dev zlib1g-dev \
    libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm \
    libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Set-up necessary Env vars for PyEnv
ENV PYENV_ROOT $HOME/.pyenv
ENV PATH $PYENV_ROOT/shims:$PYENV_ROOT/bin:$PATH

# Install pyenv
RUN set -ex \
    && curl https://pyenv.run | bash \
    && pyenv install $PYTHON_VERSION \
    && pyenv global $PYTHON_VERSION \
    && pyenv rehash \
    && chmod -R a+w $PYENV_ROOT/shims

RUN pip install --upgrade pip==22.1.2 setuptools==63.1.0 wheel

RUN pip install \
    Cython==0.29.30 \
    bump2version==1.0.1 \
    numpy==1.23.0

RUN useradd --create-home --no-log-init --shell /bin/bash --uid $BUILDER_UID builder
USER builder
WORKDIR /home/builder
