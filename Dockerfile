FROM ubuntu:latest

ARG BUILDER_UID=9999
ARG DEBIAN_FRONTEND=noninteractive

ENV TZ=Australia/Hobart
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV PATH /home/builder/.local/bin:$PATH
ENV PYTHON_VERSION 3.5.2

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    rm -rf /var/lib/apt/lists/*

RUN add-apt-repository ppa:rael-gc/rvm && apt-get update

RUN if [ X"$PYTHON_VERSION" = X"3.5.2" ]; \
        then apt-get install -y libssl1.0-dev; \
        else apt-get install -y  libssl-dev; \
    fi

RUN apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
    libmagic1 \
    libudunits2-dev \
    python3-dev \
    wget \
    libffi-dev \
    # Pyenv pre-requisites
    make zlib1g-dev libbz2-dev libreadline-dev \
    libsqlite3-dev wget curl llvm libncurses5-dev \
    libncursesw5-dev xz-utils tk-dev libffi-dev \
    liblzma-dev python-openssl \
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

RUN pip install --upgrade pip==20.3.4 setuptools==50.3.2

RUN pip install \
    Cython==0.29 \
    numpy>=1.13.0 \
    bump2version==0.5.10 \
    wheel

RUN useradd --create-home --no-log-init --shell /bin/bash --uid $BUILDER_UID builder
USER builder
WORKDIR /home/builder
