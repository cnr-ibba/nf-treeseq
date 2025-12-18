#
# VERSION 0.2.1-hotfix
# DOCKER-VERSION  27.3.1
# AUTHOR:         Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
# DESCRIPTION:    A custom image with threads_arg installed
# TO_BUILD:       docker build --rm -t bunop/threads_arg .
# WITH_WAVE:      wave --containerfile Dockerfile --freeze --build-repo docker.io/bunop/threads_arg --await
# TO_RUN:         docker run --rm -ti bunop/threads_arg bash
# TO_TAG:         docker tag bunop/threads_arg:latest bunop/threads_arg:0.2.1
#

# Build stage
FROM python:3.12 AS builder

# Install build dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        build-essential \
        libboost-all-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Create virtual environment
RUN python -m venv /opt/venv

# Activate virtual environment and install packages
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir git+https://github.com/bunop/threads.git@v0.2.1-hotfix

# Final stage
FROM python:3.12-slim

LABEL maintainer="Paolo Cozzi <paolo.cozzi@ibba.cnr.it>"

ENV \
    PYTHONUNBUFFERED=1 \
    PYTHONFAULTHANDLER=1 \
    PATH="/opt/venv/bin:$PATH"

# Install only runtime dependencies
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv
