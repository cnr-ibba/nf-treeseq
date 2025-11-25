#
# VERSION 0.2.1
# DOCKER-VERSION  27.3.1
# AUTHOR:         Paolo Cozzi <paolo.cozzi@ibba.cnr.it>
# DESCRIPTION:    A custom image with threads_arg installed
# TO_BUILD:       docker build --rm -t bunop/threads_arg .
# TO_RUN:         docker run --rm -ti bunop/threads_arg bash
# TO_TAG:         docker tag bunop/threads_arg:latest bunop/threads_arg:0.2.1
#

FROM python:3.12-slim-bookworm

LABEL maintainer="Paolo Cozzi <paolo.cozzi@ibba.cnr.it>"

ENV \
    PYTHONUNBUFFERED=1 \
    PYTHONFAULTHANDLER=1

# update packages and install threads_arg
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends procps && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir threads-arg==0.2.1
