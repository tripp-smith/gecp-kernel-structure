# Build from the repository root. The provider is isolated from kernelgecp's pip environment.
FROM ubuntu:24.04@sha256:68ee22004a9109b44973313b02ff7b1de250f1de0f5ecb628d73f579aaff33ac

ENV DEBIAN_FRONTEND=noninteractive \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMBA_NUM_THREADS=1 \
    CMAKE_BUILD_PARALLEL_LEVEL=2 \
    PATH=/opt/venv/bin:/opt/green/bin:$PATH

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential ca-certificates cmake git ninja-build pkg-config \
    libeigen3-dev libfftw3-dev libgmp-dev libhdf5-dev libmpfr-dev \
    libopenblas-dev libopenmpi-dev openmpi-bin python3-dev python3-venv \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m venv /opt/venv

COPY environments/green/requirements.lock /opt/green-lock/requirements.lock
RUN /opt/venv/bin/python -m pip install --require-hashes --only-binary=:all: \
    -r /opt/green-lock/requirements.lock \
    && /opt/venv/bin/python -m pip check

COPY environments/green/sources.json environments/green/build_green.py /opt/green-lock/
RUN /opt/venv/bin/python /opt/green-lock/build_green.py

WORKDIR /workspace
COPY scripts/providers/green_smoke.py /workspace/scripts/providers/green_smoke.py

CMD ["/opt/venv/bin/python", "/workspace/scripts/providers/green_smoke.py", "--help"]
