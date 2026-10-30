# Canonical provider target: Linux x86_64, Ubuntu 24.04 at an immutable image.
FROM ubuntu:24.04@sha256:a2d43bffc46670d7108bb3c0468c67c5557c4c92e83e1e34cbff1322ad1b94dc

ENV DEBIAN_FRONTEND=noninteractive \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    QE_ROOT=/opt/qe

RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates git make gcc g++ gfortran pkg-config \
      libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libfftw3-dev \
      python3 python3-venv \
    && python3 -m venv /opt/venv \
    && dpkg-query -W > /opt/epw-system-packages.txt \
    && rm -rf /var/lib/apt/lists/*

COPY environments/epw/fetch_sources.py environments/epw/sources.lock.json /workspace/environments/epw/

RUN /opt/venv/bin/python /workspace/environments/epw/fetch_sources.py \
      --root /opt/qe \
      --lock /workspace/environments/epw/sources.lock.json

WORKDIR /opt/qe
RUN ./configure MPIF90=mpif90 CC=mpicc F90=gfortran F77=gfortran \
      BLAS_LIBS=-lblas LAPACK_LIBS=-llapack \
      --enable-parallel --disable-openmp --with-scalapack=no \
      --with-fox=no --with-libxc=no > /opt/epw-configure.log 2>&1 \
      || { tail -n 100 /opt/epw-configure.log; exit 1; }
RUN (make -j2 epw > /opt/epw-build.log 2>&1 \
      || { tail -n 100 /opt/epw-build.log; exit 1; }) \
    && test -x /opt/qe/bin/pw.x \
    && test -x /opt/qe/bin/ph.x \
    && test -x /opt/qe/EPW/bin/epw.x

COPY scripts/providers/epw_smoke.py /workspace/scripts/providers/epw_smoke.py
COPY scripts/providers/epw_anisotropic_fixture.py /workspace/scripts/providers/epw_anisotropic_fixture.py
WORKDIR /work
CMD ["/opt/venv/bin/python", "/workspace/scripts/providers/epw_smoke.py", "--help"]
