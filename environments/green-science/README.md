# Molecular scientific comparison environment

This separately locked Python 3.12 environment reads the actual Green GF2
artifacts and supplies PySCF 2.14.0 exact integrals/FCI and sparse-ir 2.1.6
native-grid interpolation. It does not replace or upgrade the native Green
container's provider environment.

```sh
uv sync --project environments/green-science --frozen
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  uv run --project environments/green-science --frozen python -m pytest -q \
  tests/test_green_integration.py tests/test_molecular_reference.py
```

Run the four actual native calculations using the cached pinned Green image:

```sh
docker build -f environments/green/Dockerfile -t kernelgecp-green:ef .
mkdir -p artifacts
docker run --rm -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -v "$PWD/artifacts:/artifacts" kernelgecp-green:ef \
  /opt/venv/bin/python /workspace/scripts/providers/green_molecular_runs.py \
  --output-dir /artifacts/native --source-revision "$(git rev-parse HEAD)"
```

Each invocation requires a fresh output directory. The per-run bound manifests
preserve native data/log/build/source versions, grids and density-fitting
integrals. The top receipt closes only the four actual native executions;
separate compression, native refinement, Dyson and model-reference gates remain.
The `Molecular Green GF2` workflow runs this native calculation on Linux and
retains raw artifacts. A missing local Docker daemon is not scientific evidence.

The complete failed and replacement comparison outputs are retained as
`experiments/data/molecular-green-6050ddd.tar.gz` and
`experiments/data/molecular-green-56c3211.tar.gz`. Their source revisions,
SHA256 values, original native artifact and acceptance boundaries are in
RESEARCH.md. Members have safe relative paths; saved request paths describe
the original host and are not portable inputs. The first campaign failed two
doubled-cutoff controls, and its consumed final queries were retired before
the second generation. Neither archive is a continuous-domain certificate.
