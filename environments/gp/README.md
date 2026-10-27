# Optional CPU GP provider runtime

This environment is separate from the package's base dependencies and the
historical benchmark environment. GPyTorch requires mpmath <=1.3; its compatible
joint runtime uses mpmath 1.3.0. Linux explicitly selects the official CPU PyTorch
index rather than downloading the CUDA runtime. The initial targets are Python
3.12 on Linux x86_64 and macOS arm64 >=14 (the latter restriction belongs to
PyTorch, not the base kernelgecp package).

```sh
uv sync --project environments/gp --frozen
uv run --project environments/gp --frozen python scripts/providers/gp_smoke.py \
  --source-revision "$(git rev-parse HEAD)" --output artifacts/gp.json
```

The probe uses actual GPyTorch distributions and low-rank operators. Solves,
log determinants and Gaussian log probabilities are checked against independent
SciPy Cholesky calculations. It is a runtime smoke check, not the EI real-data
accuracy/scaling benchmark. A failed import or numerical check exits nonzero.

The committed uv lock records artifact URLs and hashes. Execution receipts
record the actually imported versions, source revision and platform. Support is
not established by lock resolution alone; the phase PR records actual execution.
The environment uses upstream wheels and is not bundled into kernelgecp wheels.
