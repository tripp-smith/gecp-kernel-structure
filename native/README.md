# Optional reference time backend

The `kernelgecp._native` module accelerates time-kernel evaluation and cached
transpose-LU substitution at 128 bits. Enable it with
`FermionicBasis(..., reference_backend="rust")`. Python remains the default;
selection, Matsubara evaluation, NumPy fitting and array products keep their
existing implementations. Inspect the per-component basis diagnostics.

The first implementation supports ranks 1–128 and dimensionless frequency
magnitudes at most 16384. Unsupported inputs or native arithmetic failures are
explicit errors. Physical normalization and exports remain binary64. Numerical
parity is evidence, and the separate finite-point certification API supplies
explicit error budgets where it succeeds.

Build with Rust 1.92.0 and ordinary CPython 3.12 or later:

```sh
PYO3_PYTHON=/path/to/python cargo build --manifest-path native/Cargo.toml \
  --release --locked --features extension-module
```

`cargo test --manifest-path native/Cargo.toml --locked` links the test interpreter;
the `extension-module` feature is reserved for extension builds. Set
`KERNELGECP_SOURCE_REVISION` to the full checked Git commit when producing an
official artifact. Development builds report `unrecorded`; official wheel gates
require the recorded source, version and target.

Hatchling still builds the universal Python wheel and ordinary source installs.
The explicit Maturin build produces a platform wheel for the same project. Native
Linux/macOS claims require their actual installed-wheel checks. Windows retains
the universal Python fallback until a GNU/MSYS2 native lane is verified.

See [SOURCE_BUILD.md](SOURCE_BUILD.md) and
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) before distributing a native wheel.
