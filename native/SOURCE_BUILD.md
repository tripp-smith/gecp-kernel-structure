# Corresponding native sources

Generate the versioned source asset from a frozen Git commit and complete locked
Cargo vendor tree:

```sh
cargo vendor --manifest-path native/Cargo.toml --locked /absolute/path/vendor
python scripts/prepare_corresponding_source.py build --project-root . \
  --revision FULL_GIT_SHA --vendor /absolute/path/vendor --output /path/assets \
  --report /path/assets/source-report.json
```

The asset includes the checked Python/Rust/build sources, exact Cargo lock,
complete vendored crates, dependency licenses, a deterministic inventory and
per-file hashes. Its report also identifies the generated license bundle and
SBOM used by native-wheel staging. Source preparation validates the Git binding,
crate checksums, complete payload and resource limits; a directory of cached
binaries is not corresponding source.

Unpack the asset, install Rust 1.92.0 and an ordinary CPython interpreter, set
`PYO3_PYTHON`, and run `scripts/offline_native_build.sh` from its root. Cargo uses
the included vendor configuration and `--offline`. The script disables the
GMP/MPFR binary cache by default so modified C sources are actually rebuilt.
GMP/MPFR native selftests remain enabled. Native library prerequisites are a C
compiler, make, m4 and diffutils on Linux, or Xcode command-line tools on macOS.

Wheel packaging additionally needs the pinned Maturin build tool. The ordinary
PEP 517 backend remains Hatchling and can build the Python fallback without Rust.
For experiments that modify vendored LGPL sources, preserve notices and update
the affected `.cargo-checksum.json` file entries before rebuilding. Modified
sources are distinct from the original release manifest and must not be described
as the original checked artifact.

The published native wheel must link to its actual versioned source asset and
carry distribution-accurate license metadata and notices. Source inspection and
an unpacked offline build are release gates; merely generating an archive does
not establish that a wheel corresponds to it.
