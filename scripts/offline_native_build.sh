#!/usr/bin/env bash
set -euo pipefail

: "${PYO3_PYTHON:?Set PYO3_PYTHON to the ordinary CPython interpreter for this build}"
if [[ -f SOURCE_MANIFEST.json && -z "${KERNELGECP_SOURCE_REVISION:-}" ]]; then
  KERNELGECP_SOURCE_REVISION="$("$PYO3_PYTHON" -c 'import json; print(json.load(open("SOURCE_MANIFEST.json"))["source_revision"])')"
  export KERNELGECP_SOURCE_REVISION
fi
# Rebuild provided or locally modified LGPL C sources instead of reusing binaries.
export GMP_MPFR_SYS_CACHE="${GMP_MPFR_SYS_CACHE:-_}"
cargo build --manifest-path native/Cargo.toml --release --locked --offline \
  --features extension-module
