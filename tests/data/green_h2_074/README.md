# Attributed native Green fixture

These data come from the actual pinned Green 1.0.0 H2/STO-3G GF2 provider run
identified by `provider-receipt.json`: 0.74 Angstrom, beta=20 inverse Hartree,
symmetric Lowdin storage, and 83 native iterations. The original input, native
log, run settings, build provenance and provider receipt are retained unchanged.

`sim.h5` is a compact extract containing `/iter`, `/iter82`, `/iter83` and the
original grids-version attribute. Every copied numerical dataset was checked
for exact equality against the original full checkpoint file. The bound
manifest records both the original full-file hash and the extract's distinct
hash. No approximation, symmetry repair or serialized Cell deserialization was
used to create the extract.

The reader's bound v1 manifest enriches the explicit units/settings already in
the run's external provenance. `bind_green_manifest.py` performs this binding
for new complete actual runs; it refuses files absent from the source provider
receipt. The original provider readiness run does not close EF's second bond,
native-grid refinement, same-solution compression or independent model gates.
