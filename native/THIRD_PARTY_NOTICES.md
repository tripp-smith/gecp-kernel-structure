# Native dependencies

The Python package and Rust wrapper use Apache-2.0; see the repository LICENSE.
Native wheels also incorporate dependencies with additional licenses. Their
metadata must describe the bundled distribution rather than repeat Apache-2.0
alone.

The locked Rug 1.30.0 and gmp-mpfr-sys 1.7.1 crates use LGPL-3.0-or-later.
The latter supplies GMP 6.3.0 and MPFR 4.2.2. Full GNU license texts are in
`native/licenses/`; the complete corresponding-source asset also retains every
vendored crate's own license and source files. It includes upstream MPC sources
as part of the complete gmp-mpfr-sys crate, although the float-only backend does
not enable MPC.

Source preparation emits a complete third-party license bundle, a dependency
inventory and a deterministic combined native license expression from all locked
declared licenses. Native-wheel staging includes those generated materials and
binds them to the same source report. Other MIT/Apache/Unicode and platform
dependency notices must not be dropped merely because the GNU texts are present.

See SOURCE_BUILD.md for source preparation and rebuild instructions. The release
receipt supplies the actual source asset URL and hashes; no unpublished or
development artifact is represented as a completed release.
