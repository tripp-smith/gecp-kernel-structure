# Optional-provider source and notice inventory

kernelgecp's base Python source uses the repository [LICENSE](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/LICENSE).
The optional native evaluator's corresponding source, dependency SBOM and
license bundle are separate release assets described in
[native/SOURCE_BUILD.md](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/native/SOURCE_BUILD.md) and
[native/THIRD_PARTY_NOTICES.md](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/native/THIRD_PARTY_NOTICES.md).
Preserve those exact assets with the wheel that names them. For 1.2.0, consult the release manifest for the identities and bindings of
the corresponding-source and notice assets accompanying each wheel.

Scientific providers run in isolated environments and are not bundled into the
kernelgecp wheels. Their source observations, binary package metadata and actual
execution receipts have different scopes. This index links the existing
machine-readable records rather than inferring license terms from a project
name or applying one provider's license to an entire environment.

| Inventory | Recorded scope |
| --- | --- |
| [Provider sources](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/environments/provider-licenses.json) | Exact provider/dependency source commits, observed license/notice paths and file hashes; includes explicit incomplete observations |
| [Linux TRIQS environment](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/environments/locks/triqs-linux-64-conda-inventory.json) | Actual conda package builds, transitive constraints, package hashes and license metadata |
| [macOS arm64 TRIQS environment](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/environments/locks/triqs-osx-arm64-conda-inventory.json) | Separate actual platform-specific conda inventory |
| [DLR/IR comparison wheels](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/environments/locks/dlr-physics-inventory.json) | Observed installed comparison-provider wheel license files and exact identities |
| [Green native build](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/environments/locks/green-linux-build.json) | Pinned native source/build observations and actual build environment |
| [Green sources](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/environments/green/sources.json) | Exact recursively pinned Green source components |
| [EPW setup](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/environments/epw/README.md) | qe-7.5/EPW 6.0 build and execution boundaries |
| [GP runtime](https://github.com/tripp-smith/gecp-kernel-structure/blob/v1.2.0/environments/gp/README.md) | CPU-only joint environment and hash-locked upstream wheels |

The provider-source inventory explicitly records devXlib's absent standalone
license file and the observed source-header notice. It does not infer a
repository-wide license from that header. Source-file observations also do not
assert that every upstream binary distribution contains identical notices.

For exact transitive versions and MPI/HDF5 variants use the corresponding
platform/build lock and the actual run's exported conda/pip/system inventories.
Release validation must retain these records with source identities and
calculation settings. An environment assembled later is not retroactive
provenance for an earlier scientific result. If distributing the optional
runtimes themselves, assemble their complete corresponding source and notices
for those actual artifacts; the kernelgecp wheel inventory does not do that.
