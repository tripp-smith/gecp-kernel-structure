# Published kernelgecp 1.1.0

[PyPI 1.1.0](https://pypi.org/project/kernelgecp/1.1.0/) and the
[GitHub release](https://github.com/tripp-smith/gecp-kernel-structure/releases/tag/v1.1.0)
are public. Install with:

```sh
python -m pip install kernelgecp==1.1.0
```

The immutable release tag is `v1.1.0`, source commit
`9a679930db326bda2456bfa91b3256f37359e4ce`, tree
`bfab90004a3a062b55b5903daf93ee386f992bf6`. Later repository handoff changes
do not move that tag or replace any published artifact.

## Verified artifact identities

The [public release manifest](https://github.com/tripp-smith/gecp-kernel-structure/releases/download/v1.1.0/release-manifest.json)
has SHA-256 `2d2d9b2cb7cb1c35f7762cc079aa3722ad4e335434d57cbe3c98ba4894124420`.
Both public PyPI JSON and its simple file index match all four original tested
distribution hashes:

| Distribution | SHA-256 |
| --- | --- |
| `kernelgecp-1.1.0-py3-none-any.whl` | `c2b32125f0f3d54646c4bf33058ca38262b89f9c88f29608b8f05c99aa4a056e` |
| `kernelgecp-1.1.0.tar.gz` | `9a6c13c77a00f00a88e5e7b5f6995afa383452fbfcaf946f8b517da22e5ef916` |
| `kernelgecp-1.1.0-cp312-abi3-manylinux_2_28_x86_64.whl` | `537d1273c3991955b841bbc7c75dcdc2a89278e7312ce51231abec54c8a566e8` |
| `kernelgecp-1.1.0-cp312-abi3-macosx_11_0_arm64.whl` | `461bb3fab23898f30ad8658d0fa095edc492bd525c440500971c14a4b956b0ae` |

Complete native sources were public before upload:

| Source asset | SHA-256 |
| --- | --- |
| `kernelgecp-1.1.0-native-source.tar.gz` | `4fb5fdf72652e38540ef6882d7f05bd9183082d09bd0e43c0c8bb007da211a49` |
| `kernelgecp-1.1.0-native-licenses.txt` | `27d65eb649e510281db7bc0d682e54e0057146baf840eb3c5bee7e85f47dfebb` |
| `kernelgecp-1.1.0-native-dependencies.cyclonedx.json` | `99bfbb6ed41d8a6ec194ac362174d700e0022d68e1976a95b0d23201cf1fd3ab` |

The manifest also binds the source report and versioned methods/citation/release
notes. Native wheels identify that exact source revision, tree and source hash.
Full licenses and the sanitized locked-dependency inventory accompany them.

## Acceptance and publication

[PR #131](https://github.com/tripp-smith/gecp-kernel-structure/pull/131)
records the final branch, real downloaded-artifact integration, independent
reviews and exact-main checks. The native-enabled main root passes 1,726 tests
with eight optional skips, 3,635 Lean jobs and 878 permitted axiom reports,
plus typing, formatting and dependency audit. Root log SHA-256 is
`26aa302eb21d080dd3e423c9b8bf2a418187ffd186a4fe92f5f835e6c8531546`.

Successful exact-main producers:

- [Lean/Python CI 40243697353](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40243697353).
- [Pure distribution matrix 40243697368](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40243697368).
- [Original/rebuilt portable and native producer 40243697349](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40243697349).

[Publication 40244470229, attempt 2](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40244470229/attempts/2)
is fully successful. It promotes unchanged exact-main artifacts, verifies the
public corresponding source, uploads through the saved `publish.yml` / `pypi`
Trusted Publisher with PyPA attestations, and runs six clean public-index wheel
install lanes plus a separate public sdist installation. Upload timestamps are
2026-10-25 23:23:48–23:23:53 UTC. No credential was requested from the user or
stored in source.

The supported test matrix is Python 3.12/3.13 on Linux, macOS and Windows.
Native artifacts cover Linux x86_64 with glibc 2.28+ and macOS arm64 11+;
Windows and unbuilt native architectures use the universal Python path.
Installed native lanes exercise both Python and Rust, physical predictions,
actual precision certificates, portable/native archive restoration and typing.

## Retained metadata-visibility failure

In [attempt 1](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40244470229/attempts/1),
upload, public sdist and Windows 3.13 verification passed. Five other postchecks
received HTTP 404 from the JSON metadata endpoint and exhausted the original
25-second retry-delay budget while the public simple index already served the
files. First-failure log SHA-256:
`40ae836ba7622b63901f3f40554038188e2eb30d91db824086f7683866e0739f`.

After metadata became visible, exact public hashes verified and only the five
failed verification jobs were rerun. Upload/source jobs were not repeated;
all six wheel lanes and the public sdist lane pass in the combined successful
second attempt. The final repository increases future capped retry delays to
at most 135 seconds, preserving six strict attempts and all hash/version checks.
Bounded absence and delayed-success regressions use no real sleeps. The released
tag and its original workflow remain immutable; no initially all-green run is
claimed.

## Scientific scope and final repository closure

[METHODS.md](METHODS.md) preserves the exact-run versus sampled/numerical
boundaries. The Rust benefit is specific to the frozen streaming workload;
short reuse/cache hits have different costs. No continuous inexact-pivot rate,
universal fastest method or broader native support follows from publication.
[Upstream outreach](https://github.com/SpM-lab/sparse-ir-doc/issues/6#issuecomment-6400741382)
is delivered; maintainer acceptance remains pending.

[Final handoff PR #132](https://github.com/tripp-smith/gecp-kernel-structure/pull/132)
records its exact checked branch/main revisions and remote
agreement externally, so this document does not claim to contain its own final
commit. Protected unrelated user notes remain outside the committed work.
