# Practical implementation and release objective

This is the acceptance ledger for the user's application-development goal. It
does not replace the completed mathematical `SPEC.md` or weaken its theorems.
An unchecked requirement remains open even when a related phase or test passes.
DU's revision-specific root/CI/install and authorized delivery receipt is
[PR #127](https://github.com/tripp-smith/gecp-kernel-structure/pull/127).

## Required end state

| Requirement | Current evidence | State |
| --- | --- | --- |
| Compare with pydlr and sparse-ir using shared observations and external accuracy criteria | DT pilot, PR #126 | partial: pilot only |
| Separate construction, repeated fit/evaluation, time, process memory, rank and conditioning | DT measurement schema and pilot | delivered pilot; production comparison pending |
| Independently check both time and Matsubara frequency errors | DT independent time checks | frequency checks open |
| Include a small physical solver loop | none | open |
| Public reusable basis and cached factors | DU FermionicBasis/LeastSquaresPlan, merged PR #127 | delivered |
| Batching and complex/matrix-valued coefficients | DU arbitrary-axis public tests and tutorial, merged PR #127 | delivered |
| Explicit physical conventions and accuracy/precision contracts | DU API/guide and time/frequency analytic checks | implemented; continuous/rounding bridge open |
| Safe versioned serialization and short runnable tutorials | DU bounded NPZ/JSON, corruption tests and example, PR #127 | delivered |
| Preserve existing research API defaults | unchanged through DU | ongoing invariant |
| Same-workload interoperability demo with an existing DLR/IR implementation | DT comparison adapters | full demo open |
| Concrete upstream proposal and actual adoption outreach | research only | open |
| Profile and prototype the measured Rust bottleneck | DT setup/reuse timings, no Rust prototype | open |
| Keep Rust only with demonstrated end-to-end gain and numerical parity | none | open; rejection is valid only after measurement |
| Cross-platform built-wheel and outside-checkout install CI | DU main CI40143517992: six platform/version jobs pass | delivered; repeat for final release |
| Distributed typing and extension stubs if applicable | DU py.typed and installed positive/negative consumers pass | delivered; repeat for final release |
| Reproducible production benchmarks and citation-ready methods document | bounded DT pilot and research docs | open |
| Publish a tested versioned release and installable distributions | no release | open |
| PyPI publication for the practical Python module | pending Trusted Publisher saved and reloaded after user confirmation | upload and public install remain open |
| Finite-precision/inexact-pivot bridge including coefficient and solve errors | exact-run theorem and limited rounded-model pivot certificates | open |
| Implementation validated, committed, PR-merged and pushed | DT and DU delivered; DV in progress | open for full objective |

## Delivery order and gates

1. DU: public cached physical basis/fitting/evaluation, batches, serialization,
   time/frequency checks and tutorials; explicit arithmetic policy without
   changing existing GECP defaults.
2. Extend matched benchmarks and demonstrate a small solver loop and actual
   cross-library round trip. Measure matched achieved accuracy, not just equal
   library tolerance parameters; keep independent final validation untouched.
3. Profile the resulting application, prototype the dominant eligible Rust
   operation, and retain or reject it from measured end-to-end evidence.
4. Deliver the numerical bridge, wheel/install matrix, methods/citation material,
   release artifacts and publication. Seek upstream adoption with a focused
   independently useful proposal and link the actual external evidence.

Steps can share work, but a phase receipt is not proof that unrelated rows are
complete. Publication/adoption and broad platform claims require their actual
external results. GitHub releases alone do not establish a PyPI publication.

## Current external prerequisites

GitHub repository access is available. On 2026-10-24 EDT the user explicitly
confirmed adding the PyPI pending Trusted Publisher for project `kernelgecp`,
repository `tripp-smith/gecp-kernel-structure`, workflow `publish.yml`, environment
`pypi`. The signed-in PyPI management page displays it after reload. No package
has been uploaded, and this does not reserve the name or prove publication.
The matching gated workflow/environment and tested release remain to be
delivered. Never place upload tokens in chat or source files.
