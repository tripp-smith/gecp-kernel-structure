# Practical implementation and release objective

This is the acceptance ledger for the user's application-development goal. It
does not replace the completed mathematical `SPEC.md` or weaken its theorems.
An unchecked requirement remains open even when a related phase or test passes.
DU's revision-specific receipt is [PR #127](https://github.com/tripp-smith/gecp-kernel-structure/pull/127);
DV's full root/platform/install and merged-main closure is in
[PR #128](https://github.com/tripp-smith/gecp-kernel-structure/pull/128).
All controlled requirements of the completed DU–EA application scope are delivered in
[kernelgecp 1.1.0 on PyPI](https://pypi.org/project/kernelgecp/1.1.0/) and the
[GitHub release](https://github.com/tripp-smith/gecp-kernel-structure/releases/tag/v1.1.0).
The [release receipt](RELEASE_RECEIPT.md) records the immutable source, tested
file hashes and actual publication checks. Final repository handoff acceptance
is recorded in [PR #132](https://github.com/tripp-smith/gecp-kernel-structure/pull/132).
The maintainer-invited [executable notebook PR #10](https://github.com/SpM-lab/sparse-ir-doc/pull/10)
is submitted and locally verified; [PR #133](https://github.com/tripp-smith/gecp-kernel-structure/pull/133)
owns this repository's final adoption-record acceptance. Current upstream
review, CI approval and merge decisions remain with its maintainers.
The separately authorized [production integration extension](#production-integrations-12)
is in progress; its required 1.2.0 capabilities are not delivered by 1.1.0.

## Required end state

| Requirement | Current evidence | State |
| --- | --- | --- |
| Compare with pydlr and sparse-ir using shared observations and external accuracy criteria | DV frozen A retained failure; B12/12 passing, independent audit, merged PR #128 | delivered bounded comparison |
| Separate construction, repeated fit/evaluation, time, process memory, rank and conditioning | DV two-worker matched cohort B and explicit scopes | implemented and checked; single-machine evidence |
| Independently check both time and Matsubara frequency errors | DV95-digit all-forward and solver frequency/direct-integral checks | delivered numerical evidence; not a continuous certificate |
| Include a small physical solver loop | DV Bethe loop passes with all three providers | implemented and checked; noninteracting scope |
| Public reusable basis and cached factors | DU FermionicBasis/LeastSquaresPlan, merged PR #127 | delivered |
| Batching and complex/matrix-valued coefficients | DU arbitrary-axis public tests and tutorial, merged PR #127 | delivered |
| Explicit physical conventions and accuracy/precision contracts | DU API/guide; DW merged finite-query arithmetic/reference separation | delivered; continuous inexact rate remains a research extension |
| Safe versioned serialization and short runnable tutorials | DU bounded NPZ/JSON, corruption tests and example, PR #127 | delivered |
| Preserve existing research API defaults | unchanged through DU | ongoing invariant |
| Same-workload interoperability demo with an existing DLR/IR implementation | DV three-provider executable example, merged PR #128 | delivered |
| Concrete upstream proposal and actual adoption outreach | proposal, maintainer invitation and [actual executable notebook PR #10](https://github.com/SpM-lab/sparse-ir-doc/pull/10); source/schema/execution/full-book checks pass under upstream lock | contribution submitted; upstream review/CI approval/merge externally controlled |
| Profile and prototype the measured Rust bottleneck | profiling, actual pilot, DX optional backend and exact-main platform acceptance in PR #130 | delivered |
| Keep Rust only with demonstrated end-to-end gain and numerical parity | DX frozen six-process cohort: 50.15% application, 47.76% process-wall reduction; paired differences <4.26e-15, physical <2.57e-11 | delivered for streaming workload; small reuse/cache-hit limits explicit |
| Cross-platform built-wheel and outside-checkout install CI | DY exact-main producer and six public OS/Python install lanes plus public sdist installation; release receipt | delivered for the recorded platform matrix |
| Distributed typing and extension stubs if applicable | DY installed positive/negative typing consumers pass for pure/native distributions and public installs | delivered |
| Reproducible production benchmarks and citation-ready methods document | DV/DW/DX frozen evidence; published METHODS and schema-validated CITATION with explicit claim boundaries | delivered |
| Publish a tested versioned release and installable distributions | public v1.1.0 release from exact-main source 9a679930db326bda2456bfa91b3256f37359e4ce; three wheels and original sdist | delivered |
| PyPI publication for the practical Python module | four public files match the tested manifest; workflow 40244470229 attempt 2 passes all public installs | delivered |
| Finite-precision/inexact-pivot bridge including coefficient and solve errors | DW portable repeated six-case artifact, actual finite-query fit and finite-grid true-pivot budgets, two conditional Lean theorems, merged PR #129 | delivered bounded bridge; continuous inexact rate remains a research extension |
| Implementation validated, committed, PR-merged and pushed | DT through DZ delivered, PR #131 publication and PR #132 exact-main closure; invited upstream PR submitted | release/handoff delivered; final adoption-record acceptance in PR #133 |

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
complete. Actual public installation now closes publication. Native wheel
targets are Linux x86_64 (manylinux 2.28) and macOS arm64 (deployment tag 11.0).
Their public installation and execution pass on hosted Ubuntu/macOS runners
with Python 3.12/3.13; Windows uses the universal Python fallback. No native
Windows, Linux arm64 or macOS x86_64 wheel is claimed. Broader adoption still
requires maintainer action.

## Publication and external follow-up

GitHub repository access is available. On 2026-10-24 EDT the user explicitly
confirmed adding the PyPI pending Trusted Publisher for project `kernelgecp`,
repository `tripp-smith/gecp-kernel-structure`, workflow `publish.yml`, environment
`pypi`. The user subsequently explicitly authorized publication on 2026-10-25.
That saved publisher uploaded version 1.1.0 on 2026-10-25 at
23:23:48–23:23:53 UTC. The original tested bytes are public, including complete
corresponding native source, notices and dependency inventory published on
GitHub before upload. No user-supplied upload token or password was requested
or stored in the repository.

The [publication receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/131#issuecomment-6404917491)
preserves the initial five postchecks that exhausted their short PyPI JSON
visibility window with HTTP 404. Only those failed verification jobs were
rerun; upload was not repeated. Attempt 2 passes all checks, and independent
public JSON/simple-index hashes match all four tested distributions. Final
handoff and a bounded readiness-backoff correction are tracked in PR #132;
released assets remain unchanged. The public upstream proposal is delivered,
while maintainer acceptance remains pending.

## Production integrations 1.2

The user authorized implementation of the complete [production integration
contract](PRODUCTION_INTEGRATIONS.md). This is additional application scope;
the original mathematical goals and published 1.1.0 acceptance above remain
frozen. Production support means tested public adapters and real-provider
workflows, with a closed-loop DMFT integration. Green and EPW internal solvers
remain native. Upstream acceptance is not a delivery gate.

| Phase | Mandatory output | Current acceptance |
| --- | --- | --- |
| EB | Locked optional runtimes; actual provider imports and minimal calculations | all six platform/provider jobs pass at `d5f58a1`; replay locks/notices and final exact-source/main acceptance in [PR #134](https://github.com/tripp-smith/gecp-kernel-structure/pull/134) |
| EC | Frequency fitting, archive v2/v1 reader, actual TRIQS adapters | complete, merged/pushed main `64e9f34`; full root and all four exact-main workflows pass ([receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/135#issuecomment-6436634601)) |
| ED | Deterministic interacting Bethe IPT with independent references | complete, merged/pushed main `7530ec7`; full root and all four exact-main workflows pass ([receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/136#issuecomment-6445705017)) |
| EE | Actual CT-HYB DMFT and independent-replica uncertainty | `7200ec0` complete fresh Linux campaign and all seven independent scientific gates pass; final source-bound replay, root/exact-head CI and delivery remain open in PR #141 |
| EF | Native molecular GF2 runs and supported Green importer | complete, merged/pushed main `07c4115`; root and all exact-main workflows pass ([receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/137#issuecomment-6481863146)) |
| EG | EPW MgB2 native output/compression/restart interoperability | delivered at main `a268f24`, root and all five exact-main workflows pass; explicit model-plus-lossless native route passes, pure-fit route remains failed ([receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/138#issuecomment-6518307400)) |
| EH | Actual covariance-aware MaxEnt including EE observations | normalized Mac synthetic/calibration pass; fresh measured GECP/flat fails RMS 2.0232; Linux provider failure and predeclared 256-chain candidate remain open in PR #142 |
| EI | Column-oracle covariance factorization, SciPy/CPU GPyTorch, real data | complete, merged/pushed main `f6ba348`; root and all five exact-main workflows pass ([receipt](https://github.com/tripp-smith/gecp-kernel-structure/pull/139#issuecomment-6493607721)) |
| EJ | Verified 1.2.0 public distributions and final handoff | pending all domain gates |

EB/EC/ED/EF/EG/EI are delivered, approximately 70% of the new integration program by completed acceptance obligations; EE/EH/EJ remain open. Each phase requires its named
acceptance evidence, root verification, exact-revision CI, authorized merge,
and pushed-main verification. A runtime lane that skips because its promised
provider is missing does not close EB or its dependent application. Failed
mandatory cases remain open; no universal speedup is required.

Remaining delivery order follows the current dependencies. Finish
EE exact-head/root/merge/main delivery after its accepted Linux campaign and
independent audit. EG is delivered with failed routes retained. Complete EH's
frozen measured-data comparison, independent audit and real-provider CI before
its delivery. EJ then audits every original domain obligation, updates guides,
compatibility/migration/notices and handoff metadata, tests installed release
distributions, publishes 1.2.0, verifies public installations and synchronizes
exact-main checks. Green and EPW internal solver acceleration stays outside scope.

The active goal retains that full end state. Bounded Terra workers own individual
implementations, independent audits or evidence packages; the primary agent owns
canonical launches, shared numerical decisions, acceptance and merges. Workers
return artifact paths, hashes, commands, outcomes and open obligations. Failed
worker calls do not count as reviews, and completed subsets do not redefine the
goal. Long scientific runs use durable checkpoints and confirmed live handles.

EG's supported route, reproduction and preserved failed coordinate audit are
recorded in [EPW_INTEROPERABILITY.md](EPW_INTEROPERABILITY.md#supported-route-acceptance).
The generic basis and pure-fit native approaches failed and remain unsupported;
their failed evidence is preserved separately from the delivered native archive.
