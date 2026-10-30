# EPW interoperability evidence

The current supported archive route meets the original scientific gates; see
[supported route acceptance](#supported-route-acceptance). The preceding 17 K
fit-only sections below preserve the failed engineering path and its historical
pending statements, rather than extending the original release objective.

## Historical 17 K frequency representation

A frozen finite-grid frequency-cross experiment and its independent calibration
and final audits pass locally at source revision
`80ea4bd774842d02ad5c3205af8eebb0577287de`. The native QE 7.5 / EPW 6 data
were generated at `13915c807e3e82cb071af6bd44dd9c3a6bbb640c` in
[Linux CI run 40646021445](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40646021445).
The dedicated [native restart run 40854443322](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40854443322)
failed its real-axis comparison. Both native reads and the negative control
completed as expected. The fitted-restart obligation and delivery remain open.
The root verification
at `80ea4bd` stopped at formatting before Python tests; it is not a passing root
verification. These final rows were first consumed at `80ea4bd`; later formatting
or documentation does not create fresh blind data.

The coarse MgB2 native fixture uses `wscut=1 eV`, with 140 states and 109
positive Matsubara frequencies at 17 K; the 20 K and 30 K files are also
retained. The 17 K native file SHA-256 is
`14479ca5f9f5760d9b0b5553c794d8443dd16122159134dcae512c3f69a9defd`.
The model uses explicit real Eliashberg Nambu parity and a 1,201-pole dictionary
with an engineering cutoff of 13 eV. This dictionary cutoff is not the EPW
physical cutoff. The result is finite-data representation evidence, not a
material prediction, a GF2 or EPW solver acceleration result, or a universal
speedup. No analytic or rounding certificate is claimed.

## Frozen calibration and held-out checks

Calibration uses frequency indices 36, 40, 44; indices 45, 49, 53 were
reserved for final validation. Rank 48 failed the predeclared calibration
margin; rank 60 passed and was frozen. The rank-48 fit and candidate arrays
are retained alongside rank 60. Construction uses 256-bit precision; fitting
and evaluation use binary64 with `rcond=1e-12`.

The independent audit gives these maxima on the final held-out frequencies:

| Quantity | Measured maximum | Limit |
| --- | ---: | ---: |
| Scaled Green error | 1.5534051217667377e-9 | 1e-6 |
| Z error | 3.5248493013284588e-9 | 1e-5 |
| Gap error, eV | 2.059880239470291e-11 | 1e-6 |
| Shift chi error, eV | 1.256443636632465e-8 | 1e-6 |

The additional global-scaled Green discrepancy is 5.715563358805093e-10 and
imaginary defect is 8.326672684688674e-17. All three fresh-process scientific
arrays are identical; the independently audited final cardinal condition
maximum is 5.54836 and recovery condition maximum is approximately 1. The
independent audit, calibration and case receipts retain all metrics, controls,
conditioning and failure details. The actual EPW `imag_read`/Padé restart
comparison, including raw and fitted data plus wrong-frequency negative
control, must pass in Linux CI before that interoperability claim is made;
it is not a self-consistency iteration.

## Timings and retained evidence

The three fresh local cases report total times 1.042, 1.074 and 1.059 seconds
and peak RSS 104,071,168, 102,809,600 and 99,483,648 bytes respectively.
Their fit JSON plus NPZ are 218,763, 218,762 and 218,761 bytes. For comparison,
the original ASCII table is 1,541,361 bytes, its gzip encoding 340,551 bytes,
and the raw binary arrays 368,232 bytes. These are storage and timing
observations for this configuration, not end-to-end solver speedups. The
fit archive contains coefficients, selected design and metadata; the evidence
bundle separately retains the full original native data and fitted export,
so the fit-only size is not the storage cost of the entire validation package.

The [evidence archive](experiments/data/epw-frequency-80ea4bd.tar.gz) has
SHA-256 `a111f8f8b192a2b9a5e0822008ad86a9662cf01552ec878e65fccc074812144b`.
It contains all calibration candidates and raw arrays, freeze, all three case
archives and native exports, both original independent audit receipts, and the
native receipt, input, log and three temperature tables. The native solver's
full 140 MB working directory is not included. The
[summary](data/epw/summary.json), [calibration audit](data/epw/calibration-audit.json),
[campaign audit](data/epw/campaign-audit.json), [freeze](data/epw/freeze.json)
and [native receipt](data/epw/native-receipt.json) are available separately.
Historical absolute paths inside the original receipts preserve provenance;
summary artifact links are repository-relative. A future audit must supply
extracted paths and a matching source/runtime instead of interpreting historical
`/tmp` paths as portable. Prior unsuccessful time-cardinal and deflation
research routes remain recorded in RESEARCH.md.

The [failed restart comparison](data/epw/lossy-restart-failure.json) retains
all unchanged channel budgets and counts after repairing only padded-header
parsing. The [witnesses](data/epw/lossy-restart-witnesses.json) retain each
channel's worst budget ratio, and the
[failed-run archive](experiments/data/epw-lossy-restart-80ea4bd.tar.gz) retains
native inputs, logs, the original receipt and comparison. Its SHA-256 is
`4a8ee83dffd3e1728acb702416303a63e1769937654cb3f693ccf1513a82f1d6`.
Full real-axis tables remain in artifact 12888469873, digest
`07df1620aeda4aed67643dc559e8c268666cb0014d7d7fce77984df15fecbacf`,
and in the downloaded local artifact; durable release-asset preservation is
still required. Default-order Padé output is highly sensitive in this case.
A bounded native order-calibration study is specified in RESEARCH.md;
20 K is reserved for fresh validation after a successful freeze. No tolerance
relaxation or lossless fallback closes the failed fitted-restart obligation.

<a id="development-checkpoint-native-continuation-and-archival"></a>

## Historical development checkpoint: native continuation and archival

The following paragraphs describe the state at that checkpoint; later results
are recorded in [supported route acceptance](#supported-route-acceptance).

The fixed seven-order study at `942dd85` (run 40859909319) is terminal
inconclusive: all seven candidates fail, while the native baseline and compiled
helper rechecks pass. The complete [order-study receipt](data/epw/pade-order-study.json)
records the results and provenance. Its SHA-256 is
`e4053cd975c80f5dd7750d1ded46dafdae6316d78571d9c36bb6b91ef6d7e51f`;
artifact 12891476500 has digest
`15cd18b08c87d3776df8e7a42e6aa4c065ef6b2064be9f950ce7020541541ad8`.

A separate development diagnostic scales the observed fitted-minus-raw
perturbation by 1, 0.1, 0.01, 0.001 and 0 at fixed native order 98. It is
truth-assisted feasibility work, not an available compression method. The
17 K formatted endpoints were checked against fitted and raw inputs respectively;
eight fake-provider tests pass, but no native sensitivity result exists yet.

The explicit `EPWRestartArchive` combines a model with lossless XOR corrections
to reconstruct original field bits and validates the baseline hash. It does not
silently fall back or refit. The full 17 K bundle is 581,819 bytes versus
340,551 bytes for gzip of the original table: this is no gzip storage advantage
and does not demonstrate a pure fitted-restart pass. A candidate 20 K driver,
auditor and comparator support are under development; no 20/30 K fields have
been fitted or scored. Earlier synthetic checks passed 58 tests; cross-platform,
full-root and native gates for the new code remain open.

The full EPW science workflow runs on pushes to main or manually on a frozen
branch; a separately scoped sensitivity workflow runs on relevant pushes. This
avoids repeating the known 17 K failure on every development push. The required
native fitted-restart gate remains mandatory.

## Historical frozen 20 K archive checkpoint and coordinate audit correction

All three producer cases at `a9c2c39` pass the predeclared pure-model held-out
budgets and restore the original native file byte for byte. The first independent
audit fails on cardinal conditioning because it reconstructs coordinates using
`pi*kB*T`, while the declared binary64 API uses `beta=1/(kB*T)` followed by
`pi/beta`. At mode zero the former is three ULP away from a selected node and
has condition approximately 1.43e17; the exact archived node has condition 1.
The archived full mesh has maximum condition 10.1618 under the independent
MP256 product check. Nearby queries remain unsafe and must be rejected.

The verifier correction requires bit-identical declared training coordinates,
retains the independent product formulation and unchanged 1e4 cap, and has a
regression for both exact-node acceptance and adjacent-coordinate rejection.
It changes no fitted model, controls, data, or tolerance. Reauditing the consumed
20 K outputs is a verifier-correction replay, not fresh blind validation. The
[first failed audit](data/epw/restart20-a9c2c39-first-audit.json) and
[producer/archive checkpoint](data/epw/restart20-a9c2c39-checkpoint.json)
preserve this history; the native source corpus is also retained in the original
frequency evidence archive described above.

The complete 20 K model-plus-XOR payload is 524,604 bytes versus 1,300,981
bytes ASCII and 291,358 bytes gzip. Corrected native values encode original
information and are explicitly not fit predictions. The attributed frozen
fixture in `tests/data/epw/restart20-a9c2c39` checks cross-platform decoding and
native byte equality; compatibility requires those tests to execute on each
promised platform. A dedicated native archive workflow now requires independent
20 K auditing, unmodified EPW raw/restored comparison and the wrong-frequency
negative control. Those native gates, the corrected audit, final verification,
exact-head CI and delivery remain open. The failed 17 K fit-only restart claim
is not replaced by a successful archive decode.

## Supported route acceptance

Original EG scientific acceptance is met at source `1b5d60a`: the shared
Cauchy-frequency GECP model passes the three predeclared 20 K pure-model
holdouts and independent audit, and an explicit model-plus-lossless-XOR archive
restores the native input. In [Linux run 40979924989](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40979924989),
unmodified EPW reads raw and restored inputs and all six output channels agree
with zero error over 630,000 rows (4,500 frequencies by 140 states); the
wrong-frequency control exits 1. The native result SHA-256 is
`a1abe291a8067890b2dd9b2f958614201ca4a2d09b1e37c984d2c2cd618709ce` and
the audit SHA-256 is
`0d6c08a369bf22708d127f2f344ca9babb15ff3431bced6797d4851480b8f03b`.
The [native receipt](data/epw/native-archive-linux.json),
[independent audit](data/epw/native-archive-linux-audit.json),
[acceptance record](data/epw/supported-route-acceptance.json), and
[curated evidence](experiments/data/epw-restart20-linux-1b5.tar.gz) bind the
results. The complete workflow artifact 12944763024 is 73,163,061 bytes,
SHA-256 `5b5449aa21f40f40feb2e6a00672b6078c51075824943b828b16dea4cfb2811f`;
it is retained locally pending EJ asset publication.

This is reproduction of the earlier 20 K run, not a fresh blind validation.
The first failed independent audit and its coordinate-only verifier correction
remain documented above. The Mac 20 K payload is 524,604 bytes versus 1,300,981
ASCII and 291,358 gzip bytes; the Mac 17 K payload is 581,819 versus 1,541,361
ASCII and 340,551 gzip bytes. Linux 20 K bundles are 524,562–524,563 bytes versus 1,300,981 ASCII and
292,867 gzip bytes (model-only 218,121–218,122 bytes).
These data do not show a gzip storage advantage or universal performance.

The 17 K pure-fit native comparison failed. In the truth-assisted sensitivity
run [40978168520](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40978168520),
all four nonzero residual scales (1, 0.1, 0.01, 0.001) fail, while the exact
zero-copy control passes. No pure-fit native acceptance follows. The earlier
plan to use generic `FermionicBasis`/`FermionicFit` and the later pure-fit-only
native gate are superseded failed engineering approaches, not the original EG
objective; neither is supported or silently counted as passed. Accuracy, noise
and native budgets are unchanged. Full root and exact-head CI passed for
`1b5d60a`; the documentation revision still requires its own root and CI,
followed by merge and main verification before delivery.

The historical fit-only workflow is manual and named as a development study.
The supported native archive workflow remains automatic on main; moving the
failed experiment out of automatic production checks does not alter its
recorded failure or the required native archive comparison.
