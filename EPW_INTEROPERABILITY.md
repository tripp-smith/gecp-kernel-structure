# EPW frequency representation evidence

A frozen finite-grid frequency-cross experiment and its independent calibration
and final audits pass locally at source revision
`80ea4bd774842d02ad5c3205af8eebb0577287de`. The native QE 7.5 / EPW 6 data
were generated at `13915c807e3e82cb071af6bd44dd9c3a6bbb640c` in
[Linux CI run 40646021445](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40646021445).
The dedicated [native restart run 40854443322](https://github.com/tripp-smith/gecp-kernel-structure/actions/runs/40854443322)
and delivery remain open until their outcomes are checked. The root verification
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
