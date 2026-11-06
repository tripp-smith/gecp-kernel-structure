# Troubleshooting production integrations

Start with the installed package version, imported file path, provider version
and the source/configuration hashes recorded by the calculation. Reproduce
with that source and its exact environment. Historical absolute artifact paths
identify the original machine; relocate extracted inputs explicitly.

| Symptom | Diagnosis and next action |
| --- | --- |
| `import kernelgecp` succeeds but an adapter raises `ImportError` | Install its isolated optional runtime. Native scientific providers are not default pip dependencies. Ordinary imports must remain available without them. |
| TRIQS mesh/statistics or beta validation fails | Use the tested TRIQS 4 interface and actual mesh coordinates. Preserve ordered blocks and distinct 0+ / beta- endpoints. Do not repair data by relabeling its mesh. |
| Frequency fit is rank deficient or poorly conditioned | Inspect retained SVD rank, cutoff and conditioning. Validate independent frequencies/time points. Positive-only samples do not imply conjugate completion or accurate unsampled tails. |
| A molecular checkpoint exists but import rejects convergence | Supply the actual provider receipt and bound manifest, complete final/source iteration pair, native log and input files. Exit zero or a checkpoint alone is insufficient. |
| EPW reconstructed imaginary-axis values look good but native restart fails | Check the supported model-plus-lossless archive route. The retained pure-fit route failed real-axis Padé comparison despite small imaginary-axis errors. Do not relax restart budgets or substitute a copied raw file as a fit-only success. |
| Native EPW order, units or chemical potential disagree | Preserve state ordinal and native frequency order. Supply the declared eV convention and an explicit offset when chemical potential changes. Reject missing provenance rather than infer it. |
| MaxEnt reports a singular solve or exhausts iterations | Retain the exception cause, complete trial history, input covariance and explicit solver controls. The native minimizer can reduce its initial damping. A different damping is a new declared configuration requiring calibration, not an automatic retry. |
| A continuation passes training but fails held-out RMS | Keep it failed. Compression error, noise estimation and model/continuation error are distinct. Interpolation does not create independent observations; propagate the complete covariance and forward operator. |
| Cholesky callback is rejected | Check length, real finite values, deterministic selected columns, diagonal agreement and the assumed symmetric PSD kernel. Callback checks do not certify unobserved kernel entries. |
| A tiny-noise shifted solve fails | Inspect diagnostics and the independently recomputed residual. Do not form an inverse or change noise to hide the failure. Recalibrate a declared model if a different noise is scientifically justified. |
| A campaign times out or has partial files | Observe the existing process/job first. Preserve logs and checkpoint state. Missing observations, orphan files and budget exhaustion remain failed/inconclusive; never interpret an empty result as zero error. |

Run the root `./scripts/verify.sh` for repository changes. Dedicated provider
jobs must execute their promised runtime; a skip cannot close a domain gate.
For canonical runs, a green wrapper is insufficient when scientific stages
were skipped. Check every required independent audit and final case inventory.

For an issue report include a small reproducible input, version/platform,
command, numerical configuration, source/environment hashes, failure status
and relevant diagnostics. Keep credentials and unrelated local files out of
reports. Existing frozen numerical archives are pickle-free; do not deserialize
an unsupported provider object to bypass validation.
