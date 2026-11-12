# GECP project token usage and API-equivalent cost

**Through the kernelgecp 1.2.0 handoff, the project used 3,569,732,986 recorded tokens, with an estimated API-equivalent model cost of $1,106.** The estimate uses OpenAI Standard text rates checked on November 12, 2026. Usage with a recorded token breakdown accounts for 69.0% of the total and costs about $825 at those rates. The remaining 31.0% has thread-level totals and an estimated cost of about $281. The account invoice and any usage outside the recoverable records are unknown.

The project record spans August 15 to November 12, 2026. The cutoff is `2026-11-12T14:46:47.911Z`, when the cost-report request began, after the public release and final handoff. The accounting includes project work before that request and excludes this report. [Final handoff and publication evidence](https://github.com/tripp-smith/gecp-kernel-structure/pull/144#issuecomment-6569479579).

## Project approach

We developed the formal and numerical results together. The original scope paired Lean proofs and axiom audits with reproducible experiments on GECP and the fermionic kernel. The application work turned those results into a Python library, compared it against DLR and sparse IR, and added Rust for the measured streaming bottleneck. The production release added native workflows across correlated materials, quantum chemistry, superconductivity, spectroscopy, and covariance approximation. We kept the phase contracts and acceptance evidence in the repository and used targeted subagents for implementation and independent review. [Research handoff](FINAL_HANDOFF.md) · [Application acceptance ledger](PRODUCTIZATION.md).

The accounting captures model calls for planning, implementation, experiment troubleshooting, verification, reviews, and release coordination where those calls appear in the local records. An experiment running in GitHub Actions can consume compute without adding model tokens. Its logs can become model input when an agent reads them. The dollar figures cover model tokens and exclude workstation costs, GitHub Actions compute and storage, separately priced API tools, and account-specific subscription or credit charges.

## Recorded token consumption

| Work | Tokens with a breakdown | Thread totals only | Total recorded |
| --- | ---: | ---: | ---: |
| Primary threads | 2,243,706,302 | 209,496,982 | 2,453,203,284 |
| Subagents | 157,985,316 | 820,926,550 | 978,911,866 |
| Automatic approval review | 61,165,609 | 76,452,227 | 137,617,836 |
| **Total** | **2,462,857,227** | **1,106,875,759** | **3,569,732,986** |

The detailed portion contains **2,455,146,788 input tokens** and **7,710,439 output tokens**. Of the input, **2,408,258,816 tokens were cached reads** and **46,887,972 were uncached**. That is a **98.09% cache hit rate by input tokens**. The records show zero cache-write tokens in this portion. The **3,755,065 reasoning tokens** are already included in output. Repeated processing of shared instructions, conversation history, and tool results can generate many input tokens over a long project. Cached reads still have a price. [OpenAI usage fields](https://developers.openai.com/api/docs/guides/agents-api/observability) · [Prompt-caching accounting](https://developers.openai.com/api/docs/guides/prompt-caching).

## API-equivalent estimate

| Runtime model label | Recorded tokens | Standard USD per million, input / cached / output | Cost from detailed usage | Estimated from thread totals | Total estimate |
| --- | ---: | ---: | ---: | ---: | ---: |
| GPT-5.6 Sol | 1,139,133,660 | $4.00 / $0.40 / $20.00 | $475.17 | $106.03 | $581.20 |
| GPT-6.1 Sol | 1,878,596,992 | $2.00 / $0.10 / $10.00 | $327.10 | $141.12 | $468.21 |
| GPT-6.1 Terra | 414,384,498 | $1.00 / $0.05 / $5.00 | $13.82 | $22.42 | $36.24 |
| Automatic approval review | 137,617,836 | GPT-6.1 Terra rates as a proxy | $8.94 | $11.17 | $20.10 |
| **Total** | **3,569,732,986** | | **$825.02** | **$280.74** | **$1,105.76** |

Rates are for Standard short-context text. Requests above 272,000 input tokens use twice the input and cached-read rates and 1.5 times the output rate. The detailed estimate applies that rule per recorded request. GPT-5.6 Sol uses the promotional rates listed on the pricing date. The recorded labels identify the model configuration available to this audit. Exact serving snapshots, billed tiers, and account charges are unavailable. Each total is rounded after calculation. [GPT-5.6 Sol](https://developers.openai.com/api/docs/models/gpt-5.6-sol) · [GPT-6.1 Sol](https://developers.openai.com/api/docs/models/gpt-6.1-sol) · [GPT-6.1 Terra](https://developers.openai.com/api/docs/models/gpt-6.1-terra).

The local `codex-auto-review` label has no corresponding public rate in the pricing sources used here. Automatic review makes separate model calls. The estimate applies Terra rates as an explicit proxy. Pricing that row at Sol rates would put the full estimate at **$1,125.87**. [OpenAI automatic-review accounting](https://developers.openai.com/api/docs/guides/agents-api/sessions/approvals).

### Sensitivity to missing token categories

The central estimate prices each thread-only token at the observed average cost per total token for its indexed model label. The observed model mix includes the recorded shares of input, output, cached reads, and long requests. Early threads and missing subagent sessions may differ from that mix. Holding the detailed estimate fixed gives these scenarios for the 1.107 billion thread-only tokens.

| Assumption for thread-only records | Thread-only estimate | Project estimate |
| --- | ---: | ---: |
| Observed cache and context mix by model | $280.74 | **$1,105.76** |
| 99% of input cached | $254.91 | $1,079.93 |
| 95% of input cached | $364.05 | $1,189.07 |
| 90% of input cached | $500.47 | $1,325.49 |

The fixed-cache scenarios retain the observed output and context-tier proportions by model and change only the cache hit rate of thread-only input. These scenarios show sensitivity to assumptions. Actual spending could fall outside them if the missing usage had a different context mix, model, cache policy, or service tier. In the measured portion alone, applying the published long-context rule added **$107.05** over short-context prices.

## How the records were reconciled

1. Select project threads from the local Codex thread index (`state_5.sqlite`, `threads`) using the repository checkout or worktree path and the Git origin. Check parent and attachment links for additional indexed threads. Include primary work, delegated sessions, and automatic approval review through the cutoff. Count each thread ID once. The index contains 206 in-scope threads, including 65 with a retained rollout and 141 with only an indexed aggregate.
2. For retained rollouts, difference successive `total_token_usage` snapshots in `token_count` events and ignore repeated snapshots. Match each increment to its turn model and check it against `last_token_usage`. The audit found 14,139 positive increments, all matching the last-usage record. For threads with no retained rollout, use the index's `tokens_used` value once. The retained project rollouts start on October 7. August usage is represented by the earlier thread's index total, without its original category breakdown. Threads with zero recorded usage contribute zero to the recorded sum, with their actual usage unverified.
3. Subtract cached reads and cache writes from input to get ordinary input. Multiply ordinary input, cached reads, and output by their applicable prices per million. Reasoning is part of output. Recorded cache writes were zero. For thread-only usage, divide the observed model cost by the observed total tokens and multiply by the indexed total. Use the latest indexed model label for those threads because their model changes cannot be reconstructed. [OpenAI token accounting](https://developers.openai.com/api/docs/guides/token-counting) · [OpenAI pricing](https://developers.openai.com/api/docs/pricing).
4. Keep the goal counters as historical milestones. The v1 handoff recorded **1,801,132 aggregate goal tokens** and the original research goal closed at **24,738,403 aggregate goal tokens**. The research goal updates occur inside the same project session records used for this audit. The v1 work is represented by the earlier indexed project thread. Neither goal counter is added to the session ledger. The counters have different accounting scopes and their historical estimates remain in the handoff. [Implementation-run metadata](FINAL_HANDOFF.md#implementation-run-metadata).

This reconstruction can be repeated from local metadata and token counters without copying prompts or tool output into the repository. A billing reconciliation would require per-request usage and billed rates from the account, including the missing rollout periods and the treatment of automatic review.
