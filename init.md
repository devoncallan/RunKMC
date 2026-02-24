## Working with This Repo and User

### Communication cadence

- Stay in sync: deliver findings first, wait for a confirmation (“please implement”) before touching code.
- Separate analysis from action: run an analysis/plan phase, only implement after explicit approval.
- Summarize assumptions: spell out what you think the fix is before editing; this keeps alignment and gives the user a checkpoint.
- Show your work: when reviewing files, run the `sed`/`rg` commands so the command log proves you read them—no silent scans.
- Match tone and detail: concise but thorough explanations with references to actual code and outputs are preferred.

### Code changes

- Respect existing conventions: reuse existing patterns (structs, helpers) instead of dropping in “from scratch” code unless approved.
- Extract small helpers before large edits: incremental, reversible changes are easier to review.
- Words like “fix” or “implement” should trigger a pause and confirmation. The user decides when implementation begins.

### Profiling/reporting

- When asked to “analyze,” assume read-only mode and return findings only. Offer a proposed fix, but wait for confirmation before editing.

### Summary for new agents

If you stick to: diagnose first, show evidence (command output), wait for explicit approval, and align with the repo’s style, you’ll stay in sync with the user and avoid rework.
EOF
'
## Working with This Repo and User

### Communication cadence

- Stay in sync: deliver findings first, wait for a confirmation (“please implement”) before touching code.
- Separate analysis from action: run an analysis/plan phase, only implement after explicit approval.
- Summarize assumptions: spell out what you think the fix is before editing; this keeps alignment and gives the user a checkpoint.
- Show your work: when reviewing files, run the `sed`/`rg` commands so the command log proves you read them—no silent scans.
- Match tone and detail: concise but thorough explanations with references to actual code and outputs are preferred.

### Code changes

- Respect existing conventions: reuse existing patterns (structs, helpers) instead of dropping in “from scratch” code unless approved.
- Extract small helpers before large edits: incremental, reversible changes are easier to review.
- Words like “fix” or “implement” should trigger a pause and confirmation. The user decides when implementation begins.

### Profiling/reporting

- When asked to “analyze,” assume read-only mode and return findings only. Offer a proposed fix, but wait for confirmation before editing.

### Summary for new agents

If you stick to: diagnose first, show evidence (command output), wait for explicit approval, and align with the repo’s style, you’ll stay in sync with the user and avoid rework.
