# Project Instructions for AI Contributors

1. Keep dependencies to a minimum. Lightweight is the right weight: prefer base R and existing project dependencies when they are sufficient, and add a dependency only when it materially simplifies or strengthens the work.
2. Every commit containing AI-authored work must identify each contributing AI as a co-author using a valid Git trailer. For Codex, add:

   `Co-authored-by: OpenAI Codex <codex@openai.com>`

3. Never commit generated simulation data, fitted model weights, caches, or rendered intermediate files. Commit scripts, configuration, tests, documentation, and deliberately small fixtures only.
4. Keep the prototype runnable in phases: simulation first, ML training second, and prediction/reporting third.
5. Use native R `torch` for the ML surrogate. Do not introduce Python or `reticulate` unless the project owner explicitly requests it.
6. Before committing, run the relevant tests and confirm the commit message contains the required co-author trailer.

