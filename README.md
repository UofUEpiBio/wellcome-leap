# AMPLIFY Computational Prototype

This repository contains a toy computational prototype for the AMPLIFY Wellcome Leap proposal. It demonstrates how a mechanistic agent-based disease-transmission simulation can generate training data for a missing-data-aware machine-learning surrogate.

The prototype has three parts:

1. simulate two antibiotic-resistant organism scenarios with `epiworldR::ModelSEIRCONN()`;
2. learn individual realized secondary infections, \(R_i\), from an agent feature \(X\) and the organism scenario using native R `torch`; and
3. expose a prediction wrapper that works when both inputs are available or when only \(X\) or scenario is available.

The transmission probability is modeled as

\[
P(\text{transmission}_{is})=
\operatorname{logit}^{-1}(\alpha+\beta X_i+\delta S_s),
\]

where \(X_i\sim U(0,1)\) and \(S_s\) distinguishes the two organism scenarios. Scenario parameters will be calibrated toward early-epidemic reproduction numbers near 1.5 and 2.0.

## Project status

The scientific and implementation plan is documented in [PROTOTYPE_PLAN.md](PROTOTYPE_PLAN.md). Implementation proceeds in three independently testable phases: simulation and calibration, masked R `torch` training, and prediction/reporting.

Generated simulation data and trained model artifacts are intentionally excluded from version control.

## Simulation quick start

From the repository root:

```sh
Rscript scripts/01_smoke_test.R
Rscript scripts/02_calibrate.R 200 4
Rscript scripts/03_run_production.R 10000 4 100
```

The calibration command accepts the number of replicates per candidate and worker count. The production command accepts replicates per scenario, workers, and batch size. Outputs are written under the ignored `data/derived/` directory.

## Repository documents

- [FullProposal.md](FullProposal.md): grant narrative and scientific context.
- [PROTOTYPE_PLAN.md](PROTOTYPE_PLAN.md): modeling decisions, calibration design, ML design, and implementation gates.
- [AGENTS.md](AGENTS.md): durable instructions for AI contributors.
