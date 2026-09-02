# Generated data

Simulation and calibration scripts write generated RDS files to `data/derived/`.
That directory is ignored by Git. Recreate outputs with the scripts in `scripts/`.

The optional multiscale workflow writes separated truth, observation,
fitted-profile, and ABM layers under `data/derived/multiscale/`. These generated
objects are also ignored and must not be committed.
