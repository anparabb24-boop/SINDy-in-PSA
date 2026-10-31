# IEEE 9-Bus MBAM + SINDy Work

This folder contains exploratory IEEE 9-bus data, generator-model, MBAM, and SINDy work. Start with [PROJECT_STATUS.md](PROJECT_STATUS.md) for the workspace inventory and current limitations.

## Prepare the user-designated healthy segment

The user-designated healthy interval is 5–10 s in the current 10-second exports. The earlier A–G event near 2.0–2.3 s is excluded. The script separates the fast waveform grid from slower generator and P/Q channels, preserves their native sampling, extracts one-cycle phasors and sequence components, and creates a temporal split:

- Fit segment: 5.0 <= t < 8.0 s
- Same-run temporal validation segment: 8.0 <= t <= 10.0 s

Run from the workspace root:

```bash
python sindy_in_PSA/prepare_healthy_data.py
```

Outputs are written to `prepared_data/healthy_5s_onward/`; the dataset decisions and limitations are in [DATASET_REPORT.md](DATASET_REPORT.md). The overview figure is `prepared_data/healthy_5s_onward/healthy_5s_onward_overview.png`.

Prepared signals are currently in raw export units. Per-unit scaling is deferred until machine and bus bases, phase convention, and current sign are verified. The temporal holdout comes from the same simulation run and does not replace an independent healthy validation simulation.

## Physical generator baseline

After preparing the data, run the two-axis/subtransient generator fit and same-run temporal validation:

```bash
python sindy_in_PSA/mbam_package/validate_generator_models.py
```

Results are written to `mbam_package/validation_results/`. The current provisional run meets the stated held-out NRMSE criterion for all outputs, but G1 and G3 parameter optimizers hit their evaluation limits. Review [VALIDATION_REPORT.md](mbam_package/validation_results/VALIDATION_REPORT.md) before starting physical MBAM.

## Network DAE status

The extracted branch topology is in [IEEE9_connections_loadflow.csv](IEEE9_connections_loadflow.csv). The network algebraic implementation is in [network_dae.py](network_dae.py); verified branch R/X/charging and transformer tap data are still needed. See [NETWORK_DAE_REPORT.md](NETWORK_DAE_REPORT.md).

## Current research boundary

`mbam_package/` contains a single-generator parameterized model and geodesic/FIM routines. It does not implement the full network-coupled 9-bus DAE. Read the package scripts and [DATASET_REPORT.md](DATASET_REPORT.md) before fitting or interpreting model parameters. Several older scripts still use absolute paths and guessed fault windows.
