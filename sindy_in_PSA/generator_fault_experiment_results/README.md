# Generator-level robust STLSQ fault-detection pilot

This is a single-event, generator-terminal, binary detection pilot. It is not an implementation of MBAM or constrained SINDy and does not modify the IEEE 9-bus DAE or physical parameters.

## Method

Each generator has its own healthy-only model. The source CSV combines waveform and slow streams; the existing `split_merged_streams` function is reused. Waveforms remain at 1,440 Hz and slow channels remain at native 60 Hz. No slow-channel interpolation is performed. Required channels and missing values are audited; same-time fast/slow row collisions in the merged source are deliberately resolved by the shared splitter, then each separated stream is checked for repeated timestamps and its sampling interval. Phasor endpoint windows and one-cycle margins are excluded.

Features are `[omega/w0-1, theta_rel, Pm-Pm_fit_median, Efd-Efd_fit_median]`; targets are `d(theta_rad)/dt` and `d(omega/w0)/dt`. Theta is interpreted as degrees from its observed range, converted to radians and unwrapped. Derivatives use a centered 7-point, quadratic Savitzky-Golay filter on native slow samples (100 ms span), which smooths noise but blurs rapid fault edges. Each target has an independent robust STLSQ model: standardized columns, Huber IRLS (delta 1.5), ridge 1e-4, sequential threshold 0.08, 15 support updates maximum.

The dimensionless score is `sqrt(mean_j((y_dot_j - yhat_dot_j)^2 / s_j^2))`, where `s_j` is the robust MAD scale of the corresponding healthy 5–8 s derivative target. For each generator, threshold is the empirical higher-method 99th percentile of healthy calibration scores at 8–10 s, after excluding the 50 ms centered-derivative support at the fit/calibration boundary and the final 1/60 s endpoint. Alarm requires two consecutive 60 Hz scores above threshold, approximately one cycle. No fault samples enter fit, feature scaling, or threshold selection. Fault evaluation excludes 50 ms at each measured event boundary to exceed the requested one-cycle guard and ensure the centered derivative stencil stays inside the event. Derivatives are centered and use three future slow samples, so the score is for offline diagnosis and is not a causal real-time detector as implemented.

The comparison baseline is `|V2|/|V1|` from one-cycle phasors of the generator-terminal waveforms. Its threshold is independently calibrated by the same 99th-percentile rule on that generator's healthy 8–10 s waveform. This baseline is not a SINDy result.

## Input audit

| Generator | Healthy raw/wave/slow rows | Fault raw/wave/slow rows | Slow dt (healthy/fault) | Fault slow missing channels |
|---|---:|---:|---:|---|
| G1 | 15000/14401/601 | 15000/14401/601 | 0.01666667/0.01666667 s | healthy {'Va': 0, 'Vb': 0, 'Vc': 0, 'Ia': 0, 'Ib': 0, 'Ic': 0, 'theta': 14399, 'w': 14399, 'Efd': 14399, 'Pmpu': 14399}; fault {'Va': 0, 'Vb': 0, 'Vc': 0, 'Ia': 0, 'Ib': 0, 'Ic': 0, 'Pm_pu': 14399, 'theta': 14399, 'w': 14399, 'Excitation V': 14399} |
| G2 | 15000/14401/601 | 15000/14401/601 | 0.01666667/0.01666667 s | healthy {'Va': 0, 'Vb': 0, 'Vc': 0, 'Ia': 0, 'Ib': 0, 'Ic': 0, 'theta': 14399, 'w': 14399, 'Efd': 14399, 'Pmpu': 14399}; fault {'Va': 0, 'Vb': 0, 'Vc': 0, 'Ia': 0, 'Ib': 0, 'Ic': 0, 'theta': 14399, 'w': 14399, 'excitation V': 14399, 'Pmpu': 14399} |
| G3 | 15000/14401/601 | 15000/14401/601 | 0.01666667/0.01666667 s | healthy {'Va': 0, 'Vb': 0, 'Vc': 0, 'Ia': 0, 'Ib': 0, 'Ic': 0, 'theta': 14399, 'w': 14399, 'Efd': 14399, 'Pmpu': 14399}; fault {'Va': 0, 'Vb': 0, 'Vc': 0, 'Ia': 0, 'Ib': 0, 'Ic': 0, 'theta': 14399, 'w': 14399, 'excitation V': 14399, 'Pmpu': 14399} |

## Event boundary evidence

The measured event interval is estimated from negative-to-positive sequence voltage ratio exceedance against the independent healthy 8–10 s q99 threshold; each one-cycle phasor is centered and has one-cycle support. Event evaluation uses a 50 ms guard at each boundary, exceeding one cycle and covering the centered derivative stencil.

| Generator | Threshold | Measured start (s) | Measured end (s) | Exceedance phasors |
|---|---:|---:|---:|---:|
| G1 | 0.00117298 | 2.008333 | 2.291667 | 18 |
| G2 | 0.00117304 | 2.008333 | 2.291667 | 18 |
| G3 | 0.00117301 | 2.008333 | 2.291667 | 18 |

## Results

| Generator | q99 threshold | Healthy FA rate | Healthy score median / p95 | Fault score median / p95 | Fault alarm fraction | Detected | Delay (s) | V2/V1 fault median | V2/V1 alarm fraction |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| G1 | 1.4616 | 0.862% | 0.8388 / 1.389 | 960.2 / 1483 | 100.000% | True | 0.075 | 0.8113 | 100.000% |
| G2 | 1.1891 | 0.862% | 0.3138 / 1.095 | 4949 / 6452 | 100.000% | True | 0.075 | 0.8822 | 100.000% |
| G3 | 7.5138 | 0.862% | 3.456 / 7.142 | 2819 / 4450 | 100.000% | True | 0.075 | 0.8439 | 100.000% |

## Active standardized coefficients

The table lists the fitted coefficients in standardized feature/target coordinates; the intercept is included. Inactive library terms are zero. `model_coefficients.csv` also stores feature means/scales and target normalization needed to reproduce predictions.

| Generator | Target | Active terms (standardized coefficient) |
|---|---|---|
| G1 | omega_pu_dot_s | intercept=-0.168433, omega_dev=1.6138, theta_rel_rad=0.857939, Pm_dev=1.88989, Efd_dev=-1.01827 |
| G1 | theta_dot_rad_s | intercept=0.784781, omega_dev=1.5286 |
| G2 | omega_pu_dot_s | intercept=-0.196084, omega_dev=1.93276, theta_rel_rad=-1.21127, Pm_dev=1.03897, Efd_dev=-1.8479 |
| G2 | theta_dot_rad_s | intercept=0.625836, omega_dev=1.65407 |
| G3 | omega_pu_dot_s | intercept=-0.190282, omega_dev=0.84833, theta_rel_rad=5.79746, Pm_dev=1.95084, Efd_dev=-4.97121 |
| G3 | theta_dot_rad_s | intercept=0.816442, omega_dev=1.9102, theta_rel_rad=0.232171, Efd_dev=-0.258472 |

Per-sample predictions/residuals are in `G*_healthy_predictions.csv` and `G*_fault_predictions.csv`. Separate voltage baseline series are `G*_voltage_ratio_baseline.csv`; plots are `G*_residual_detection.png`. `generator_metrics.csv`, `event_summary.csv`, and `input_audit.json` provide machine-readable summaries. Run from the workspace root with `python sindy_in_PSA/generator_fault_experiment.py`.

## Integrity and limitations

The fault rows are used only for forward prediction and event scoring. There is one healthy run and one fault run per generator; the three generator records are correlated observations from one system event, not independent fault cases. Healthy calibration is a temporal segment of the same run as training. The event interval is operationally inferred from fault-run negative-sequence ratio crossing its healthy-calibrated threshold; therefore the voltage-ratio baseline confirms this same electrical event signature and is not independent ground truth. Boundary timing has one-cycle phasor uncertainty. Savitzky-Golay smoothing blurs onset and low-rate slow channels limit temporal resolution. A successful alarm therefore supports only a single-event offline pilot claim; it does not establish general real-time fault detection across fault locations, operating points, or independent simulations. The voltage baseline is terminal-only and not bus localization.

No physical parameters or DAE code are changed. No MBAM or constrained SINDy is implemented.
