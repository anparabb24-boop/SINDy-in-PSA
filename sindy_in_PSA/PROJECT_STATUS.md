# Project Status: IEEE 9-Bus MBAM + Constrained SINDy

**Assessment date:** 2026-10-31  
**Scope:** Workspace inventory refreshed after addition of `mbam_package` and 10-second simulation exports.

## Current state

The project has separate 10-second records in `healthy_datasets/G*_10s_H.csv` and earlier post-fault records in `fault datasets/G*_Fag_10s.csv`. The newly added healthy records remain balanced through the earlier 2.0–2.5 s event interval and are used as a separate no-fault candidate for trajectory transfer. Exact source SLX provenance is unconfirmed. Network parameters and wiring were extracted from the available `/Users/Arya/Desktop/IEEE9BusSystem.slx`; the supplied audit also describes a distinct 10-second `IEEE9BusSystem(1).slx` variant.

Step 3 now has source-derived network tables, an analyst-selected 100 MVA reference base, positive-sequence Ybus, and coupled generator/network DAE residual assembly. The Ybus agrees with the stored load-flow snapshot to 0.16065 MVA maximum nodal mismatch. The previous generator fit was transferred without refitting to the new no-fault run. Measured generator currents give low KCL residuals in a boundary-driven network reconstruction; model-generated currents leave larger residuals. The full closed-loop DAE is not yet validated to completion.

## Inventory

### Data and artifacts

- New 10-second generator records, located in `fault datasets/`: `G1_Fag_10s.csv`, `G2_Fag_10s.csv`, and `G3_Fag_10s.csv`. Per user instruction, use only 5–10 s as healthy.
- Separate no-fault candidate records: `healthy_datasets/G1_10s_H.csv`, `G2_10s_H.csv`, and `G3_10s_H.csv`.
- The earlier files `healthy datasets/G1_H.csv`, `G2_H.csv`, and `G3_H.csv` are absent from the current working tree; Git reports them deleted. Preserve this workspace state unless the user asks to restore them.
- Fault generator records: `fault datasets/G1_Fag.csv`, `G2_Fag.csv`, and `G3_Fag.csv`.
- New 10-second bus waveform records, located in `fault datasets/`: `IVPQ_B2_H.csv`, `IVPQ_B7_H.csv`, and `IVPQ_B8_H.csv`. Per user instruction, use only 5–10 s as healthy.
- The earlier `IVPQ_B1_H.csv`, `B4_H.csv`, and `B5_H.csv` files are absent from the current working tree; Git reports them deleted.
- Fault bus waveform records: `fault datasets/IVPQ_B1_Fag.csv`, `B4_Fag.csv`, and `B5_Fag.csv`.
- Load-flow table: `IEEE9BUS_loadflow_results.csv` (16 rows, 15 columns; 23 missing cells across fields that appear inapplicable for some block types).
- Existing figure: `ieee9_healthy_baseline_overview.png`.
- Extracted network topology and load-flow branch flows: `IEEE9_connections_loadflow.csv`. The source workbook is still deleted in the working tree; this CSV was extracted from its Git reference without restoring it.

The previous preparation script documents the post-fault 5–10 s split. New clean-run checks separately process the added `healthy_datasets` files without changing their raw inputs.

### Dataset observations

- The new 10-second generator files contain `t`, three phase voltages/currents, and four slower channels. G1 uses `Pm_pu`, `theta`, `w`, and `Excitation V`; G2/G3 use `Pmpu`, `theta`, `w`, and `excitation V`. Stream separation yields 14,401 waveform samples at 1,440 Hz and 601 slow-channel samples at 60 Hz.
- The new bus files use `t`, three phase voltages, three phase currents, and `P, Q`. Stream separation yields 14,401 waveform samples and 601 valid P/Q samples per bus. The combined 15,000 rows are not a single uniform-rate stream.
- The older 1-second records remain in the fault-data folder. They have 61 generator samples at 1/60 s and 1,501 bus samples with six duplicate timestamps; `P`/`Q` each have 1,440 missing values.
- The `fault datasets/G*_Fag_10s.csv` records show a strong A-G event near 2.0–2.3 s; their 5–10 s segments are post-fault recovery. The new `healthy_datasets/G*_10s_H.csv` records are a separate no-fault candidate, but exact source model version is unconfirmed.
- The earlier 0.3–0.5 s interval is not supported by the observed event in these new 10-second records. The visible event is near 2.0–2.3 s.
- `mbam_package/prep10.py` assumes 60 Hz, generator voltage/power bases, pole-pair values, current sign, and a phasor convention. These are explicit code assumptions, not yet verified from a Simulink model or export metadata.
- The measurements are mixed-rate (60 Hz slower channels versus approximately 1.44 kHz waveforms). Exact alignment, units, phase convention, and bases still need verification.

## Existing code assessment

- `ieee9_healthy_baseline_sindy.py` implements robust regression and a healthy/fault residual workflow, but its configured paths point to an old `/Users/Arya/Desktop/vs_code/...` location. Its `HEALTHY_BUSES` configuration also includes a fault B1 file. The configured healthy training window and known fault interval are assumptions that need data verification.
- `ieee9_robust_sindy_fim.py` implements waveform preprocessing, one-cycle phasors, symmetrical components, robust SINDy, and FIM-guided pruning. Its header correctly says this is not formal MBAM. Its paths also point to the old `vs_code` folder and its configured fault interval is 0.30–0.50 s.
- `sindy_9bus.py` and `sindy_in_psa01.py` are related multi-bus residual/FIM-pruning experiments. They use hard-coded paths and call the FIM eigen-analysis “MBAM-style”; this is not geodesic continuation to a physical model boundary.
- `smib_robust_sindy_fim.py` is a robust SINDy/FIM-pruning SMIB demo, and explicitly says it is not full MBAM. It is not an identified physical IEEE 9-bus model.
- `SINDy_inPSA.py`, `02.py`, `03.py`, `04.py`, `nongaussian.py`, and `nongaussian01.py` contain parameter-identification, synthetic-response, or table-generation experiments. Their MBAM labels/refinements need auditing before research claims; the visible routines are not evidence of physical geodesic MBAM.
- `ybus.py` currently loads the load-flow CSV and generator records, initializes a short vector from generator channels, and leaves `compute_algebraic_residuals` empty. It is a DAE scaffold, not a functioning Ybus-coupled DAE. It also uses an absolute local path.
- `mbam_package/mbam_core.py` defines a six-state two-axis/subtransient machine model with measured terminal voltage/angle, excitation, and mechanical power as inputs. It fits electrical/rotor outputs and implements finite-difference sensitivities, an FIM, geodesic continuation, and candidate limits (time constants to small/large values; damping to zero). This is physical-model MBAM machinery, distinct from earlier feature-pruning heuristics, but limiting reductions and model accuracy still need evidence from held-out data.
- `mbam_package/prep10.py` handles phasor extraction, sequence components, normalization, and resampling/alignment to slow generator channels. Its bases, phase convention, and current sign are hard-coded assumptions that need validation. `load10.py` and `real_run.py` reference obsolete `/mnt/user-data/attachments/...` and `/tmp/claude-0/mbam` paths, so the workflow is not portable/runnable in the current workspace without path fixes.
- `mbam_package/mbam_run.py` provides nominal machine parameters and an initial fitting / MBAM reduction loop. These parameters are embedded in code; provenance and correspondence to actual simulation blocks have not been verified. The package does not implement network algebraic equations or the complete 9-bus DAE.
- `mbam_package/synth_test.py` is a synthetic-data fit demonstration, not evidence of performance on independent healthy measurements.
- `prepare_healthy_data.py` separates the merged fast and slow streams, uses the user-designated 5–10 s segment, creates 5–8 s fitting and 8–10 s same-run temporal validation splits, and extracts one-cycle phase/sequence phasors. It retains waveform samples and does not interpolate slow channels.
- `mbam_package/validate_generator_models.py` fits the full physical machine model on 5–8 s and evaluates a continuous simulation on 8–10 s. All 12 held-out output NRMSE values are at or below the provisional 0.15 criterion (maximum 0.1436 for G3 active power). G2's terminal phasor power sign agrees with logged B2 P/Q to median relative complex-power error 0.00062. The optimizer converged for G2; G1 and G3 hit their function-evaluation limits, so their parameter estimates remain provisional.
- `network_dae.py` now provides a positive-sequence Ybus and constant-impedance-load algebraic residual implementation, but requires explicit branch R/X/charging/tap data. The load-flow workbook has topology, lengths, and point flows only; it does not provide those branch parameters. Generator-current coupling and dynamic network validation are not yet implemented.
- Several similarly named or scratch scripts exist. The `sindy_9bus.py`/`sindy_in_psa01.py` and `03.py`/`04.py` pairs differ rather than being byte-identical; they should be compared semantically before consolidation. `tempCodeRunnerFile.py` is an editor scratch file.
- The workspace-level `README.md` describes Lorenz/GCN experiments, not the IEEE 9-bus project.

## Missing research inputs and deliverables

The exact `IEEE9BusSystem(1).slx` version is absent; the unnumbered model and supplied parameter audits are available. Branch and transformer parameters are now extracted. The selected transformer phase convention has a strong operating-point consistency check but no isolated Simscape phasor experiment. Constrained-SINDy validation and held-out fault metrics remain open.

## Workspace state and preservation notes

At refresh time, Git reports the earlier generator and bus files and load-flow XLSX as deleted, while the load-flow CSV, new 10-second data, and `mbam_package` are untracked. `sindy_in_PSA/ybus.py` remains modified. The new 10-second files are in `fault datasets/`; source inputs were not moved or rewritten. Prepared outputs are under `prepared_data/healthy_5s_onward/`.

## Recommended next priority

Step 1 data preparation is complete for the post-fault interval and the new no-fault candidate. Step 2's prior same-run fit remains provisional; the clean-run transfer metrics are in `mbam_package/validation_results/independent_healthy_run_generator_transfer_metrics.csv`. Step 3 has the Ybus and DAE assembly but remains in progress because fitted generator currents do not close the independent network reconstruction closely enough and transmission-bus voltages are not directly logged. See `STEP3_VALIDATION_REPORT.md` for results and the completion gate.
