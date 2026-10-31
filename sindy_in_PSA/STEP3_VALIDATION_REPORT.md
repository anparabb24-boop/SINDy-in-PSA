# Step 3 — IEEE 9-Bus Network DAE

## Status

**Step 3 is not complete: the coupled positive-sequence DAE integrates and solves its algebraic equations, but it does not yet pass independent healthy-trajectory validation.** The unchanged fitted generator models have systematic terminal phase errors on the separate 5–10 s healthy candidate, and no independent transmission-bus voltage channels are present. The complete coupled run and metrics are documented below; no model parameters were retuned to reduce residuals.

The requested `IEEE9BusSystem(1).slx` was not available locally; `/Users/Arya/Desktop/IEEE9BusSystem.slx` is present. The supplied handoff reports the unnumbered model has StopTime 0.1 s and the `(1)` model has StopTime 10 s, with changed operating-point tags. It reports shared major network parameterization, which matches the directly inspected line, transformer, and load settings. Exact 10-second run provenance remains unverified. The audited Bus 2 display tag (1.070 pu) differs from its PV subsystem label (1.025 pu, 163 MW); the former is reference metadata, not the PV setpoint.

## Topology and load interpretation

The SLX wiring gives lines 4–5, 4–6, 5–7, 6–9, 7–8, and 8–9; transformers 4–1, 7–2, and 9–3. Generator connections are G1@Bus1→TF4-1→Bus4, G2@Bus2→TF7-2→Bus7, and G3@Bus3→TF9-3→Bus9. Load A is Wye at Bus5 (125 MW/50 MVAr); Load B is Wye at Bus6 (90 MW/30 MVAr); Load C is Wye at Bus8 (100 MW/35 MVAr).

The additional block is wired at Bus6 and is Wye connected. Its effective component structure is `R` (source block `component_structure_PQ=1`): a 50 MW resistive load with zero reactive consumption. Although `Qpos=30 MVAr` is stored, it is inactive for the selected R-only structure. The model load-flow snapshot independently records approximately 50.93 MW and zero MVAr at 1.009 pu voltage, consistent with a 50 MW fixed impedance. Effective total load: **365 MW + 115 MVAr**, including **140 MW + 30 MVAr at Bus6**. The attached handoff's 145 MVAr total counts the inactive 30 MVAr field; the SLX component structure and load-flow snapshot support 115 MVAr.

The Wye load blocks parameterize fixed RLC components from rated powers. At 60 Hz their positive-sequence current is linear in voltage: `Iload=conj(Srated_pu)/|Vnom_pu|² * V`. They are constant impedance in `ieee9_dae.py`, not ideal constant-PQ current sinks. [MathWorks Wye-Connected Load](https://uk.mathworks.com/help/simscape-electrical/ref/wyeconnectedload.html)

## Bases and per-unit conversion

The physical Simscape network does not declare a common MVA base. The DAE adopts **100 MVA as an analyst-selected reference base**, supported by the three transformer ratings and supplied handoff; it is not a discovered Simulink setting. Voltage bases are 16.5, 18, and 13.8 kV line-to-line at buses 1–3 and 230 kV at buses 4–9.

`Zbase=VLL_base²/Sbase`; line series impedance pu is `Zohm/Zbase`, and line charging pu is `Bsiemens*Zbase`. Transformer winding pu impedances are summed on their 100 MVA rating and scaled by `(Sbase/Srated)*(Vwinding_base/Vbus_base)²`; the selected voltage bases match, so this scale is one. Magnetizing branch pu admittance is `0.002−j0.002` at each HV winding node. Generator differential equations retain their machine bases; terminal current converts to network base as `Ipu_system=Ipu_machine*Srated/Sbase` because voltage bases match. The machine impedance base-conversion helper is in `ieee9_dae.py`.

## Line conversions

For the balanced positive-sequence network, the three-phase line model uses `R1=R−Rm`, `L1=L−M`, and `C1=Cg+3Cl`; each line is stamped as a π equivalent with half the total charging at each end. This is an explicit positive-sequence reduction, not a phase-domain fault model. [MathWorks Transmission Line (Three-Phase)](https://www.mathworks.com/help/simscape-electrical/ref/transmissionlinethreephase.html)

| Branch | R1 (Ω) | X1 (Ω) | Total B1 (S) | R1 pu | X1 pu | B1 pu |
|---|---:|---:|---:|---:|---:|---:|
| 4–5 | 5.2900 | 44.9656 | 0.00033271 | 0.01000 | 0.08500 | 0.17601 |
| 4–6 | 8.9930 | 48.6677 | 0.00029867 | 0.01700 | 0.09200 | 0.15800 |
| 5–7 | 16.9280 | 85.1698 | 0.00057846 | 0.03200 | 0.16100 | 0.30600 |
| 6–9 | 20.6310 | 89.9312 | 0.00067674 | 0.03900 | 0.17000 | 0.35799 |
| 7–8 | 4.4965 | 38.0893 | 0.00028169 | 0.00850 | 0.07200 | 0.14901 |
| 8–9 | 6.2950 | 53.3216 | 0.00039509 | 0.01190 | 0.10080 | 0.20900 |

## Transformer phase and Ybus

All transformers use Yn on winding 1 (HV) and Delta 1 o'clock on winding 2 (LV). MathWorks defines Delta 1 o'clock as a 30-degree lag relative to the connected wye voltage. With W1 wired to buses 4/7/9 and W2 to buses 1/2/3, `VLV/VHV=∠−30°`. `build_ybus.py` uses the equivalent HV-side tap angle `+30°` under its convention `VHV/a=VLV`. [MathWorks Two-Winding Transformer (Three-Phase)](https://www.mathworks.com/help/simscape-electrical/ref/twowindingtransformerthreephase.html)

The tap sign was checked against the stored operating-point snapshot using bus voltage tags, machine powers, rated-power fixed-impedance loads, and reconstructed Ybus. At +30° tap angle, maximum complex nodal mismatch is **0.16065 MVA** and RMS mismatch is **0.12584 MVA**. Reversing the tap sign gives **1859.86 MVA** maximum and **1491.94 MVA** RMS mismatch. This supports the selected terminal convention; it is not a controlled Simscape transformer experiment.

Bus order: `[1,2,3,4,5,6,7,8,9]`. Common-base Ybus (pu, rounded to 4 decimals):

```text
[ 0.0603-17.3609j, 0, 0, 8.6282+15.0651j, 0, 0, 0, 0, 0 ]
[ 0, 0.0512-15.9998j, 0, 0, 0, 0, 7.9556+13.8819j, 0, 0 ]
[ 0, 0, 0.0582-17.0646j, 0, 0, 0, 0, 0, 8.4819+14.8075j ]
[-8.7327+15.0048j, 0, 0, 3.3696-39.3106j, -1.3652+11.6039j, -1.9422+10.5107j, 0, 0, 0 ]
[ 0, 0, 0, -1.3652+11.6039j, 2.5527-17.3380j, 0, -1.1876+5.9751j, 0, 0 ]
[ 0, 0, 0, -1.9422+10.5107j, 0, 3.2242-15.8409j, 0, 0, -1.2820+5.5882j ]
[ 0, -8.0443+13.8307j, 0, 0, -1.1876+5.9751j, 0, 2.8578-35.4469j, -1.6170+13.6975j, 0 ]
[ 0, 0, 0, 0, 0, 0, -1.6170+13.6975j, 2.7722-23.3031j, -1.1551+9.7846j ]
[ 0, 0, -8.5828+14.7493j, 0, 0, -1.2820+5.5882j, 0, -1.1551+9.7846j, 2.4974-32.1559j ]
```

Full precision is in `mbam_package/validation_results/ybus_reference_check.json`. Phase shifting makes Ybus nonsymmetric. Transformer magnetizing branches and line charging are included. Zero-sequence behavior, saturation, and fault switching are outside this balanced positive-sequence model.

## Healthy-run trajectory validation

The new `healthy_datasets/G1_10s_H.csv`, `G2_10s_H.csv`, and `G3_10s_H.csv` differ from the earlier `G*_Fag_10s.csv` faulted records. New terminal waveforms remain balanced around 2.0–2.5 s; the old records show large phase imbalance. The new 5–10 s segment is used as a separate no-fault candidate. Rotor angle still drifts, so it is dynamic data rather than an equilibrium. Exact source SLX version for this run is unconfirmed.

Previously fitted generator parameters were reused without refitting. The single-machine model was driven by measured terminal voltage and scored against measured-current-derived P/Q, speed, and rotor-to-terminal angle. RMSE is machine pu for P/Q/speed and radians for angle:

| Generator | P RMSE (pu) | Q RMSE (pu) | Speed RMSE (pu) | δ−θ RMSE (rad) |
|---|---:|---:|---:|---:|
| G1 | 0.012384 | 0.013115 | 7.75e−6 | 0.000329 |
| G2 | 0.005540 | 0.006035 | 3.15e−5 | 0.003517 |
| G3 | 0.008089 | 0.006165 | 5.41e−5 | 0.005069 |

The boundary-driven network reconstruction used measured generator terminal voltages and currents, solved the six unmeasured transmission-bus voltages at each of 301 samples, and checked the three generator-bus KCL residuals. With measured currents, RMS residuals were 0.000452, 0.000631, and 0.000414 system pu at buses 1–3 (0.0452, 0.0631, and 0.0414 MVA). Reconstructed mean/min/max voltage magnitudes were B4 1.02172 [1.02161,1.02187], B5 0.99175 [0.99166,0.99189], B6 1.00010 [0.99994,1.00025], B7 1.02279 [1.02270,1.02293], B8 1.01180 [1.01165,1.01195], and B9 1.02854 [1.02828,1.02872] pu. These are algebraic reconstructions, not direct sensor comparisons.

When fitted generator models supply currents instead, RMS KCL residuals rise to 0.04292, 0.01841, and 0.01269 system pu at buses 1–3. The fitted machine trajectories do not yet close the full network DAE accurately on this separate run. Detailed metrics are in `independent_healthy_run_generator_transfer_metrics.csv`, `healthy_network_reconstruction_metrics.csv`, and `healthy_network_reconstructed_bus_voltage_ranges.csv`.

## DAE and checks

`ieee9_dae.py` assembles `dx/dt=f(x,V,u,p)` and `0=Ybus V+Iload(V)−Igen(x,V)`. Machine equations stay on their own bases; current and P/Q outputs convert to the common base. `build_ybus.py` handles line π sections and transformer series, phase tap, and magnetizing branches.

Eight unit tests pass. `validate_ybus_reference.py` compares Ybus with the stored operating-point snapshot and both tap signs. `validate_healthy_run_transfer.py` scores unchanged fitted models on the new no-fault data. `validate_healthy_network_reconstruction.py` reconstructs network voltages and reports algebraic residuals without claiming direct voltage measurements.

**Completion gate:** topology, effective load configuration, selected 100 MVA reference-base conversions, tap sign, full Ybus, algebraic equations, and component tests are checked. Step 3 remains **incomplete** until the generator models reproduce the new run closely enough for a coupled trajectory and transmission-bus measurements or another independent network trajectory verify the reconstructed bus voltages.

## Coupled DAE validation follow-up (5–10 s no-fault candidate)

`run_coupled_ieee9_dae.py` integrated all 18 generator states over the independent healthy candidate interval, solving all nine complex bus voltages from KCL at every RHS evaluation. Previously fitted parameters were fixed. Measured `Efd` and `Pmpu` channels were supplied as time-varying inputs using the existing fitted `kf`/`km` scales; this checks the coupled electrical and machine response under recorded controls, not a free-running governor/AVR model. It used `rtol=2e-6`, `atol=2e-8`, maximum step 1/600 s. A tighter `rtol=2e-8`, `atol=2e-10` run produced a maximum saved-trajectory difference of `4.59e-13` in the output columns, so the reported mismatch is not explained by ODE tolerance. There were 18,002 ODE function evaluations, 18,303 algebraic solves, zero failed solves, and maximum algebraic root infinity residual `4.81e-13 pu`.

That very small KCL residual is the nonlinear algebraic solver's residual because voltages are solved from KCL; it is not an independent validation score. The buswise RMS and maximum KCL magnitude from the saved coupled trajectory are:

| Bus | RMS KCL (pu) | Max KCL (pu) | Healthy voltage status |
|---:|---:|---:|---|
| 1 | 4.02e−14 | 1.51e−13 | Measured generator terminal |
| 2 | 1.49e−14 | 6.40e−14 | Measured generator terminal |
| 3 | 1.20e−14 | 5.01e−14 | Measured generator terminal |
| 4 | 4.23e−15 | 1.25e−14 | Reconstructed; no independent channel |
| 5 | 1.84e−15 | 5.18e−15 | Reconstructed; no independent channel |
| 6 | 1.48e−15 | 3.64e−15 | Reconstructed; no independent channel |
| 7 | 4.05e−15 | 1.33e−14 | Reconstructed; no independent channel |
| 8 | 2.17e−15 | 5.67e−15 | Reconstructed; no independent channel |
| 9 | 3.63e−15 | 1.16e−14 | Reconstructed; no independent channel |

Generator model predictions were compared at the generator terminals with measured positive-sequence channels. Current is reported as exported current; the CSV current sensors point into the machines, so the measurement phasor is negated. The measured-current network reconstruction's low residuals and the fitted-current reconstruction's larger residuals therefore have a consistent sign convention. The transformer current ratio is already represented by winding bases and the Ybus transformer stamp; no extra ratio is applied to generator terminal currents. For each machine, `Ipu_system = Ipu_machine * Srated/100 MVA`, with matching generator terminal voltage bases. These checks and the operating-point tap-sign check give no evidence that a sign, phase-tap, or base conversion explains the coupled discrepancy.

| Generator | Complex current RMSE (system pu) | Magnitude RMSE (system pu) | Current phase RMSE (rad) | Voltage complex RMSE (pu) | Voltage magnitude RMSE (pu) | Voltage phase RMSE (rad) |
|---|---:|---:|---:|---:|---:|---:|
| G1 | 0.06892 | 0.00722 | 0.06861 | 0.06570 | 0.00188 | 0.06325 |
| G2 | 0.10693 | 0.00504 | 0.06113 | 0.06466 | 0.00148 | 0.06315 |
| G3 | 0.05884 | 0.00399 | 0.06275 | 0.06481 | 0.00129 | 0.06331 |

| Generator | P RMSE (machine pu) | Q RMSE (machine pu) | Speed RMSE (pu) | Rotor-relative angle RMSE (rad) |
|---|---:|---:|---:|---:|
| G1 | 0.002792 | 0.006063 | 1.024e−4 | 9.153e−4 |
| G2 | 0.002124 | 0.003562 | 1.094e−4 | 0.003625 |
| G3 | 0.003124 | 0.003974 | 1.172e−4 | 0.004897 |

The complex voltage and current residuals are dominated by phase, while magnitude and P/Q errors are materially smaller. The phase drift is quantitatively linked to the generator speed error: predicted terminal-voltage phase errors begin at 0.0001, 0.0033, and 0.0030 rad for G1–G3, then reach approximately 0.1388 rad for all three at 10 s. Linear phase-error slopes are 0.03001, 0.03038, and 0.03038 rad/s. The respective mean speed errors are `8.04e-5`, `7.91e-5`, and `8.03e-5 pu`; multiplying by synchronous electrical speed `2π·60 rad/s` predicts 0.03033, 0.02984, and 0.03029 rad/s. The agreement identifies accumulated rotor-speed bias in the reused generator trajectories as the main source of the common electrical phase drift. The present results do not isolate whether its upstream cause is the fitted mechanical-power scaling, inertia/damping, other generator-model structure, or a discrepancy in logged speed/control channels; no parameter was changed to force agreement.

Current phase largely follows the voltage phase: the mean of `(predicted current phase error − predicted voltage phase error)` is 0.0117 rad for G1 and −0.0044/−0.0041 rad for G2/G3. Together with the close magnitude and P/Q reconstruction, this indicates the dominant current phasor mismatch is inherited from the rotating terminal reference as the rotor trajectories drift, rather than a separate current sign or winding phase error. The sensor-to-export sign was checked against measured complex power; measured sensors point into the generator while the network injection is outward. The bus-side transformer ratio and winding orientation are already embodied in the transformer stamp and voltage-base selection. No extra generator current ratio or 30-degree rotation belongs between the generator terminal measurement and its connected LV bus. The algebraic root residual is far below the larger measured-current versus fitted-current KCL discrepancy, and the tight-tolerance run reproduces the same trajectory, excluding algebraic and ODE tolerance as its cause.

Rotor angle is compared relative to terminal-voltage phase with the phase branch aligned at the initial sample, which avoids treating arbitrary `2π` wrapping as physical error. Generator terminal voltages/currents are measured; Bus 4–9 voltages are reconstructed in the network solve, not measured or independently validated. Thus the electrical KCL equations and numerical solver pass their internal check, but the independent coupled model comparison does not pass the user's completion gate.

Artifacts: `mbam_package/validation_results/coupled_dae_healthy_candidate_5s_10s.csv`, `coupled_dae_healthy_candidate_metrics.csv`, and `coupled_dae_solver_diagnostics.json`; the tighter run is saved with `_tight` suffixes. The diagnostics explicitly record that no parameter refit occurred.

## Rotor-speed bias investigation

The diagnostic reuses the fitted parameter file unchanged and compares two simulations on the separate no-fault candidate: (1) the single-machine fit with measured terminal voltage imposed and (2) the coupled DAE with terminal voltage solved from network KCL. Both initialize rotor speed from the measured value at 5 s. The first is a generator-model/input check; only the second tests network-generator coupling.

| Generator | Model run | Mean speed error (pu) | Median (pu) | Std. dev. (pu) | RMSE (pu) |
|---|---|---:|---:|---:|---:|
| G1 | Measured terminal V imposed | −1.54e−7 | −5.92e−8 | 7.74e−6 | 7.75e−6 |
| G1 | Coupled DAE | 8.04e−5 | 8.81e−5 | 6.34e−5 | 1.024e−4 |
| G2 | Measured terminal V imposed | −1.34e−6 | 1.64e−6 | 3.15e−5 | 3.15e−5 |
| G2 | Coupled DAE | 7.91e−5 | 9.76e−5 | 7.55e−5 | 1.094e−4 |
| G3 | Measured terminal V imposed | 5.01e−8 | 4.22e−6 | 5.41e−5 | 5.41e−5 |
| G3 | Coupled DAE | 8.03e−5 | 9.97e−5 | 8.53e−5 | 1.172e−4 |

The initial coupled speed error is exactly zero for each generator because `omega(5 s)` is explicitly initialized from the logged speed. It appears at the next sample (16.7 ms later): G1 `+8.82e−6 pu`, G2 `−1.08e−5 pu`, G3 `−4.19e−5 pu`. The first 0.1 s mean errors are respectively `+2.35e−5`, `−2.82e−5`, and `−1.32e−4 pu`, so there is no common positive offset immediately at initialization. The error is transient and time-varying early, then becomes a steadily positive drift: over elapsed intervals 1–3 s, means are `6.17e−5`, `6.23e−5` pu for G2/G3 (G1 `5.44e−5`); over 3–5 s they are `1.46e−4`, `1.47e−4`, and `1.48e−4 pu`. Linear error slopes over the full five seconds are `4.18e−5`, `4.57e−5`, and `4.73e−5 pu/s`. Thus it accumulates after initialization rather than being a starting-speed bias. The measured-voltage simulations do not show this common accumulation.

| Generator | Initial error | First next sample | Mean 0–0.1 s | Mean 0.1–1 s | Mean 1–3 s | Mean 3–5 s | Linear slope (pu/s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| G1 | 0 | 8.82e−6 | 2.35e−5 | −9.45e−7 | 5.44e−5 | 1.46e−4 | 4.18e−5 |
| G2 | 0 | −1.08e−5 | −2.82e−5 | −1.92e−5 | 6.17e−5 | 1.47e−4 | 4.57e−5 |
| G3 | 0 | −4.19e−5 | −1.32e−4 | −3.54e−6 | 6.23e−5 | 1.48e−4 | 4.73e−5 |

### Speed reference, pole pairs, and angle integration

The source SLX machine blocks specify `nPolePairs=20` for G1 and `1` for G2/G3. The logged `w` values have mechanical angular-speed magnitudes near their rated synchronous speeds, so the correct conversions are `omega_sync_mech=2π·60/p` and `w_pu=w/omega_sync_mech`. The logged `theta` is mechanical rotor angle in degrees; electrical rotor angle is `theta_e=p·unwrap(theta·π/180)`. The physical swing-equation angle rate is therefore `d(delta_e)/dt=2π·60·(w_pu−1)`, which is exactly the convention in `mbam_core.py` and `ieee9_dae.py`.

This mapping was checked directly against the 5–10 s Simulink channels. RMS differences between the derivative of converted electrical rotor angle and `2π·60·(w_pu−1)` are `2.86e−5 rad/s` for G1, `7.66e−5 rad/s` for G2, and `5.00e−5 rad/s` for G3. The agreement validates the 20-pole-pair conversion for G1 (not 10 pole pairs) and the one-pair conversions for G2/G3. It is orders of magnitude smaller than the coupled phase drift near `0.03 rad/s`; a speed-reference, pole-pair, degree/radian, or angle-integration scale error is not its cause.

The original logged rotor-angle changes and model rotor-angle changes are compared after removing only the arbitrary initial angle branch/offset. The measured-voltage generator simulations stay close to the logged angle shape (RMSE `0.000329`, `0.00352`, and `0.00507 rad` for G1–G3). The coupled angle errors at 10 s grow to `0.1516`, `0.1491`, and `0.1514 rad`, consistent with integrating their positive late-window speed errors. Trajectory plots are in `mbam_package/validation_results/rotor_speed_bias_trajectories.png` and `rotor_speed_bias_errors.png`.

### Supported source of the coupled-only bias

The speed equation and speed/angle unit mapping do not generate a common bias when the same fitted generator models receive the measured terminal-voltage trajectories. The bias appears only after the network algebraic solver is allowed to choose generator terminal voltages from KCL. At the initial 5 s algebraic solve, the coupled terminal-voltage magnitudes already differ from the measured values by `−0.00304`, `+0.00224`, and `+0.00112 pu` for G1–G3; their phase differences are `+0.00010`, `+0.00329`, and `+0.00299 rad`. The coupled initial electrical-power errors are `−0.00906`, `+0.00531`, and `+0.00500 machine pu` active, and `−0.00916`, `+0.00764`, and `+0.00560 machine pu` reactive. The internal machine states are initialized from fitted `eq0/ed0` and a fitted angle offset, while network voltages are solved freely; those states are not reinitialized to match measured terminal current and the measured-voltage network reconstruction simultaneously.

Before releasing the bus voltages in the DAE, the initialized fitted machine at the **measured** terminal voltage already has current-phasor errors of `0.0178`, `0.0254`, and `0.0464 machine pu` at G1–G3 (phase errors `+0.0409`, `−0.0185`, and `−0.0290 rad`). This directly confirms initial generator current/state inconsistency; the network algebraic solve then shifts terminal voltages to satisfy KCL for those model injections.

This makes **generator/network initial-condition consistency and transfer of fitted internal states** the best-supported origin of the coupled-only discrepancy. It is not a nonzero initial `w` error. The generator equations can reproduce logged speed under imposed measured voltage, but their fitted internal state/current and active-power response do not produce the same initial algebraic operating point when released into the network. The positive average speed bias then accumulates during the network-coupled run. The fitted parameter estimates also have limitations: the original Step 2 report records that G1 and G3 optimizer runs hit their evaluation limits, and all three parameter sets were estimated on the earlier post-disturbance dataset. These facts leave parameter-transfer error as a possible upstream contributor to the initial-state mismatch; the data do not uniquely separate fitted parameter error from initialization mismatch or residual network-model error.

The 10 s endpoint in the no-fault candidate has a small separate signal-processing edge effect: the final 60 Hz phasor window is centered before the last slow-channel timestamp, and the interpolation holds the final phasor instead of extrapolating it. The measured-voltage-only prediction consequently has a last-point speed excursion. Removing that one endpoint changes the coupled mean errors only from `8.04e−5/7.91e−5/8.03e−5` to `8.01e−5/7.88e−5/8.00e−5 pu`; it cannot explain the multi-second drift.

For a second trajectory check, the same fitted models were run over 5–10 s from the older post-disturbance records. Their speed errors remain near zero in both the 5–8 s fit segment and 8–10 s temporal holdout:

| Generator | Segment | Mean error (pu) | Median error (pu) | Std. dev. (pu) |
|---|---|---:|---:|---:|
| G1 | 5–8 s fit segment | −6.95e−7 | −9.22e−7 | 9.94e−6 |
| G1 | 8–10 s holdout | −5.31e−7 | 1.20e−7 | 1.21e−5 |
| G2 | 5–8 s fit segment | 3.87e−6 | 2.23e−6 | 4.72e−5 |
| G2 | 8–10 s holdout | −1.75e−6 | −1.37e−6 | 1.56e−5 |
| G3 | 5–8 s fit segment | 2.27e−6 | 1.31e−7 | 2.91e−5 |
| G3 | 8–10 s holdout | 8.17e−7 | 3.71e−7 | 2.70e−5 |

This old segment follows a disturbance and is not an independent no-fault validation. The three new `healthy_datasets/G*_10s_H.csv` files are generator channels from one separate no-fault candidate simulation, not three independent system runs. No second independent no-fault Simulink trajectory was found in the workspace, so the coupled-only cause is localized to the network-release/initial-condition interface, but not yet isolated to a unique physical parameter or independently replicated across separate no-fault runs.

Detailed samples and summaries are saved in `rotor_speed_bias_timeseries.csv`, `rotor_speed_bias_summary.csv`, `rotor_speed_angle_unit_checks.csv`, and `rotor_speed_bias_legacy_run_summary.csv` under `mbam_package/validation_results/`. No physical parameter was changed during this investigation. **Step 3 remains incomplete**: the coupled speed and rotor-angle trajectories still diverge from the independent candidate, and there is no independent healthy transmission-bus voltage trajectory to validate the network solution.
