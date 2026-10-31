"""Trace measured versus fitted rotor speed/angle on the no-fault candidate.

Uses the fitted parameter file without changing it. Outputs per-sample and
summary diagnostics plus a figure contrasting measured-voltage single-machine
transfer and fully network-coupled trajectories.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
PKG = ROOT / "mbam_package"
OUT = PKG / "validation_results"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(PKG))

from mbam_core import NAMES, SGProblem
from mbam_run import NOMINAL
from validate_healthy_run_transfer import load_new_run, OMEGA, POLES
from validate_generator_models import load_generator_data as load_legacy_run


def fitted_parameters():
    table = pd.read_csv(OUT / "fitted_parameters.csv")
    params = {}
    for g in (1, 2, 3):
        p = NOMINAL[g].copy()
        for row in table.loc[table.generator == f"G{g}"].itertuples(index=False):
            p[row.parameter] = float(row.estimate)
        params[g] = p
    return params


def main():
    pars = fitted_parameters()
    coupled = pd.read_csv(OUT / "coupled_dae_healthy_candidate_5s_10s.csv")
    rows, summaries, angle_checks = [], [], []
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True, constrained_layout=True)
    error_fig, error_axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True, constrained_layout=True)
    legacy_rows = []

    for g in (1, 2, 3):
        d = load_new_run(g)
        t = d["t"]
        c = pars[g]["c"]
        prob = SGProblem(t, d["V"], d["th"], d["vf"], d["Pm"], d["P"], d["Q"],
                         d["w"], d["dmt"], pars[g], OMEGA, salient=(g == 1), ns=1)
        prediction = prob.simulate(np.array([pars[g][name] for name in NAMES]))
        w_fit = prediction[2]
        # In the SG model, predicted rotor angle is terminal phasor angle plus
        # delta-theta. Remove the fitted angular calibration offset for a
        # direct comparison with the logged rotor-angle channel.
        delta_fit = d["th"] + prediction[3] - c
        delta_meas = d["dmt"] + d["th"]
        delta_meas = np.unwrap(delta_meas)
        delta_fit = np.unwrap(delta_fit)
        delta_fit_change = delta_fit - delta_fit[0]
        delta_meas_change = delta_meas - delta_meas[0]

        cpart = coupled[["t_s", f"G{g}_w_pred_pu", f"G{g}_delta_state_rad"]]
        if not np.allclose(cpart.t_s, t, rtol=0, atol=1e-9):
            raise ValueError(f"G{g}: coupled and source data time grids differ")
        w_coupled = cpart[f"G{g}_w_pred_pu"].to_numpy()
        delta_coupled = np.unwrap(cpart[f"G{g}_delta_state_rad"].to_numpy() - c)
        delta_coupled_change = delta_coupled - delta_coupled[0]

        # Verify the source channel meaning: theta is mechanical degrees,
        # w is mechanical rad/s, and electrical angle uses n pole pairs.
        df = pd.read_csv(ROOT / "healthy_datasets" / f"G{g}_10s_H.csv")
        slow = df.dropna(subset=["theta", "w"]).drop_duplicates("t").sort_values("t")
        raw_t = slow.t.to_numpy(float)
        theta_mech = np.unwrap(np.deg2rad(slow.theta.to_numpy(float)))
        w_mech = slow.w.to_numpy(float)
        omega_sync_mech = OMEGA / POLES[g]
        wpu_raw = w_mech / omega_sync_mech
        theta_e = POLES[g] * theta_mech
        dtheta_e = np.gradient(theta_e, raw_t)
        speed_implied_delta_dot = OMEGA * (wpu_raw - 1.0)
        healthy_mask = (raw_t >= 5.0) & (raw_t <= 10.0)
        angle_checks.append({
            "generator": f"G{g}", "pole_pairs_from_slx": POLES[g],
            "omega_sync_mechanical_rad_s": omega_sync_mech,
            "theta_e_dot_vs_omega_delta_rmse_rad_s": float(np.sqrt(np.mean(
                (dtheta_e[healthy_mask] - speed_implied_delta_dot[healthy_mask]) ** 2))),
            "theta_e_dot_mean_rad_s": float(np.mean(dtheta_e[healthy_mask])),
            "omega_delta_mean_rad_s": float(np.mean(speed_implied_delta_dot[healthy_mask])),
        })

        models = {"measured_voltage_fit": w_fit, "coupled_dae_fit": w_coupled}
        error_vectors = {}
        for model, wpred in models.items():
            err = wpred - d["w"]
            error_vectors[model] = err
            slope = float(np.polyfit(t - t[0], err, 1)[0])
            summary = {
                "generator": f"G{g}", "model": model, "n": len(err),
                "mean_error_pu": float(np.mean(err)),
                "median_error_pu": float(np.median(err)),
                "std_error_pu": float(np.std(err)),
                "rmse_error_pu": float(np.sqrt(np.mean(err**2))),
                "min_error_pu": float(np.min(err)), "max_error_pu": float(np.max(err)),
                "initial_error_pu": float(err[0]),
                "first_sample_after_init_error_pu": float(err[1]),
                "first_0p1s_mean_error_pu": float(np.mean(err[(t-t[0]) <= .1])),
                "linear_error_slope_pu_per_s": slope,
                "mean_error_0_to_0p1s_pu": float(np.mean(err[(t-t[0]) <= .1])),
                "mean_error_0p1_to_1s_pu": float(np.mean(err[(t-t[0] > .1) & (t-t[0] <= 1)])),
                "mean_error_1_to_3s_pu": float(np.mean(err[(t-t[0] > 1) & (t-t[0] <= 3)])),
                "mean_error_3_to_5s_pu": float(np.mean(err[(t-t[0] > 3) & (t-t[0] <= 5)])),
                "predicted_angle_drift_rate_rad_s_from_mean_speed_error": float(OMEGA*np.mean(err)),
            }
            summaries.append(summary)
            for k in range(len(t)):
                rows.append({"generator": f"G{g}", "t_s": t[k], "elapsed_s": t[k]-t[0],
                             "w_measured_pu": d["w"][k],
                             "w_measured_mechanical_rad_s": d["w"][k]*omega_sync_mech,
                             "w_measured_electrical_rad_s": d["w"][k]*OMEGA,
                             "w_fitted_pu": wpred[k], "speed_error_pu": err[k],
                             "model": model})

        # Angle summaries are referenced to their own t=5 s values so only
        # accumulated trajectory error, not fitted offset c, is compared.
        angle_models = {
            "measured_voltage_fit": delta_fit_change,
            "coupled_dae_fit": delta_coupled_change,
        }
        for model, angle_pred in angle_models.items():
            ae = angle_pred - delta_meas_change
            ss = next(s for s in summaries if s["generator"] == f"G{g}" and s["model"] == model)
            ss["rotor_angle_change_rmse_rad"] = float(np.sqrt(np.mean(ae**2)))
            ss["rotor_angle_change_final_error_rad"] = float(ae[-1])

        # The older post-disturbance record is included as a separate
        # trajectory check, but is not an independent no-fault simulation.
        legacy, _ = load_legacy_run(g)
        legacy_prob = SGProblem(legacy["t"], legacy["V"], legacy["th"], legacy["vf"],
                                legacy["Pm"], legacy["P"], legacy["Q"], legacy["w"],
                                legacy["dmt"], pars[g], OMEGA, salient=(g == 1), ns=1)
        legacy_w = legacy_prob.simulate(np.array([pars[g][name] for name in NAMES]))[2]
        legacy_err = legacy_w - legacy["w"]
        for split, mask in (("5_to_8_fit_segment", legacy["t"] < 8.0),
                            ("8_to_10_temporal_holdout", legacy["t"] >= 8.0)):
            es = legacy_err[mask]
            legacy_rows.append({
                "generator": f"G{g}", "dataset": "legacy_post_disturbance_5s_onward",
                "split": split, "n": int(mask.sum()),
                "mean_error_pu": float(np.mean(es)), "median_error_pu": float(np.median(es)),
                "std_error_pu": float(np.std(es)), "rmse_error_pu": float(np.sqrt(np.mean(es**2))),
                "min_error_pu": float(np.min(es)), "max_error_pu": float(np.max(es)),
            })

        axw, axa = axes[g-1]
        axw.plot(t, d["w"], color="black", linewidth=1.5, label="Simulink logged")
        axw.plot(t, w_fit, label="Fitted machine, measured V")
        axw.plot(t, w_coupled, linestyle="--", label="Fitted machine, coupled DAE")
        axw.set_ylabel(f"G{g} speed (pu)")
        axw.grid(alpha=.25)
        axw.legend(fontsize=8, loc="best")
        axa.plot(t, delta_meas_change, color="black", linewidth=1.5, label="Simulink rotor angle change")
        axa.plot(t, delta_fit_change, label="Fitted machine angle change")
        axa.plot(t, delta_coupled_change, linestyle="--", label="Coupled DAE angle change")
        axa.set_ylabel(f"G{g} Δ rotor angle (rad)")
        axa.grid(alpha=.25)
        axa.legend(fontsize=8, loc="best")

        axe_w, axe_a = error_axes[g-1]
        axe_w.plot(t, error_vectors["measured_voltage_fit"], label="Fitted machine, measured V")
        axe_w.plot(t, error_vectors["coupled_dae_fit"], linestyle="--", label="Fitted machine, coupled DAE")
        axe_w.axhline(0, color="black", linewidth=.8)
        axe_w.set_ylabel(f"G{g} speed error (pu)")
        axe_w.grid(alpha=.25)
        axe_w.legend(fontsize=8, loc="best")
        fit_angle_error = delta_fit_change-delta_meas_change
        coupled_angle_error = delta_coupled_change-delta_meas_change
        axe_a.plot(t, fit_angle_error, label="Fitted machine, measured V")
        axe_a.plot(t, coupled_angle_error, linestyle="--", label="Fitted machine, coupled DAE")
        axe_a.axhline(0, color="black", linewidth=.8)
        axe_a.set_ylabel(f"G{g} rotor angle error (rad)")
        axe_a.grid(alpha=.25)
        axe_a.legend(fontsize=8, loc="best")

    axes[-1, 0].set_xlabel("Time (s)")
    axes[-1, 1].set_xlabel("Time (s)")
    fig.suptitle("Measured and fitted rotor speed / angle, independent no-fault candidate")
    fig.savefig(OUT / "rotor_speed_bias_trajectories.png", dpi=180)
    plt.close(fig)
    error_axes[-1, 0].set_xlabel("Time (s)")
    error_axes[-1, 1].set_xlabel("Time (s)")
    error_fig.suptitle("Speed and rotor-angle errors on the independent no-fault candidate")
    error_fig.savefig(OUT / "rotor_speed_bias_errors.png", dpi=180)
    plt.close(error_fig)
    pd.DataFrame(rows).to_csv(OUT / "rotor_speed_bias_timeseries.csv", index=False)
    pd.DataFrame(summaries).to_csv(OUT / "rotor_speed_bias_summary.csv", index=False)
    pd.DataFrame(angle_checks).to_csv(OUT / "rotor_speed_angle_unit_checks.csv", index=False)
    pd.DataFrame(legacy_rows).to_csv(OUT / "rotor_speed_bias_legacy_run_summary.csv", index=False)
    print(pd.DataFrame(summaries).to_string(index=False, float_format=lambda x: f"{x:.7g}"))
    print("\nAngle/speed unit checks")
    print(pd.DataFrame(angle_checks).to_string(index=False, float_format=lambda x: f"{x:.7g}"))
    print("\nLegacy post-disturbance source trajectory (not independent no-fault)")
    print(pd.DataFrame(legacy_rows).to_string(index=False, float_format=lambda x: f"{x:.7g}"))


if __name__ == "__main__":
    main()
