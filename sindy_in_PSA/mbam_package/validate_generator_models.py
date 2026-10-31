"""Fit the full physical generator model on 5–8 s; evaluate on 8–10 s.

This is Step 2 (baseline physical-model validation), before MBAM reduction.
Run after prepare_healthy_data.py:
    python sindy_in_PSA/mbam_package/validate_generator_models.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_DIR = Path(__file__).resolve().parent
DATA_DIR = ROOT / "prepared_data" / "healthy_5s_onward"
OUTPUT_DIR = PACKAGE_DIR / "validation_results"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(PACKAGE_DIR))

from mbam_core import DESC, IDX, NAMES, SGProblem  # noqa: E402
from mbam_run import NOMINAL  # noqa: E402
from prepare_healthy_data import positive_sequence_phasors  # noqa: E402

F0 = 60.0
OMEGA_B = 2.0 * np.pi * F0
TRAIN_END = 8.0
SN = {1: 247.5e6, 2: 230.0e6, 3: 128.0e6}
VN_LL = {1: 16.5e3, 2: 18.0e3, 3: 13.8e3}
POLE_PAIRS = {1: 20, 2: 1, 3: 1}
GEN_SLOW_COLS = {
    1: {"pm": "Pm_pu", "theta": "theta", "w": "w", "vf": "Excitation V"},
    2: {"pm": "Pmpu", "theta": "theta", "w": "w", "vf": "excitation V"},
    3: {"pm": "Pmpu", "theta": "theta", "w": "w", "vf": "excitation V"},
}


def complex_column(df, prefix):
    return df[f"{prefix}_re"].to_numpy(float) + 1j * df[f"{prefix}_im"].to_numpy(float)


def interpolate_complex(t_new, t, z):
    return np.interp(t_new, t, z.real) + 1j * np.interp(t_new, t, z.imag)


def load_generator_data(gen):
    wave = pd.read_csv(DATA_DIR / f"G{gen}_waveforms_raw_units.csv")
    slow = pd.read_csv(DATA_DIR / f"G{gen}_slow_channels_raw_units.csv")
    ph = positive_sequence_phasors(wave, "voltage")
    iph = positive_sequence_phasors(wave, "current")
    ph = ph.merge(iph, on="t_center_s", validate="one_to_one")
    t_ph = ph["t_center_s"].to_numpy(float)
    t = slow["t"].to_numpy(float)
    c = GEN_SLOW_COLS[gen]
    v1 = complex_column(ph, "voltage_positive")
    i1 = complex_column(ph, "current_positive")
    v1_t = interpolate_complex(t, t_ph, v1)
    i1_t = interpolate_complex(t, t_ph, i1)

    # The current sign follows generator export convention. G2's calculated
    # terminal P/Q is independently checked against the B2 logged P/Q below.
    s_mva = -3.0 * v1_t * np.conj(i1_t)
    p_pu = s_mva.real / SN[gen]
    q_pu = s_mva.imag / SN[gen]
    w_pu = slow[c["w"]].to_numpy(float) / (OMEGA_B / POLE_PAIRS[gen])
    delta = np.unwrap(np.deg2rad(slow[c["theta"]].to_numpy(float))) * POLE_PAIRS[gen]
    theta_v = np.unwrap(np.angle(v1_t))
    dmt = delta - theta_v
    keep = (t >= 5.0) & (t <= 10.0)

    data = {
        "t": t[keep],
        "V": np.abs(v1_t[keep]) / (VN_LL[gen] / np.sqrt(3.0)),
        "th": theta_v[keep],
        "vf": slow[c["vf"]].to_numpy(float)[keep],
        "Pm": slow[c["pm"]].to_numpy(float)[keep],
        "P": p_pu[keep],
        "Q": q_pu[keep],
        "w": w_pu[keep],
        "dmt": dmt[keep],
    }
    return data, ph


def evaluate_predictions(prob, p):
    pred = list(prob.simulate(p))
    # c is the fitted sensor/model angle offset. Return the rotor-angle output
    # to the measured reference frame before scoring against measured dmt.
    pred[3] = pred[3] - p[IDX["c"]]
    observed = (prob.y["P"], prob.y["Q"], prob.y["w"], prob.y["dmt"])
    names = ("P_pu", "Q_pu", "w_pu", "delta_minus_theta_rad")
    train = prob.fit_mask
    validation = ~train
    records = []
    for name, y, yhat in zip(names, observed, pred):
        for split, mask in (("train", train), ("validation", validation)):
            err = yhat[mask] - y[mask]
            rmse = float(np.sqrt(np.mean(err**2)))
            scale = float(np.ptp(y[train]))
            nrmse = rmse / max(scale, 1e-9)
            records.append(
                {"output": name, "split": split, "n": int(mask.sum()), "rmse": rmse,
                 "normalization_scale_train_range": scale, "nrmse": nrmse}
            )
    return pred, pd.DataFrame(records)


def main():
    if not DATA_DIR.exists():
        raise FileNotFoundError("Run prepare_healthy_data.py before this validation script")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_metrics = []
    param_rows = []
    plot, axes = plt.subplots(3, 4, figsize=(15, 10), constrained_layout=True)
    output_names = ("P", "Q", "w", "dmt")
    output_labels = ("P (pu)", "Q (pu)", "Speed (pu)", "delta - theta (rad)")
    result_log = {}

    for gen in (1, 2, 3):
        print(f"G{gen}: assembling data and fitting physical model...", flush=True)
        d, ph = load_generator_data(gen)
        fit_mask = d["t"] < TRAIN_END
        base = dict(NOMINAL[gen])
        base.update(
            kf=1.8 / max(float(d["vf"][0]), 1e-9),
            km=1.0,
            c=0.0,
            eq0=1.1,
            ed0=0.0 if gen == 1 else 0.3,
        )
        p0 = np.array([base[n] for n in NAMES], dtype=float)
        prob = SGProblem(
            d["t"], d["V"], d["th"], d["vf"], d["Pm"], d["P"], d["Q"], d["w"], d["dmt"],
            base, OMEGA_B, salient=(gen == 1), ns=1, sig_frac=0.005, fit_mask=fit_mask,
        )
        nuisance = [n for n in ("c", "kf", "km", "eq0", "ed0", "H2", "D") if n not in prob.fixed]
        p_stage1, stage1_cost, stage1 = prob.fit(p0, free=nuisance, max_nfev=120)
        pfit, cost, solution = prob.fit(
            p_stage1, max_nfev=400, xtol=1e-8, ftol=1e-8, gtol=1e-8
        )
        predictions, metrics = evaluate_predictions(prob, pfit)
        metrics.insert(0, "generator", f"G{gen}")
        all_metrics.append(metrics)
        for name, value in zip(NAMES, pfit):
            param_rows.append(
                {"generator": f"G{gen}", "parameter": name, "description": DESC.get(name, name),
                 "estimate": float(value)}
            )

        result_log[f"G{gen}"] = {
            "success": bool(solution.success),
            "message": str(solution.message),
            "fit_cost": float(cost),
            "fit_nfev": int(solution.nfev),
            "stage1_success": bool(stage1.success),
            "stage1_nfev": int(stage1.nfev),
            "stage1_cost": float(stage1_cost),
            "training_samples": int(fit_mask.sum()),
            "validation_samples": int((~fit_mask).sum()),
        }
        print(
            f"G{gen}: optimizer success={solution.success}, nfev={solution.nfev}, "
            f"training cost={cost:.4g}", flush=True
        )
        y_values = (prob.y["P"], prob.y["Q"], prob.y["w"], prob.y["dmt"])
        for col, (name, label, y, yhat) in enumerate(zip(output_names, output_labels, y_values, predictions)):
            ax = axes[gen - 1, col]
            ax.plot(d["t"], y, label="measured/derived", linewidth=1.2)
            ax.plot(d["t"], yhat, label="physical model", linewidth=1.1, linestyle="--")
            ax.axvline(TRAIN_END, color="black", linestyle=":", linewidth=1)
            ax.set_title(f"G{gen}: {label}")
            ax.grid(alpha=0.25)
            if gen == 3:
                ax.set_xlabel("Time (s)")
            if col == 0:
                ax.set_ylabel("Value")
            if gen == 1 and col == 0:
                ax.legend(fontsize=8)

    plot.savefig(OUTPUT_DIR / "healthy_generator_model_validation.png", dpi=160)
    plt.close(plot)
    metric_frame = pd.concat(all_metrics, ignore_index=True)
    metric_frame.to_csv(OUTPUT_DIR / "metrics.csv", index=False)
    pd.DataFrame(param_rows).to_csv(OUTPUT_DIR / "fitted_parameters.csv", index=False)
    (OUTPUT_DIR / "run_summary.json").write_text(json.dumps(result_log, indent=2), encoding="utf-8")

    # Independent logged P/Q cross-check for G2 against bus B2.
    bus = pd.read_csv(ROOT / "fault datasets" / "IVPQ_B2_H.csv")
    bus = bus.loc[(bus.t >= 5.0) & bus.P.notna() & bus.Q.notna()].sort_values("t")
    g2_ph = pd.read_csv(DATA_DIR / "G2_sequence_phasors_raw_units.csv")
    gv = complex_column(g2_ph, "voltage_positive")
    gi = complex_column(g2_ph, "current_positive")
    gs = -3.0 * gv * np.conj(gi)
    bus_power_estimate = np.interp(bus.t, g2_ph.t_center_s, gs.real) + 1j * np.interp(
        bus.t, g2_ph.t_center_s, gs.imag
    )
    bus_power_measured = bus.P.to_numpy(float) + 1j * bus.Q.to_numpy(float)
    sign_check = {
        "bus": "B2",
        "generator": "G2",
        "samples": int(len(bus)),
        "median_relative_complex_power_error": float(
            np.median(np.abs(bus_power_estimate - bus_power_measured) / np.maximum(np.abs(bus_power_measured), 1e-9))
        ),
        "median_P_measured_W": float(np.median(bus_power_measured.real)),
        "median_P_from_generator_phasors_W": float(np.median(bus_power_estimate.real)),
        "current_sign_used": "negative of 3*V1*conj(I1), matching generator export at G2/B2",
    }
    (OUTPUT_DIR / "G2_B2_power_sign_check.json").write_text(json.dumps(sign_check, indent=2), encoding="utf-8")

    report = [
        "# Physical Generator Model Validation (Step 2)",
        "",
        "The full two-axis/subtransient physical model was fit on 5.0 <= t < 8.0 s and evaluated on 8.0 <= t <= 10.0 s. Only the fitting segment contributes to the optimizer residuals and robust residual scales. The validation segment is simulated continuously from the fitted model initial state without parameter fitting.",
        "",
        "This is a same-run temporal validation. It does not establish generalization to an independent healthy simulation. The P/Q targets are derived from the terminal positive-sequence phasors, not independently measured for all three generator terminals. G2's derived power is cross-checked against the logged B2 P/Q channels.",
        "",
        "Per-unit bases used here are the nominal generator ratings in the package (247.5, 230, and 128 MVA) and nominal terminal voltage bases (16.5, 18, and 13.8 kV line-to-line). Speed is divided by nominal synchronous mechanical speed using pole-pair counts (20, 1, 1); rotor angle in degrees is converted to electrical radians using the same counts. These mappings are provisional pending confirmation against the source Simulink blocks.",
        "",
        "A provisional trajectory acceptance criterion for this run is validation NRMSE <= 0.15 for each output, normalized by that output's training range. Passing this criterion supports predictive fit on the same-run holdout only; it does not establish parameter identifiability or independent-run generalization.",
        "",
        f"Current result: {'all outputs pass' if metric_frame.loc[metric_frame.split == 'validation', 'nrmse'].max() <= 0.15 else 'one or more outputs fail'} the provisional criterion; maximum validation NRMSE is {metric_frame.loc[metric_frame.split == 'validation', 'nrmse'].max():.5f}. The optimizer converged for {', '.join(k for k,v in result_log.items() if v['success']) or 'none'}; parameters for other generators are provisional because they hit the evaluation limit.",
        "",
        "## Fit status",
        "",
        "```json",
        json.dumps(result_log, indent=2),
        "```",
        "",
        "## Metrics",
        "",
        "`metrics.csv` contains RMSE and NRMSE. NRMSE is normalized by the training-segment range for each output.",
        "",
        "```text",
        metric_frame.to_string(index=False, float_format=lambda x: f"{x:.5g}"),
        "```",
        "",
        "## Limitations",
        "",
        "- Bus B2 P/Q independently supports the selected export current sign for G2; equivalent logged terminal P/Q cross-checks are unavailable for G1 and G3 in the current data bundle.",
        "- Excitation and mechanical-power channel units/scales, nominal machine parameters, and generator voltage/power bases still require confirmation from the original model/export configuration.",
        "- Do not run MBAM reduction or report a validated reduced model unless the full-model baseline is acceptable on the held-out interval.",
        "",
        "The fit/prediction trajectories are in `healthy_generator_model_validation.png`; parameter estimates are in `fitted_parameters.csv`.",
        "",
    ]
    (OUTPUT_DIR / "VALIDATION_REPORT.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Validation outputs written to {OUTPUT_DIR}")
    print(metric_frame.to_string(index=False, float_format=lambda x: f"{x:.5g}"))
    print("G2/B2 power sign check:", json.dumps(sign_check))


if __name__ == "__main__":
    main()
