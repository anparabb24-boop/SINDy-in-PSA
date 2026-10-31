"""Healthy-only robust STLSQ models and generator-level binary fault pilot.

Run: python sindy_in_PSA/generator_fault_experiment.py

This is a data-driven residual experiment. It does not modify or parameterize
the physical DAE, and it is neither MBAM nor constrained SINDy.
"""
from __future__ import annotations

from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generator_fault_experiment_results"
sys.path.insert(0, str(ROOT))
from prepare_healthy_data import split_merged_streams, positive_sequence_phasors  # noqa: E402

FIT = (5.0, 8.0)
CAL = (8.0, 10.0)
NOMINAL_EVENT = (2.0, 2.3)
F0 = 60.0
THRESHOLD_QUANTILE = 0.99
STLS_THRESHOLD = 0.08
HUBER_DELTA = 1.5
RIDGE = 1e-4
FEATURES = ["omega_dev", "theta_rel_rad", "Pm_dev", "Efd_dev"]
TARGETS = ["theta_dot_rad_s", "omega_pu_dot_s"]
GEN = {
    1: {"healthy": "healthy_datasets/G1_10s_H.csv", "fault": "fault datasets/G1_Fag_10s.csv", "p": "Pmpu", "pf": "Pm_pu", "e": "Efd", "ef": "Excitation V"},
    2: {"healthy": "healthy_datasets/G2_10s_H.csv", "fault": "fault datasets/G2_Fag_10s.csv", "p": "Pmpu", "pf": "Pmpu", "e": "Efd", "ef": "excitation V"},
    3: {"healthy": "healthy_datasets/G3_10s_H.csv", "fault": "fault datasets/G3_Fag_10s.csv", "p": "Pmpu", "pf": "Pmpu", "e": "Efd", "ef": "excitation V"},
}


def robust_scale(x):
    x = np.asarray(x, float)
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def split_data(path, pcol, ecol):
    """Use the shared stream splitter; preserve native slow samples and audit input."""
    raw = pd.read_csv(path)
    cols = ["t", "Va", "Vb", "Vc", "Ia", "Ib", "Ic", "theta", "w", pcol, ecol]
    if missing := sorted(set(cols) - set(raw.columns)):
        raise ValueError(f"{path}: missing columns {missing}")
    # Audit source before stream separation.
    raw["t"] = pd.to_numeric(raw.t, errors="coerce")
    if raw.t.isna().any():
        raise ValueError(f"{path}: missing/non-numeric timestamps")
    if not np.isfinite(raw.t).all():
        raise ValueError(f"{path}: nonfinite timestamps")
    # Source merge intentionally places fast and slow records at identical
    # timestamps. The shared splitter disambiguates these by populated channels.
    raw_collision_count = int(raw.t.duplicated(keep=False).sum())
    fast, slow = split_merged_streams(path, ("theta", "w", pcol, ecol))
    if slow.t.duplicated().any() or fast.t.duplicated().any():
        raise ValueError(f"{path}: duplicate timestamps after split")
    raw.attrs["timestamp_collision_rows"] = raw_collision_count
    return raw, fast, slow


def derivative(slow, pcol, ecol, w0, theta0):
    """7-point quadratic Savitzky-Golay derivative on the native 60 Hz grid.

    theta is treated as a degree-valued exported angle (confirmed by its
    multi-hundred-unit range); convert to radians, unwrap, then differentiate.
    The centered 7-point stencil spans 100 ms, so it smears sharp transients.
    """
    t = slow.t.to_numpy(float)
    dt = float(np.median(np.diff(t)))
    if not np.allclose(np.diff(t), dt, rtol=1e-5, atol=1e-9) or not (1/61 < dt < 1/59):
        raise ValueError(f"slow stream is not uniform native ~60 Hz: median dt={dt}")
    theta = np.unwrap(np.deg2rad(slow.theta.to_numpy(float)))
    omega = slow.w.to_numpy(float) / w0
    smooth_theta = savgol_filter(theta, 7, 2, mode="interp")
    theta_dot = savgol_filter(theta, 7, 2, deriv=1, delta=dt, mode="interp")
    omega_dot = savgol_filter(omega, 7, 2, deriv=1, delta=dt, mode="interp")
    features = pd.DataFrame({
        "t": t,
        "omega_dev": omega - 1.0,
        "theta_rel_rad": smooth_theta - theta0,
        "Pm_dev": slow[pcol].to_numpy(float) - float(np.median(slow[pcol][(t >= FIT[0]) & (t < FIT[1])])),
        "Efd_dev": slow[ecol].to_numpy(float) - float(np.median(slow[ecol][(t >= FIT[0]) & (t < FIT[1])])),
        "theta_dot_rad_s": theta_dot,
        "omega_pu_dot_s": omega_dot,
    })
    return features


def huber_fit(A, y):
    weights = np.ones(len(y))
    penalty = np.eye(A.shape[1]); penalty[0, 0] = 0
    coef = np.zeros(A.shape[1])
    for _ in range(40):
        sw = np.sqrt(weights)
        c = np.linalg.solve((A * sw[:, None]).T @ (A * sw[:, None]) + RIDGE * penalty,
                            (A * sw[:, None]).T @ (y * sw))
        residual = y - A @ c
        scale = max(robust_scale(residual), 1e-10)
        u = np.abs(residual) / (HUBER_DELTA * scale)
        nw = np.where(u > 1, 1 / np.maximum(u, 1e-12), 1.0)
        if np.linalg.norm(c-coef) < 1e-8*(1+np.linalg.norm(coef)):
            coef = c; break
        coef, weights = c, nw
    return coef


def fit_stlsq(X, y):
    """Robust STLSQ in standardized feature/target coordinates."""
    mu, sd = X.mean(0), X.std(0)
    valid = sd > 1e-9; sd[~valid] = 1
    Z = (X-mu)/sd; Z[:, ~valid] = 0
    ymu, ysd = float(np.median(y)), max(robust_scale(y), 1e-7)
    yn = (y-ymu)/ysd
    A = np.column_stack([np.ones(len(Z)), Z])
    support = np.r_[True, valid]
    for _ in range(15):
        ids = np.flatnonzero(support)
        c = huber_fit(A[:, ids], yn)
        keep = np.abs(c) >= STLS_THRESHOLD
        keep[0] = True
        ns = np.zeros_like(support); ns[ids[keep]] = True
        if np.array_equal(ns, support): break
        support = ns
    ids = np.flatnonzero(support); c = huber_fit(A[:, ids], yn)
    coef = np.zeros(A.shape[1]); coef[ids] = c
    return {"mu": mu, "sd": sd, "valid": valid, "ymu": ymu, "ysd": ysd, "coef": coef, "support": support}


def predict(model, X):
    Z = (X-model["mu"])/model["sd"]
    Z[:, ~model["valid"]] = 0
    return (np.column_stack([np.ones(len(Z)), Z]) @ model["coef"])*model["ysd"] + model["ymu"]


def residual_metric(residual, scales):
    """Dimensionless Euclidean norm of the two target residuals scaled by
    healthy-fit robust derivative scales; identical for cal and fault data."""
    return np.sqrt(np.sum((residual / scales[None, :])**2, axis=1) / residual.shape[1])


def calibrate_threshold(scores):
    scores = np.asarray(scores, float)
    if not len(scores) or not np.isfinite(scores).all(): raise ValueError("invalid calibration scores")
    return float(np.quantile(scores, THRESHOLD_QUANTILE, method="higher"))


def make_voltage_baseline(g, healthy_fast, fault_fast):
    h = positive_sequence_phasors(healthy_fast, "voltage")
    f = positive_sequence_phasors(fault_fast, "voltage")
    # Discard endpoint/edge windows whose full cycle is not strictly inside data.
    h = h[(h.t_center_s >= CAL[0]) & (h.t_center_s < CAL[1] - 1/F0)]
    f = f[(f.t_center_s >= 1/F0) & (f.t_center_s < 10.0 - 1/F0)]
    hr = h.voltage_negative_magnitude.to_numpy()/np.maximum(h.voltage_positive_magnitude.to_numpy(), 1e-9)
    fr = f.voltage_negative_magnitude.to_numpy()/np.maximum(f.voltage_positive_magnitude.to_numpy(), 1e-9)
    threshold = calibrate_threshold(hr)
    return h, f, hr, fr, threshold


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    all_metrics, event_evidence, coeff_rows = [], [], []
    record = {}
    for gen, cfg in GEN.items():
        hp = ROOT / cfg["healthy"]; fp = ROOT / cfg["fault"]
        if not hp.exists() or not fp.exists(): raise FileNotFoundError(f"G{gen}: {hp} or {fp} missing")
        hraw, hw, hs = split_data(hp, cfg["p"], "Efd")
        fraw, fw, fs = split_data(fp, cfg["pf"], cfg["ef"])
        # Existing G2/G3 healthy files use Efd too.
        if cfg["e"] != "Efd":
            # healthy file schema is Efd; update alias for consistent feature definition
            hs = hs.rename(columns={"Efd": cfg["e"]})
        if hs.t.between(FIT[0], FIT[1], inclusive="left").sum() < 100:
            raise ValueError(f"G{gen}: insufficient healthy fit samples")
        w0 = float(np.median(hs.loc[hs.t.between(FIT[0], FIT[1], inclusive="left"), "w"]))
        theta0 = float(np.deg2rad(hs.loc[hs.t.between(FIT[0], FIT[1], inclusive="left"), "theta"].iloc[0]))
        h = derivative(hs.rename(columns={cfg["e"]: "Efd"}), cfg["p"], "Efd", w0, theta0)
        f = derivative(fs.rename(columns={cfg["ef"]: "Efd", cfg["pf"]: cfg["p"]}), cfg["p"], "Efd", w0, theta0)
        tr = (h.t >= FIT[0]) & (h.t < FIT[1])
        # Three sample intervals equal the centered derivative's half-support
        # (50 ms), and exceed the requested one-cycle guard at 60 Hz.
        derivative_guard = 3.0 * float(np.median(np.diff(h.t)))
        cal = (h.t >= CAL[0]+derivative_guard) & (h.t < CAL[1]-1/F0)
        # Model on 5-8s healthy only. No fault sample enters any fitting/calibration step.
        Xh = h[FEATURES].to_numpy(float); Xf = f[FEATURES].to_numpy(float)
        models = [fit_stlsq(Xh[tr], h.loc[tr, target].to_numpy(float)) for target in TARGETS]
        preds_h = np.column_stack([predict(m, Xh) for m in models])
        preds_f = np.column_stack([predict(m, Xf) for m in models])
        yh = h[TARGETS].to_numpy(float); yf = f[TARGETS].to_numpy(float)
        scales = np.maximum(np.array([robust_scale(yh[tr, j]) for j in range(2)]), 1e-8)
        rh, rf = yh-preds_h, yf-preds_f
        sh, sf = residual_metric(rh, scales), residual_metric(rf, scales)
        threshold = calibrate_threshold(sh[cal.to_numpy()])
        # One-cycle debounce: at 60 Hz, two consecutive samples are approximately one cycle.
        alarm = sf > threshold
        sustained = np.convolve(alarm.astype(int), np.ones(2, int), mode="valid") == 2
        first = int(np.flatnonzero(sustained)[0]) if sustained.any() else None
        # Event evidence independently derives boundaries from voltage sequence ratio threshold.
        ph_h, ph_f, hr, fr, vth = make_voltage_baseline(gen, hw, fw)
        event_candidate = fr > vth
        active = event_candidate & (ph_f.t_center_s >= 1.5) & (ph_f.t_center_s <= 3.0)
        if active.any():
            ti = ph_f.t_center_s.to_numpy()
            idx = np.flatnonzero(active); measured_start, measured_end = float(ti[idx[0]]), float(ti[idx[-1]])
        else:
            measured_start, measured_end = NOMINAL_EVENT
        # Fault event detection: first two-sample residual alarm in measured event interior.
        measured_window = (f.t >= measured_start+derivative_guard) & (f.t <= measured_end-derivative_guard)
        fi = np.flatnonzero(measured_window.to_numpy())
        detected_idx = next((i for i in fi if i+1 < len(alarm) and alarm[i] and alarm[i+1]), None)
        delay = float(f.t.iloc[detected_idx+1]-measured_start) if detected_idx is not None else np.nan
        cal_scores = sh[cal.to_numpy()]
        fault_scores = sf[measured_window.to_numpy()]
        false_rate = float(np.mean(cal_scores > threshold))
        metrics = {
            "generator": f"G{gen}", "fit_samples": int(tr.sum()), "calibration_samples": int(cal.sum()),
            "fault_event_start_s": measured_start, "fault_event_end_s": measured_end,
            "threshold_q": THRESHOLD_QUANTILE, "residual_threshold": threshold,
            "healthy_calibration_median": float(np.median(cal_scores)), "healthy_calibration_p95": float(np.quantile(cal_scores,.95)),
            "healthy_calibration_p99": float(np.quantile(cal_scores,.99)), "healthy_false_alarm_rate": false_rate,
            "fault_residual_median": float(np.median(fault_scores)) if len(fault_scores) else np.nan,
            "fault_residual_p95": float(np.quantile(fault_scores,.95)) if len(fault_scores) else np.nan,
            "fault_alarm_fraction": float(np.mean(fault_scores > threshold)) if len(fault_scores) else np.nan,
            "detected": detected_idx is not None, "detection_delay_s": delay,
            "voltage_ratio_threshold": vth, "voltage_ratio_cal_median": float(np.median(hr)),
            "voltage_ratio_fault_median": float(np.median(fr[(ph_f.t_center_s >= measured_start+1/F0)&(ph_f.t_center_s <= measured_end-1/F0)])),
            "voltage_ratio_event_detected": bool(active.any()),
            "voltage_ratio_false_alarm_rate": float(np.mean(hr > vth)),
            "voltage_ratio_alarm_fraction_event": float(np.mean(event_candidate[(ph_f.t_center_s >= measured_start+1/F0)&(ph_f.t_center_s <= measured_end-1/F0)])),
            "healthy_source": str(hp.relative_to(ROOT)), "fault_source": str(fp.relative_to(ROOT)),
        }
        all_metrics.append(metrics)
        hpred = pd.DataFrame({"t": h.t, "theta_dot_true": yh[:,0], "theta_dot_pred": preds_h[:,0], "theta_dot_residual": rh[:,0],
                              "omega_dot_true": yh[:,1], "omega_dot_pred": preds_h[:,1], "omega_dot_residual": rh[:,1], "residual_score": sh,
                              "partition": np.where(tr,"fit",np.where(cal,"calibration","unused")), "alarm": sh>threshold})
        fpred = pd.DataFrame({"t": f.t, "theta_dot_true": yf[:,0], "theta_dot_pred": preds_f[:,0], "theta_dot_residual": rf[:,0],
                              "omega_dot_true": yf[:,1], "omega_dot_pred": preds_f[:,1], "omega_dot_residual": rf[:,1], "residual_score": sf,
                              "partition": np.where(measured_window,"fault_event", "outside_event"), "alarm": alarm})
        hpred.to_csv(OUT/f"G{gen}_healthy_predictions.csv", index=False)
        fpred.to_csv(OUT/f"G{gen}_fault_predictions.csv", index=False)
        for j, target in enumerate(TARGETS):
            m=models[j]
            for name, c, active_term in zip(["intercept", *FEATURES],m["coef"],m["support"]):
                coeff_rows.append({"generator":f"G{gen}","target":target,"term":name,"coefficient_standardized":float(c),"active":bool(active_term),"feature_mean":None if name=="intercept" else float(m["mu"][FEATURES.index(name)]),"feature_scale":None if name=="intercept" else float(m["sd"][FEATURES.index(name)]),"target_median":m["ymu"],"target_robust_scale":m["ysd"]})
        # Voltage baseline series saved as explicitly separate from SINDy residuals.
        pd.DataFrame({"t_center_s":ph_f.t_center_s,"negative_positive_voltage_ratio":fr,"threshold":vth,"alarm":event_candidate}).to_csv(OUT/f"G{gen}_voltage_ratio_baseline.csv",index=False)
        fig, ax=plt.subplots(figsize=(10,4)); ax.plot(h.t[cal],sh[cal.to_numpy()],label="healthy calibration 8–10 s"); ax.plot(f.t[measured_window],sf[measured_window.to_numpy()],label="fault event"); ax.axhline(threshold,color="black",ls="--",label="healthy q99 threshold"); ax.set(xlabel="Time (s)",ylabel="dimensionless residual score",title=f"G{gen} healthy-model residual detector"); ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(OUT/f"G{gen}_residual_detection.png",dpi=160); plt.close(fig)
        record[gen]={"healthy_rows_raw":len(hraw),"fault_rows_raw":len(fraw),"healthy_timestamp_collision_rows":hraw.attrs["timestamp_collision_rows"],"fault_timestamp_collision_rows":fraw.attrs["timestamp_collision_rows"],"healthy_wave_rows":len(hw),"healthy_slow_rows":len(hs),"fault_wave_rows":len(fw),"fault_slow_rows":len(fs),"healthy_slow_dt_s":float(np.median(np.diff(hs.t))),"fault_slow_dt_s":float(np.median(np.diff(fs.t))),"healthy_missing_by_channel":{c:int(hraw[c].isna().sum()) for c in hraw.columns if c != "t"},"fault_missing_by_channel":{c:int(fraw[c].isna().sum()) for c in fraw.columns if c != "t"}}
        event_evidence.append({"generator":f"G{gen}","sequence_ratio_threshold":vth,"event_start_s":measured_start,"event_end_s":measured_end,"threshold_exceedance_samples":int(active.sum())})
    pd.DataFrame(coeff_rows).to_csv(OUT/"model_coefficients.csv",index=False)
    pd.DataFrame(all_metrics).to_csv(OUT/"generator_metrics.csv",index=False)
    pd.DataFrame(event_evidence).to_csv(OUT/"event_summary.csv",index=False)
    (OUT/"input_audit.json").write_text(json.dumps(record,indent=2))
    write_report(all_metrics,event_evidence,record,coeff_rows)
    print(pd.DataFrame(all_metrics).to_string(index=False))
    print(f"\nSaved experiment outputs to {OUT}")


def write_report(metrics, events, audit, coeff_rows):
    m=pd.DataFrame(metrics); e=pd.DataFrame(events)
    lines=["# Generator-level robust STLSQ fault-detection pilot","",
           "This is a single-event, generator-terminal, binary detection pilot. It is not an implementation of MBAM or constrained SINDy and does not modify the IEEE 9-bus DAE or physical parameters.","",
           "## Method","",
           "Each generator has its own healthy-only model. The source CSV combines waveform and slow streams; the existing `split_merged_streams` function is reused. Waveforms remain at 1,440 Hz and slow channels remain at native 60 Hz. No slow-channel interpolation is performed. Required channels and missing values are audited; same-time fast/slow row collisions in the merged source are deliberately resolved by the shared splitter, then each separated stream is checked for repeated timestamps and its sampling interval. Phasor endpoint windows and one-cycle margins are excluded.","",
           "Features are `[omega/w0-1, theta_rel, Pm-Pm_fit_median, Efd-Efd_fit_median]`; targets are `d(theta_rad)/dt` and `d(omega/w0)/dt`. Theta is interpreted as degrees from its observed range, converted to radians and unwrapped. Derivatives use a centered 7-point, quadratic Savitzky-Golay filter on native slow samples (100 ms span), which smooths noise but blurs rapid fault edges. Each target has an independent robust STLSQ model: standardized columns, Huber IRLS (delta 1.5), ridge 1e-4, sequential threshold 0.08, 15 support updates maximum.","",
           "The dimensionless score is `sqrt(mean_j((y_dot_j - yhat_dot_j)^2 / s_j^2))`, where `s_j` is the robust MAD scale of the corresponding healthy 5–8 s derivative target. For each generator, threshold is the empirical higher-method 99th percentile of healthy calibration scores at 8–10 s, after excluding the 50 ms centered-derivative support at the fit/calibration boundary and the final 1/60 s endpoint. Alarm requires two consecutive 60 Hz scores above threshold, approximately one cycle. No fault samples enter fit, feature scaling, or threshold selection. Fault evaluation excludes 50 ms at each measured event boundary to exceed the requested one-cycle guard and ensure the centered derivative stencil stays inside the event. Derivatives are centered and use three future slow samples, so the score is for offline diagnosis and is not a causal real-time detector as implemented.","",
           "The comparison baseline is `|V2|/|V1|` from one-cycle phasors of the generator-terminal waveforms. Its threshold is independently calibrated by the same 99th-percentile rule on that generator's healthy 8–10 s waveform. This baseline is not a SINDy result.","",
           "## Input audit","",
           "| Generator | Healthy raw/wave/slow rows | Fault raw/wave/slow rows | Slow dt (healthy/fault) | Fault slow missing channels |","|---|---:|---:|---:|---|"]
    for g,a in audit.items(): lines.append(f"| G{g} | {a['healthy_rows_raw']}/{a['healthy_wave_rows']}/{a['healthy_slow_rows']} | {a['fault_rows_raw']}/{a['fault_wave_rows']}/{a['fault_slow_rows']} | {a['healthy_slow_dt_s']:.8f}/{a['fault_slow_dt_s']:.8f} s | healthy {a['healthy_missing_by_channel']}; fault {a['fault_missing_by_channel']} |")
    lines += ["","## Event boundary evidence","","The measured event interval is estimated from negative-to-positive sequence voltage ratio exceedance against the independent healthy 8–10 s q99 threshold; each one-cycle phasor is centered and has one-cycle support. Event evaluation uses a 50 ms guard at each boundary, exceeding one cycle and covering the centered derivative stencil.","","| Generator | Threshold | Measured start (s) | Measured end (s) | Exceedance phasors |","|---|---:|---:|---:|---:|"]
    for _,r in e.iterrows(): lines.append(f"| {r.generator} | {r.sequence_ratio_threshold:.6g} | {r.event_start_s:.6f} | {r.event_end_s:.6f} | {r.threshold_exceedance_samples} |")
    lines += ["","## Results","","| Generator | q99 threshold | Healthy FA rate | Healthy score median / p95 | Fault score median / p95 | Fault alarm fraction | Detected | Delay (s) | V2/V1 fault median | V2/V1 alarm fraction |","|---|---:|---:|---:|---:|---:|---|---:|---:|---:|"]
    for _,r in m.iterrows(): lines.append(f"| {r.generator} | {r.residual_threshold:.5g} | {r.healthy_false_alarm_rate:.3%} | {r.healthy_calibration_median:.4g} / {r.healthy_calibration_p95:.4g} | {r.fault_residual_median:.4g} / {r.fault_residual_p95:.4g} | {r.fault_alarm_fraction:.3%} | {r.detected} | {r.detection_delay_s:.6g} | {r.voltage_ratio_fault_median:.4g} | {r.voltage_ratio_alarm_fraction_event:.3%} |")
    lines += ["","## Active standardized coefficients","","The table lists the fitted coefficients in standardized feature/target coordinates; the intercept is included. Inactive library terms are zero. `model_coefficients.csv` also stores feature means/scales and target normalization needed to reproduce predictions.","","| Generator | Target | Active terms (standardized coefficient) |","|---|---|---|"]
    cdf=pd.DataFrame(coeff_rows)
    for (gen,target), grp in cdf[cdf.active].groupby(["generator","target"], sort=True):
        terms=", ".join(f"{r.term}={r.coefficient_standardized:.6g}" for r in grp.itertuples())
        lines.append(f"| {gen} | {target} | {terms} |")
    lines += ["","Per-sample predictions/residuals are in `G*_healthy_predictions.csv` and `G*_fault_predictions.csv`. Separate voltage baseline series are `G*_voltage_ratio_baseline.csv`; plots are `G*_residual_detection.png`. `generator_metrics.csv`, `event_summary.csv`, and `input_audit.json` provide machine-readable summaries. Run from the workspace root with `python sindy_in_PSA/generator_fault_experiment.py`.","",
              "## Integrity and limitations","",
              "The fault rows are used only for forward prediction and event scoring. There is one healthy run and one fault run per generator; the three generator records are correlated observations from one system event, not independent fault cases. Healthy calibration is a temporal segment of the same run as training. The event interval is operationally inferred from fault-run negative-sequence ratio crossing its healthy-calibrated threshold; therefore the voltage-ratio baseline confirms this same electrical event signature and is not independent ground truth. Boundary timing has one-cycle phasor uncertainty. Savitzky-Golay smoothing blurs onset and low-rate slow channels limit temporal resolution. A successful alarm therefore supports only a single-event offline pilot claim; it does not establish general real-time fault detection across fault locations, operating points, or independent simulations. The voltage baseline is terminal-only and not bus localization.","",
              "No physical parameters or DAE code are changed. No MBAM or constrained SINDy is implemented."]
    (OUT/"README.md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__": main()
