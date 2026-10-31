"""Separate and prepare the user-designated healthy segment (t >= 5 s).

The source CSVs contain a fast waveform stream and a slower stream merged into
one table. At slower-channel timestamps, those streams can produce two nearly
identical timestamps. This script separates the streams before sorting or
splitting so it does not discard waveform samples or interpolate slow states.

Run from any working directory with:
    python sindy_in_PSA/prepare_healthy_data.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
SOURCE_DIR = ROOT / "fault datasets"
OUTPUT_DIR = ROOT / "prepared_data" / "healthy_5s_onward"
HEALTHY_START = 5.0
TRAIN_END = 8.0
FAULT_FREQUENCY_HZ = 60.0
TIME_KEY_DECIMALS = 9
SAMPLES_PER_CYCLE = 24

GEN_FILES = {
    1: ("G1_Fag_10s.csv", ("Pm_pu", "theta", "w", "Excitation V")),
    2: ("G2_Fag_10s.csv", ("Pmpu", "theta", "w", "excitation V")),
    3: ("G3_Fag_10s.csv", ("Pmpu", "theta", "w", "excitation V")),
}
BUS_FILES = {k: f"IVPQ_B{k}_H.csv" for k in (2, 7, 8)}
WAVE_COLS = ("Va", "Vb", "Vc", "Ia", "Ib", "Ic")
PHASE_COLS = {
    "voltage": ("Va", "Vb", "Vc"),
    "current": ("Ia", "Ib", "Ic"),
}


def split_merged_streams(path: Path, slow_cols: tuple[str, ...]):
    """Return a uniform fast waveform stream and an un-interpolated slow stream."""
    df = pd.read_csv(path)
    needed = {"t", *WAVE_COLS, *slow_cols}
    missing = sorted(needed.difference(df.columns))
    if missing:
        raise ValueError(f"{path.name}: missing required columns: {missing}")
    for col in needed:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["t", *WAVE_COLS]).copy()
    df["_time_key"] = df["t"].round(TIME_KEY_DECIMALS)

    # Fast records are the rows without a slow measurement when a timestamp
    # collision occurs. A singleton row belongs to both streams.
    group_sizes = df.groupby("_time_key")["t"].transform("size")
    has_slow = df[list(slow_cols)].notna().any(axis=1)
    fast = df.loc[(group_sizes == 1) | ~has_slow, ["t", *WAVE_COLS]].copy()
    fast["_time_key"] = fast["t"].round(TIME_KEY_DECIMALS)
    if fast["_time_key"].duplicated().any():
        collisions = fast.loc[fast["_time_key"].duplicated(False), "_time_key"].unique()
        raise ValueError(
            f"{path.name}: ambiguous waveform rows at rounded times {collisions[:5]}"
        )

    slow = df.dropna(subset=list(slow_cols)).loc[:, ["t", *slow_cols]].copy()
    slow["_time_key"] = slow["t"].round(TIME_KEY_DECIMALS)
    if slow["_time_key"].duplicated().any():
        raise ValueError(f"{path.name}: repeated slow-channel timestamps")

    fast = fast.sort_values("t").reset_index(drop=True)
    slow = slow.sort_values("t").reset_index(drop=True)
    fast = fast.drop(columns="_time_key")
    slow = slow.drop(columns="_time_key")
    return fast, slow


def positive_sequence_phasors(wave: pd.DataFrame, channel_group: str):
    """One-cycle least-squares RMS phasors and Fortescue sequence components."""
    phases = PHASE_COLS[channel_group]
    t = wave["t"].to_numpy(float)
    dt = np.diff(t)
    dt0 = float(np.median(dt))
    expected = 1.0 / (FAULT_FREQUENCY_HZ * SAMPLES_PER_CYCLE)
    if not np.allclose(dt, expected, rtol=1e-6, atol=1e-10):
        raise ValueError(
            f"Expected a uniform {FAULT_FREQUENCY_HZ*SAMPLES_PER_CYCLE:g} Hz "
            f"waveform grid; {channel_group} median dt={dt0:.12g} s"
        )
    n = SAMPLES_PER_CYCLE
    omega = 2.0 * np.pi * FAULT_FREQUENCY_HZ
    # One complete cycle, non-overlapping windows, centered timestamps.
    starts = np.arange(0, len(wave) - n + 1, n)
    center_times = t[starts] + n / (2.0 * FAULT_FREQUENCY_HZ * SAMPLES_PER_CYCLE)
    fit_t = t[:n]
    design = np.column_stack(
        [np.ones(n), np.cos(omega * (fit_t - fit_t[0])), np.sin(omega * (fit_t - fit_t[0]))]
    )
    phasors = {}
    inverse = np.linalg.pinv(design)
    for phase in phases:
        y = wave[phase].to_numpy(float)
        windows = sliding_window_view(y, n)[starts]
        coef = windows @ inverse.T
        # RMS phasor with local cycle start as reference; report the
        # corresponding absolute-time phase by rotating the coefficient.
        local = (coef[:, 1] - 1j * coef[:, 2]) / np.sqrt(2.0)
        phasors[phase] = local * np.exp(-1j * omega * t[starts])

    a = np.exp(2j * np.pi / 3.0)
    A, B, C = (phasors[p] for p in phases)
    seq = {
        "zero": (A + B + C) / 3.0,
        "positive": (A + a * B + a**2 * C) / 3.0,
        "negative": (A + a**2 * B + a * C) / 3.0,
    }
    result = {"t_center_s": center_times}
    for phase, values in phasors.items():
        result[f"{phase}_re"] = values.real
        result[f"{phase}_im"] = values.imag
        result[f"{phase}_magnitude"] = np.abs(values)
    for name, values in seq.items():
        result[f"{channel_group}_{name}_re"] = values.real
        result[f"{channel_group}_{name}_im"] = values.imag
        result[f"{channel_group}_{name}_magnitude"] = np.abs(values)
    return pd.DataFrame(result)


def save_splits(frame: pd.DataFrame, stem: str, time_col: str):
    healthy = frame.loc[frame[time_col] >= HEALTHY_START].copy()
    train = healthy.loc[healthy[time_col] < TRAIN_END]
    validation = healthy.loc[healthy[time_col] >= TRAIN_END]
    train.to_csv(OUTPUT_DIR / f"{stem}_train.csv", index=False)
    validation.to_csv(OUTPUT_DIR / f"{stem}_validation.csv", index=False)
    return train, validation


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    waveforms = {}
    slow_data = {}
    phasors = {}
    summary = []

    for gen, (filename, slow_cols) in GEN_FILES.items():
        fast, slow = split_merged_streams(SOURCE_DIR / filename, slow_cols)
        fast.to_csv(OUTPUT_DIR / f"G{gen}_waveforms_raw_units.csv", index=False)
        slow.to_csv(OUTPUT_DIR / f"G{gen}_slow_channels_raw_units.csv", index=False)
        save_splits(fast, f"G{gen}_waveforms", "t")
        save_splits(slow, f"G{gen}_slow_channels", "t")
        waveforms[f"G{gen}"] = fast
        slow_data[f"G{gen}"] = slow
        summary.append((filename, len(fast), len(slow), fast.t.min(), fast.t.max()))

        gphasors = positive_sequence_phasors(fast, "voltage")
        cphasors = positive_sequence_phasors(fast, "current")
        combined = gphasors.merge(cphasors, on="t_center_s", validate="one_to_one")
        combined.to_csv(OUTPUT_DIR / f"G{gen}_sequence_phasors_raw_units.csv", index=False)
        save_splits(combined, f"G{gen}_sequence_phasors", "t_center_s")
        phasors[f"G{gen}"] = combined

    for bus, filename in BUS_FILES.items():
        fast, slow = split_merged_streams(SOURCE_DIR / filename, ("P", "Q"))
        fast.to_csv(OUTPUT_DIR / f"B{bus}_waveforms_raw_units.csv", index=False)
        slow.to_csv(OUTPUT_DIR / f"B{bus}_PQ_raw_units.csv", index=False)
        save_splits(fast, f"B{bus}_waveforms", "t")
        save_splits(slow, f"B{bus}_PQ", "t")
        waveforms[f"B{bus}"] = fast
        summary.append((filename, len(fast), len(slow), fast.t.min(), fast.t.max()))

        vphasors = positive_sequence_phasors(fast, "voltage")
        iphasors = positive_sequence_phasors(fast, "current")
        combined = vphasors.merge(iphasors, on="t_center_s", validate="one_to_one")
        combined.to_csv(OUTPUT_DIR / f"B{bus}_sequence_phasors_raw_units.csv", index=False)
        save_splits(combined, f"B{bus}_sequence_phasors", "t_center_s")
        phasors[f"B{bus}"] = combined

    # Overview: slow generator speed plus positive-sequence phase voltage.
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, constrained_layout=True)
    for gen, (_, slow_cols) in GEN_FILES.items():
        s = slow_data[f"G{gen}"]
        w_col = next(c for c in slow_cols if c.lower() == "w")
        s = s[(s.t >= HEALTHY_START)]
        axes[0].plot(s.t, s[w_col] / float(s[w_col].iloc[0]), label=f"G{gen}")
    axes[0].set_ylabel("Speed / initial speed")
    axes[0].set_title("User-designated healthy segment (5–10 s)")
    axes[0].legend(ncol=3)
    for asset in ("B2", "B7", "B8"):
        p = phasors[asset]
        p = p[p.t_center_s >= HEALTHY_START]
        axes[1].plot(
            p.t_center_s,
            p["voltage_positive_magnitude"],
            label=f"{asset} positive-sequence |V1|",
        )
    axes[1].set_ylabel("Positive-sequence voltage (raw export units)")
    axes[1].set_xlabel("Time (s)")
    axes[1].legend(ncol=3)
    for ax in axes:
        ax.axvline(TRAIN_END, color="black", linestyle="--", linewidth=1)
        ax.grid(alpha=0.25)
    fig.savefig(OUTPUT_DIR / "healthy_5s_onward_overview.png", dpi=160)
    plt.close(fig)

    lines = [
        "# Prepared healthy segment: 5–10 seconds",
        "",
        "The user designated samples from 5 seconds onward as healthy. The earlier A–G disturbance near 2.0–2.3 seconds and its recovery are excluded. The 5–8 second portion is the fitting segment; 8–10 seconds is a temporal holdout from the same run, not an independent simulation validation set.",
        "",
        "## Stream separation and signal handling",
        "",
        "The source CSVs combine fast phase waveforms and slower channels. Their rows include near-coincident timestamps where the slow stream was joined to the waveform stream. Preparation groups timestamps after rounding to 9 decimal places. At collisions, the waveform-only row is retained for the fast stream and the row with recorded slow values is retained for the slow stream. Singleton rows are retained for both streams. No interpolation is applied to slow channels. Raw source CSVs are left unchanged.",
        "",
        "The phase waveform stream is retained at its native 1,440 Hz grid (24 samples per nominal 60 Hz cycle). A separate one-cycle least-squares RMS phasor series is produced at 60 samples/s, including phase phasors and zero-, positive-, and negative-sequence components. Raw exported units are retained; per-unit normalization is deferred until all machine and bus bases and current sign conventions are verified.",
        "",
        "## Files and split",
        "",
        "For each asset, `*_train.csv` contains 5.0 <= t < 8.0 s and `*_validation.csv` contains 8.0 <= t <= 10.0 s. Files with `waveforms_raw_units` retain instantaneous phase channels; generator `slow_channels_raw_units` retain generator channels; bus `PQ_raw_units` retain the recorded P/Q values; `sequence_phasors_raw_units` contain one-cycle phasors and sequence components.",
        "",
        "| Source | waveform rows kept | slower rows kept | time span (s) |",
        "| --- | ---: | ---: | ---: |",
    ]
    for filename, nfast, nslow, tmin, tmax in summary:
        lines.append(f"| `{filename}` | {nfast} | {nslow} | {tmin:.6f}–{tmax:.6f} |")
    lines += [
        "",
        "## Limitations and decisions still open",
        "",
        "- This split uses one run only. The temporal holdout does not establish generalization to an independent healthy simulation or operating point.",
        "- Source filenames and folder placement retain mixed healthy/fault naming. The user's explicit instruction defines 5–10 s as healthy for this preparation; verify the post-fault steady state before using it as the sole healthy reference.",
        "- The source's slower P/Q channels have missing values by design on the fast time grid; they are separated into their own 60 Hz files rather than filled.",
        "- The generator ratings, bus voltage bases, pole-pair conversions, and current sign convention in `mbam_package/prep10.py` are not independently verified. Prepared phasors are therefore saved in raw export units, with no guessed per-unit scaling.",
        "- The original input files contain 12 exact timestamp collisions and many additional slow/fast timestamp collisions at floating-point precision. The separated streams retain one record per rounded timestamp per stream; inspect the source export if sub-nanosecond timestamp identity is material.",
        "",
        "See `prepared_data/healthy_5s_onward/healthy_5s_onward_overview.png` for generator speed and bus positive-sequence voltage over the selected segment.",
        "",
    ]
    (ROOT / "DATASET_REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Prepared outputs written to {OUTPUT_DIR}")
    for filename, nfast, nslow, tmin, tmax in summary:
        print(f"{filename}: waveforms={nfast}, slow={nslow}, t={tmin:g}..{tmax:g} s")


if __name__ == "__main__":
    main()
