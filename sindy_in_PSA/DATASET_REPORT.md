# Prepared healthy segment: 5–10 seconds

The user designated samples from 5 seconds onward as healthy. The earlier A–G disturbance near 2.0–2.3 seconds and its recovery are excluded. The 5–8 second portion is the fitting segment; 8–10 seconds is a temporal holdout from the same run, not an independent simulation validation set.

## Stream separation and signal handling

The source CSVs combine fast phase waveforms and slower channels. Their rows include near-coincident timestamps where the slow stream was joined to the waveform stream. Preparation groups timestamps after rounding to 9 decimal places. At collisions, the waveform-only row is retained for the fast stream and the row with recorded slow values is retained for the slow stream. Singleton rows are retained for both streams. No interpolation is applied to slow channels. Raw source CSVs are left unchanged.

The phase waveform stream is retained at its native 1,440 Hz grid (24 samples per nominal 60 Hz cycle). A separate one-cycle least-squares RMS phasor series is produced at 60 samples/s, including phase phasors and zero-, positive-, and negative-sequence components. Raw exported units are retained; per-unit normalization is deferred until all machine and bus bases and current sign conventions are verified.

## Files and split

For each asset, `*_train.csv` contains 5.0 <= t < 8.0 s and `*_validation.csv` contains 8.0 <= t <= 10.0 s. Files with `waveforms_raw_units` retain instantaneous phase channels; generator `slow_channels_raw_units` retain generator channels; bus `PQ_raw_units` retain the recorded P/Q values; `sequence_phasors_raw_units` contain one-cycle phasors and sequence components.

| Source | waveform rows kept | slower rows kept | time span (s) |
| --- | ---: | ---: | ---: |
| `G1_Fag_10s.csv` | 14401 | 601 | 0.000000–10.000000 |
| `G2_Fag_10s.csv` | 14401 | 601 | 0.000000–10.000000 |
| `G3_Fag_10s.csv` | 14401 | 601 | 0.000000–10.000000 |
| `IVPQ_B2_H.csv` | 14401 | 601 | 0.000000–10.000000 |
| `IVPQ_B7_H.csv` | 14401 | 601 | 0.000000–10.000000 |
| `IVPQ_B8_H.csv` | 14401 | 601 | 0.000000–10.000000 |

## Limitations and decisions still open

- This split uses one run only. The temporal holdout does not establish generalization to an independent healthy simulation or operating point.
- Source filenames and folder placement retain mixed healthy/fault naming. The user's explicit instruction defines 5–10 s as healthy for this preparation; verify the post-fault steady state before using it as the sole healthy reference.
- The source's slower P/Q channels have missing values by design on the fast time grid; they are separated into their own 60 Hz files rather than filled.
- The generator ratings, bus voltage bases, pole-pair conversions, and current sign convention in `mbam_package/prep10.py` are not independently verified. Prepared phasors are therefore saved in raw export units, with no guessed per-unit scaling.
- The original input files contain 12 exact timestamp collisions and many additional slow/fast timestamp collisions at floating-point precision. The separated streams retain one record per rounded timestamp per stream; inspect the source export if sub-nanosecond timestamp identity is material.

See `prepared_data/healthy_5s_onward/healthy_5s_onward_overview.png` for generator speed and bus positive-sequence voltage over the selected segment.
