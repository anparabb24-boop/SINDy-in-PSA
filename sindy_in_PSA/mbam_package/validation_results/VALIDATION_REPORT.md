# Physical Generator Model Validation (Step 2)

The full two-axis/subtransient physical model was fit on 5.0 <= t < 8.0 s and evaluated on 8.0 <= t <= 10.0 s. Only the fitting segment contributes to the optimizer residuals and robust residual scales. The validation segment is simulated continuously from the fitted model initial state without parameter fitting.

This is a same-run temporal validation. It does not establish generalization to an independent healthy simulation. The P/Q targets are derived from the terminal positive-sequence phasors, not independently measured for all three generator terminals. G2's derived power is cross-checked against the logged B2 P/Q channels.

Per-unit bases used here are the nominal generator ratings in the package (247.5, 230, and 128 MVA) and nominal terminal voltage bases (16.5, 18, and 13.8 kV line-to-line). Speed is divided by nominal synchronous mechanical speed using pole-pair counts (20, 1, 1); rotor angle in degrees is converted to electrical radians using the same counts. These mappings are provisional pending confirmation against the source Simulink blocks.

A provisional trajectory acceptance criterion for this run is validation NRMSE <= 0.15 for each output, normalized by that output's training range. Passing this criterion supports predictive fit on the same-run holdout only; it does not establish parameter identifiability or independent-run generalization.

Current result: all outputs pass the provisional criterion; maximum validation NRMSE is 0.14357 (G3 active power). The optimizer converged for G2; G1 and G3 parameters remain provisional because those fits hit their evaluation limits.

## Fit status

```json
{
  "G1": {
    "success": false,
    "message": "The maximum number of function evaluations is exceeded.",
    "fit_cost": 7830.770245053511,
    "fit_nfev": 400,
    "stage1_success": true,
    "stage1_nfev": 23,
    "stage1_cost": 64769.98826621657,
    "training_samples": 180,
    "validation_samples": 121
  },
  "G2": {
    "success": true,
    "message": "`ftol` termination condition is satisfied.",
    "fit_cost": 3431.891789232716,
    "fit_nfev": 392,
    "stage1_success": true,
    "stage1_nfev": 17,
    "stage1_cost": 20713.557905389443,
    "training_samples": 180,
    "validation_samples": 121
  },
  "G3": {
    "success": false,
    "message": "The maximum number of function evaluations is exceeded.",
    "fit_cost": 3520.3096417025668,
    "fit_nfev": 400,
    "stage1_success": true,
    "stage1_nfev": 16,
    "stage1_cost": 22533.58224423095,
    "training_samples": 180,
    "validation_samples": 121
  }
}
```

## Metrics

`metrics.csv` contains RMSE and NRMSE. NRMSE is normalized by the training-segment range for each output.

```text
generator                output      split   n       rmse  normalization_scale_train_range     nrmse
       G1                  P_pu      train 180  0.0025681                          0.17023  0.015086
       G1                  P_pu validation 121   0.019551                          0.17023   0.11485
       G1                  Q_pu      train 180 0.00088287                         0.024153  0.036553
       G1                  Q_pu validation 121  0.0013735                         0.024153  0.056866
       G1                  w_pu      train 180 9.9656e-06                         0.001001  0.009956
       G1                  w_pu validation 121 1.2137e-05                         0.001001  0.012126
       G1 delta_minus_theta_rad      train 180  0.0006458                          0.02853  0.022636
       G1 delta_minus_theta_rad validation 121 0.00023382                          0.02853 0.0081955
       G2                  P_pu      train 180  0.0015022                          0.13748  0.010927
       G2                  P_pu validation 121  0.0074554                          0.13748  0.054228
       G2                  Q_pu      train 180 0.00042345                         0.025463   0.01663
       G2                  Q_pu validation 121  0.0020849                         0.025463  0.081882
       G2                  w_pu      train 180 4.7404e-05                        0.0029625  0.016001
       G2                  w_pu validation 121 1.5657e-05                        0.0029625 0.0052852
       G2 delta_minus_theta_rad      train 180 0.00057548                         0.033154  0.017358
       G2 delta_minus_theta_rad validation 121 0.00042457                         0.033154  0.012806
       G3                  P_pu      train 180  0.0011992                         0.072263  0.016595
       G3                  P_pu validation 121   0.010375                         0.072263   0.14357
       G3                  Q_pu      train 180  0.0002178                         0.016197  0.013447
       G3                  Q_pu validation 121  0.0019414                         0.016197   0.11986
       G3                  w_pu      train 180 2.9161e-05                        0.0019858  0.014685
       G3                  w_pu validation 121 2.7014e-05                        0.0019858  0.013604
       G3 delta_minus_theta_rad      train 180 0.00040105                         0.022928  0.017492
       G3 delta_minus_theta_rad validation 121 0.00035079                         0.022928    0.0153
```

## Limitations

- Bus B2 P/Q independently supports the selected export current sign for G2; equivalent logged terminal P/Q cross-checks are unavailable for G1 and G3 in the current data bundle.
- Excitation and mechanical-power channel units/scales, nominal machine parameters, and generator voltage/power bases still require confirmation from the original model/export configuration.
- Do not run MBAM reduction or report a validated reduced model unless the full-model baseline is acceptable on the held-out interval.

The fit/prediction trajectories are in `healthy_generator_model_validation.png`; parameter estimates are in `fitted_parameters.csv`.
