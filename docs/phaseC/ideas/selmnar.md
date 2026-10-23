# selmnar: Heckman probit selection model inside an EM U-step of the adapted projection

- **Queen:** representation & mathematics (round-3 directive)
- **Status:** building (this section was written BEFORE any screening run; it is the preregistration of the variants)

## Preregistration (written before the first screen, 2026-10-23)
- `impute`: the spec's STEPS 0-3, 4 EM iterations, tilt factor rho = 1 (full Heckman posterior mean).
- `impute_null`: same machinery; the selection fit sees the missing indicators permuted within each cell among ok
  proteins (`default_rng(seed + 4242)`), with the marginal probit term for every entry; E-step on the TRUE missing
  entries.
- `impute_inf`: sigma_d fixed at 1e6 (no selection) -> machinery-only control; required max |inf - champion| / sd_j
  < 1e-3 on every fast target; harness cand_raw expected 1.000 +/- 0.002.
- `impute_cal` (preregistered variant): rho in {0, 0.25, 0.5, 1} multiplies the tilt in every E-step,
  xhat = m - rho (v/s) lam. rho is chosen inside `impute` by self-simulation: using the fitted selection curve of the
  rho = 1 fit on the input, hide 10% of the observed eligible entries (eligible = ok proteins with observed fraction
  >= 0.5 in the input) with probability proportional to Phi((c_ij - x_ij)/sigma_d) (rng `default_rng(seed + 1)`),
  rerun STEPS 1-3 on that input with mu, sd, V and s2 fixed (U restarted by one ridge U-step on the observed entries),
  and pick the rho with the lowest RMSE on the hidden entries. Ties -> smaller rho. Screened as `selmnar_cal`
  (impute_cal vs impute_null).
- Promotion test (fixed): README criterion (fast cand_stack <= 0.99 and >= 0.005 below null_stack). Secondary readout
  (not a pass): MNAR-mean cand_raw <= 0.98 and >= 0.01 below the null's MNAR cand_raw (`selmnar_nullraw`, cand =
  impute_null). Credit only the margin below both null_stack and selmnar_inf.
- **Revision v2 (added after the v1/cal/inf/nullraw screens, BEFORE running v2).** v1 tilted every missing entry,
  including random dropouts, because the selection model had no MAR component (fast MCAR cand_raw 1.2721). v2 adds a
  global random-dropout probability pi: P(miss | x) = pi + (1 - pi) Phi((c - x)/sigma_d), fitted jointly (logit pi in
  [-9, 2.2], init 0.1). Closed-form E-step: xhat = m - (v/s) (1-pi) phi(alpha) / (pi + (1-pi) Phi(alpha)); rho = 1;
  everything else as v1. Null `impute_null_v2`: the same within-cell permutation null with the mixture marginal
  likelihood. Screened as `selmnar_v2` and `selmnar_v2_nullraw` (cand = impute_null_v2, for the secondary readout).
