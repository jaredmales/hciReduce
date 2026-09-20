# Step 5: rectangular-PSD precision regularization setup

## Fixed question

Does suppressing the low-eigenvalue modes of the current rectangular-PSD
covariance improve small-separation injection recovery while the independently
measured response remains fixed?

The response-smoothing comparison found gains mainly at radii 6 and 8 pixels.
The measured positive-minus-baseline amplitude increments nevertheless followed
the response-overlap prediction within a few percent. This test therefore keeps
the unsmoothed measured response and changes only how the existing PSD
covariance is inverted.

## Fixed precision policies

At every candidate pixel, the unchanged raw rectangular model uses ±20-pixel
radial pooling, ensemble-mean subtraction, a rectangular 11-by-11 periodogram
on a 21-by-21 padded grid, and the existing 0.3 flat-spectrum mixture. Write its
finite 121-by-121 covariance as

\[
C = Q\,\operatorname{diag}(\nu_i)Q^T,
\qquad
\bar v = \frac{1}{121}\sum_i \nu_i.
\]

The fixed controls and candidates are:

- **full inverse:** \(\pi_i=1/\nu_i\);
- **clipped inverse:**
  \(\pi_i=1/\max(\nu_i,c\bar v)\);
- **hard-truncated inverse:** \(\pi_i=1/\nu_i\) when
  \(\nu_i\geq c\bar v\), and \(\pi_i=0\) otherwise.

The cutoff grid is fixed at **\(c=0.5,0.75,1.0\)**. Clipping retains all 121
modes but limits their maximum precision. Hard truncation removes modes below
the cutoff. The existing 0.3 spectral mixture remains part of \(C\); this test
does not replace or retune it.

For every policy, the amplitude weight is normalized to unit response to the
same unsmoothed measured template. The report records the number of modified or
retained modes and the expected SNR efficiency relative to the full inverse if
the fitted PSD covariance were exact. That efficiency cannot exceed one by the
matched-filter optimum, so any measured recovery gain indicates covariance
model error, calibration effects, or finite-sample behavior rather than an
improvement under an exact covariance.

## Fixed data and controls

The test reuses the completed 164 baseline and 108 positive response-smoothing
analyses and performs no new P4 reduction. It retains the same saved images,
known-planet and one-lambda/D trial exclusions, covariance samples, fitted
means, candidate pixels, five-pixel searches, and production `hciAnalyze`
annular SNR.

Every new precision policy uses the frozen raw rectangular ±20 calibration
pool and receives its own maximum-null threshold before positive measurements
are read. These parent results are copied as references and must reproduce
their amplitude maps, SNR maps, thresholds, and decisions:

- full rectangular ±20 inverse with the unsmoothed response;
- identity with the unsmoothed response;
- identity with the 1.8-pixel response low pass;
- rectangular ±20 with the 2.7-pixel response low pass;
- production Gaussian FWHM 3.6.

The runner also recomputes the full rectangular inverse from the eigensystem at
every fitted pixel and requires agreement with the copied Cholesky-solve parent
map. This checks that any candidate difference comes only from the declared
eigenvalue policy.

## Outputs and interpretation

The final report includes recovery and held-out null counts, mean five-pixel
maximum and fixed-center SNR, paired amplitude throughput, paired decision
changes from the full inverse, exact-PSD efficiency, and retained/modified mode
counts. Results at radii 6 and 8 are the main diagnostic, with radii 12–24 kept
as controls. The same correlated field and previously inspected injections
make this a development comparison; it cannot select a radius-dependent
production policy by itself.

## Validation and ROC commands

After pulling the setup, run the local algebra check:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/compare_p4_step5_psd_precision.py check
```

Prepare the immutable comparison:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
  /opt/conda/envs/xpy3_13/bin/python3 \
  agents/plans/scripts/compare_p4_step5_psd_precision.py prepare \
  --comparison working/roc/p4_response_smoothing_20260920 \
  --root working/roc/p4_psd_precision_20260920 \
  --cpus 0 1 2 3 4 5 6 7 8 9 10 11
```

Launch the resumable analysis:

```sh
tmux new-session -d -s p4-psd-precision \
  "cd /home/jrmales/Source/mxApps/hciReduce && \
   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
   MPLCONFIGDIR=/tmp/p4-step5-matplotlib \
   /opt/conda/envs/xpy3_13/bin/python3 \
   agents/plans/scripts/compare_p4_step5_psd_precision.py run \
   --root working/roc/p4_psd_precision_20260920 \
   > working/roc/p4_psd_precision_20260920/driver.log 2>&1"
```

Monitor it with:

```sh
cat working/roc/p4_psd_precision_20260920/state.json
tail -n 30 working/roc/p4_psd_precision_20260920/driver.log
```

The final products will be `results.json`, `results.md`, and `comparison.png`
under the new root.
