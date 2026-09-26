# KLIP covariance-aware matched-filter test program

## Purpose

The earlier AF Lep/NACO KLIP experiments established that the measured response
is useful for localization and photometry but is a worse detection filter than
simple Gaussian smoothing when it is correlated under an identity covariance.
For the exact signal-free response, SNR was 0.825--0.902 times the unfiltered
image across the eight KL modes, while Gaussian FWHM 3.6 was 1.040--1.098 times
unfiltered. The exact and sparse response filters had nearly the same SNR, so
sparse radial response interpolation was not the principal cause.

This program tests the remaining explanation: the KLIP response's negative
lobes overlap correlated final-image residuals, and identity weighting gives
those modes the wrong detection weight. It transfers the full set of P4
covariance questions to KLIP while preserving Gaussian FWHM 3.6 as the primary
reference.

The program separates four questions:

1. Is the exact KLIP response accurate for finite positive injections at the
   contrasts used by the detection experiment?
2. Which final-image covariance estimates predict held-out KLIP residuals
   without unstable inverse weights?
3. Does covariance weighting improve injection recovery and mean SNR at a
   threshold fixed from baseline nulls?
4. Do any gains survive a fresh, unopened validation set and remain sensible
   on the known planet?

## Existing products to reuse

The archived response calculation at
`working/roc/klip_signal_free_pixel_response_20260914T232006Z` supplies:

- an optimized KLIP negative companion at separation 12.3877505 pixels,
  position angle 260.643152 degrees, and contrast 0.0045743625;
- a signal-free eight-plane baseline at KL modes
  125, 150, 175, 200, 225, 250, 300, and 350;
- exact 11-by-11 `refitDifference` response and validity stamps at all 11,192
  native pixels in the 6-to-60-pixel search annulus;
- the production sparse radial response for a response-model control; and
- the original planet-bearing science cube for the final descriptive endpoint.

The 11-pixel field remains the replay oracle. Response-convergence and
contrast-linearity tests subsequently promoted a 47-pixel footprint, so a
[new exact-response campaign](../klip-response-47-setup-20260921/README.md)
must pass before covariance screening. Its central 11 pixels must reproduce
this archive at every location and its complete stamps must reproduce 18
independent finite-difference derivatives.

The archived exact response used hciReduce commit `9bfe31e`. Before reuse, a
current no-response signal-free KLIP reduction must reproduce the archived
eight-plane baseline at a fixed numerical tolerance. If this compatibility
gate fails, the injection study must use a build of the archived commit or the
exact response must be regenerated. A response field and positive reductions
from measurably different reduction operators may not be combined.

The compatibility comparison requires identical cube dimensions, mode order,
and finite-pixel masks, plus `np.allclose` agreement with `rtol=2e-6` and
`atol=5e-7` in final-image units. These tolerances cover the previously
observed $3.95\times10^{-7}$ repeat-reduction maximum difference without
permitting a visible change in the residual operator.

## Fixed scientific contract

### Reduction and response

- Use the 621 AF Lep/NACO coadd5 frames and
  `agents/plans/scripts/klipReduce_afLepNaco_psf_response.conf`.
- Preserve the eight archived KL mode counts and their plane ordering.
- Use the validated 47-pixel exact signal-free per-pixel response as the
  primary template.
- Use the sparse candidate-avoiding radial response only as a response-model
  control.
- Retain central 11- and 31-pixel crops as fixed footprint controls. They test
  whether covariance weighting uses the nonlocal response tail or suppresses
  it; they are not substitutes for measuring that tail.
- Mode 200 is the preregistered primary endpoint because the independent
  negative-companion fit was optimized at that mode. Report every method at
  all eight modes, but do not select the best mode after seeing an injection.
- Do not maximize across mode planes unless a distinct joint-mode null
  threshold is calibrated. The mode planes are correlated reductions of the
  same frames.

For the injection reductions, retain the optimized source in `[planet]`, set
`fake.subtractPlanet=true`, and add the trial source through the vector-valued
`fake.sep`, `fake.PA`, and `fake.contrast` controls. This produces the same
signal-free baseline used by the exact response while adding one positive
source. No response is recomputed for a positive image.

The known-planet endpoint uses the original planet-bearing cube and the search
geometry in `working/analyze.conf`: lambda/D 3.6 pixels, minimum SNR radius 6,
planet exclusion radius 7, and aperture radius 3. The optimized KLIP position
defines the negative subtraction used to create signal-free injection images;
it does not replace the analysis aperture in `working/analyze.conf`.

### Detection geometry

Use nominal radii 7.5, 10, 12, 16, 20, and 24 pixels. At each radius, partition
native candidate locations before reading their scores:

- 20 maximum-null calibration locations;
- 6 development-injection locations;
- 6 validation-injection locations; and
- 6 held-out-null locations.

Choose the locations with deterministic angular maximin spacing and require
complete response, candidate, covariance-training, and annular-SNR support for
every method in a comparison. At radius 7.5 this partition uses most of the
available circumference and is necessarily correlated; it must not be
described as 38 independent resolution elements.

A geometry-only preflight of the archived exact manifest found just 12 fully
supported centers in the nominal radius-6 bin and no center whose four axial
neighbors also have complete responses. This is caused by the KLIP inner
processing boundary and the 11-by-11 response footprint. The radius-7.5 bin has
52 complete five-pixel centers and is the smallest feasible primary bin for the
common-support design. Radius 6 remains a separate boundary diagnostic for
response support and center-pixel behavior; it does not enter recovery or
threshold comparisons unless a new partial-support protocol is designed and
validated first.

For the archived 11-pixel field, the same preflight found 42, 44, 60, 91, 113,
and 160 eligible centers in the six primary bins after requiring all five
search centers to have complete responses and to lie outside the seven-pixel
fitted-planet exclusion. The 47-pixel campaign must repeat this geometry-only
count before the site partition is frozen. Site selection uses only coordinates
and validity metadata, never baseline or positive scores. If a bin cannot
support the fixed partition on the promoted footprint, the report must retain
the smaller-footprint result as a diagnostic and revise or drop that bin before
any scores are read.

Each injection is reduced separately. Its five-pixel detection search is the
center plus the four cardinal neighbors, matching the P4 study. Training for
each of those pixels excludes the union of their full response footprints.
Retain the conservative known-planet exclusion even in signal-free images so a
negative-subtraction residual cannot train the covariance.

The three positive contrasts target Gaussian-FWHM-3.6 source-only SNR values
3, 5, and 7 at mode 200. Calculate those contrasts from the baseline annular
noise and the exact processed response before any positive reduction. Freeze
all contrasts, locations, exclusions, methods, and thresholds in the manifest.
The other seven modes receive the same physical contrasts and are secondary
tests of mode dependence.

This creates 36 development sites and 36 validation sites, with three levels
at each site: 216 short KLIP reductions. At the archived 2.7--3.0 seconds per
reduction, the reduction stage should take roughly 10--12 minutes on ROC. It
does not repeat the promoted exact-response calculation after that field has
passed its separate campaign.

### Common estimator and reporting rules

Every response-based method evaluates

$$
\widehat\alpha=\frac{t^T C^{-1}(d-\widehat\mu)}
{t^T C^{-1}t},
$$

using one consistent support and coordinate order for the candidate data,
response, fitted mean, and covariance. Restrict all four objects to valid
support before solving. Report signed amplitude, conditional variance,
annularly calibrated SNR, support, fitted sample count, and solve diagnostics
separately.

For every positive image, also compare the fixed exact response with the
finite difference `(positive - signal-free baseline) / injected contrast`.
Record covariance-weighted and unweighted cosine, projection scale, and
best-scaled residual. This is a response-fidelity diagnostic only. A template
calculated from the positive image may not be used to detect that injection.

All method comparisons use common eligible support. A method that rejects a
difficult location does not receive credit for improved mean SNR or recovery.
The production annular SNR calculation and an independent oracle must agree.

## Method suite

The suite is staged so estimator choices are not multiplied into an
uninterpretable full factorial. Each stage preserves the controls needed to
identify which part changes the result.

### Permanent references

These methods appear in every injection table:

| Reference | Purpose |
| --- | --- |
| Native image, no spatial filter | Shows the starting KLIP statistic. |
| Gaussian FWHM 3.6 | Primary detection reference. |
| Exact response, identity covariance | Reproduces the known matched-filter failure with the best available response. |
| Sparse response, identity covariance | Measures response interpolation separately from covariance weighting. |
| Exact response low-pass filtered at 1.8 pixels | Tests the P4 result that suppressing fine template structure can improve recovery; candidate data remain unsmoothed. |
| Exact response low-pass filtered at 2.7 pixels | Stronger template-regularization control; candidate data remain unsmoothed. |
| Exact response with a fitted mean and isotropic covariance | Separates background subtraction from correlation weighting. |

A Gaussian development-only width sweep at 2.4, 3.0, 3.6, and 4.2 pixels is
reported but cannot replace the preregistered 3.6-pixel reference in validation.

### Stage A: response and geometry audit

Before fitting covariance:

1. verify all input and product hashes, mode ordering, coordinates, validity,
   response normalization, and signal-free baseline compatibility;
2. measure exact-response energy and negative-lobe energy by radius, including
   the fraction on the stamp border;
3. verify exact versus sparse response projection and cosine by radius and
   mode; and
4. reproduce the archived unfiltered, exact-response, sparse-response, and
   Gaussian planet SNR values.

Failure of the baseline compatibility or response-coordinate checks blocks the
injection experiment. Large response edge energy is recorded as a limitation
and triggers a larger-stamp response experiment before claiming an optimal
matched filter. The fixed trigger is more than one percent median squared
response energy on the outermost stamp border in any primary mode-200 radial
bin, or more than five percent for any selected injection location.

The maintained Stage-A runner is
[`run_klip_covariance_stage_a.py`](../../scripts/run_klip_covariance_stage_a.py).
Run both phases under the same CPU affinity so the recorded OpenMP resource
contract is replayed exactly:

```bash
root=working/roc/klip_covariance_stage_a_20260921
taskset -c 12-27 python3 agents/plans/scripts/run_klip_covariance_stage_a.py check
taskset -c 12-27 python3 agents/plans/scripts/run_klip_covariance_stage_a.py prepare "$root"
taskset -c 12-27 python3 "$root/software/run_klip_covariance_stage_a.py" run "$root"
```

If the installed executables are not the intended ROC builds, pass their
paths to `prepare` with `--klipreduce` and `--hcianalyze`. Preparation records
their hashes alongside the hashes from the archived exact-response run. The
runner stops on changed frozen inputs or CPU affinity, preserves failed task
directories for diagnosis, and reuses verified completed tasks on restart.

The canonical ROC run completed on 2026-09-21. Its `klipReduce` hash matches
the archived exact-response run and the freshly reduced signal-free cube is
bitwise identical to the archived baseline. Geometry, the independent exact
identity-filter replay, and all four planet controls pass. The response-edge
trigger fires in every primary mode-200 radial bin, so the next checkpoint is
a fixed-site larger-stamp convergence test before covariance screening chooses
a final template footprint. See the [Stage-A result](../klip-covariance-stage-a-checkpoint-20260921/README.md).

The completed
[convergence experiment](../klip-response-stamp-convergence-setup-20260921/README.md)
used 12 geometry-only sites per primary radius and 144 paired perturbation
reductions. Its external 11-pixel derivative reproduces the archive with
minimum cosine 0.999999 and maximum projection error $4.40\times10^{-5}$.
No candidate through 31 pixels clears the edge rule at every primary mode-200
radius. A 39-pixel extraction still fails at radius 10, while 47 pixels is the
first complete-support size that clears the one-percent median and five-percent
individual thresholds in all six radii and all eight modes.

The newly exposed tail is material at small separation. The square shell from
31 through 47 pixels contains median 4.79%, 7.87%, and 4.53% of the 63-pixel
mode-200 response energy at radii 7.5, 10, and 12, with individual values as
large as 11.18%. Because sizes above 31 were preregistered as diagnostics and a
47-pixel template spans about 13 $\lambda/D$, it was not promoted directly.

The completed
[contrast-linearity checkpoint](../klip-response-tail-linearity-setup-20260921/README.md)
repeated six geometry-only sites per inner radius at half and twice the
original perturbation. All gates pass. Mode-200 median tail cosine is
0.9676--0.9982 and projection is 0.9463--1.0189. Across all eight modes the
worst median cosine is 0.9673 and the worst median projection error is 5.37%.
The Richardson 47-pixel response has complete support and remains below 0.60%
median and 1.45% individual border energy. This promotes the 47-pixel
footprint.

The completed
[47-pixel exact-response campaign](../klip-response-47-setup-20260921/README.md)
regenerated all 11,192 integer search locations with the frozen archived
binary. Its baseline is bitwise identical, coordinates are unchanged, and
every central 11-pixel response exactly equals the archive. The full responses
reproduce the independent fixed-site derivatives with minimum cosine 0.999932
and maximum projection error $6.45\times10^{-5}$. This validates the larger
templates for Stage B.

### Stage B: noise-only covariance screening

Use only the signal-free baseline and deterministic angular block splits. Fit
on one set of patches and project weights onto disjoint held-out patches. Swap
the blocks and repeat. Evaluate every mode and radius, with mode 200 and radii
7.5, 10, and 12 controlling promotion.

The completed
[footprint preflight](../klip-stage-b-footprint-preflight-setup-20260923/README.md)
rejects a coupled support rule. Eleven-pixel training is viable, but 31-pixel
training yields only 0--16 accepted patches over the full radial range and
leaves at least one split empty. Forty-seven-pixel training yields no accepted
patches at any selected query. This is a geometric limit of the image and
source exclusion, not only a rank problem.

Response support and covariance-estimation support are therefore separate
axes. Screen central 11-, 31-, and 47-pixel exact responses. Direct empirical
PCA and diagonal covariance remain on 11-pixel data. The stationary PSD path
uses 11-by-11 Welch training patches with five-pixel angular and radial center
spacing, while the candidate data, response, and five-source exclusion retain
the chosen response support.

The completed
[decoupled-footprint preflight](../klip-stage-b-decoupled-preflight-setup-20260923/README.md)
passes. The 31- and 47-pixel response exclusions require a 40-pixel training
half-width at radii 7.5, 10, and 12. With the 47-pixel exclusion, the full range
still supplies at least 183--206 patches and 58--76 in each detector half at
those radii. The full-range band remains the fixed wider stability control.

Test these axes without positive injections:

| Axis | Fixed grid |
| --- | --- |
| Response support | Central 11, 31, and 47 pixels of the exact field |
| Covariance-estimation support | Direct PCA/diagonal: 11 pixels; stationary PSD: 11-pixel Welch patches applied behind each response support |
| Training radial half-width | 0, 5, 10, 20, 40, and 60 pixels; promote the narrowest split-supported band and retain the next wider band as a stability control |
| Patch coordinates | Raw; pixelwise radial-variance standardized |
| Empirical PCA | Eleven-pixel response only; retain 0, 3, 8, and every estimable mode; floor fractions 0.1, 0.3, and 1.0 |
| Diagonal covariance | Eleven-pixel response only; individual fitted pixel variances with the same fitted mean |
| PSD window | Rectangular and separable Hann on 11-by-11 training patches |
| PSD isotropic mixing | 0.1 and 0.3 of the unwindowed mean pixel variance |
| Patch-amplitude control | Post-ensemble-mean per-patch RMS normalization for selected PSD geometries |
| Fitted mean | On and off for the 11-pixel arms; a support-independent radial-mean control for larger responses |

The
[PSD extension contract](../klip-stage-b-psd-extension-check-20260923/README.md)
calculates each 11-pixel periodogram directly on the 21-, 61-, or 93-pixel
linear-lag grid of the selected response. It retains measured lags -10 through
+10 and sets longer lags to zero. Trace rescaling and isotropic mixing preserve
a positive spectrum. Zero-padded FFT convolution applies the resulting finite
block-Toeplitz covariance, and preconditioned conjugate gradients solve its
weights without a dense 961- or 2,209-component inverse. The implementation
exactly reproduces the existing 11-pixel P4 estimator and passes symmetry,
positivity, isotropic-endpoint, and residual checks at all three supports.
Candidate data and templates remain untapered.

The completed
[local mode-200 raw PSD screen](../klip-stage-b-local-noise-screen-20260923/README.md)
verified that this calculation can run from a 2.8-MB response/baseline bundle
instead of transferring the 1.6-GB field. All 1,080 frozen geometries replayed
exactly and 10,800 policy/query records completed locally in 285.8 seconds.
Hann/mixing-0.1 is the leading raw PSD family: over the controlling radii its
narrow-band candidate score variance is 2.23, 2.56, and 2.22 for response
supports 11, 31, and 47, compared with 6.16, 6.83, and 6.00 for identity.
Rectangular/mixing-0.3 gives 2.63, 3.04, and 2.80 and remains the mandatory P4
prior.

These candidates are not promoted yet. For the 11-pixel Hann/mixing-0.1 arm,
opposite-half generic patches have median measured/predicted variance 1.13,
while the fixed candidate nulls have variance 2.23. The mismatch peaks at
radius 12 and remains in center-only scores. Fitted-mean subtraction does not
remove it, and the 47-pixel response is not uniformly better than 11 pixels.
The completed
[strict radial-normalization arm](../klip-stage-b-radial-normalization-20260923/README.md)
reduces the 11-pixel Hann/mixing-0.1 primary candidate variance from 2.23 to
1.86 while leaving opposite-half variance near unity and split physical-weight
cosine at 0.990. Radius 12 remains high at 3.10, so radial scale is only a
partial explanation. The leave-site-out 3.6-pixel profile is supported at all
72 11-pixel sites, only 22 of 72 31-pixel sites, and no 47-pixel site. Larger
supports therefore cannot use this strict rule without reading the candidate
or extrapolating through an unsupported annulus. The next focused arm is
post-ensemble-mean patch-RMS normalization on the supported 11-pixel geometry,
followed by the remaining direct-covariance and mean controls.

The subsequent
[known-planet footprint audit and correction](../klip-stage-b-planet-mask-correction-20260924/README.md)
found that these first candidate screens excluded the fitted planet from
training and radial-profile estimation but did not mask outer candidate and
response pixels that entered its seven-pixel disk. At radius 12, complete
planet-clear five-query sites number 9 of 12 for the 11-pixel response, 2 of 12
for 31 pixels, and 0 of 12 for 47 pixels. The earlier 31- and 47-pixel
candidate-score comparisons are therefore superseded.

The corrected screen preserves the frozen sites, masks planet-disk coordinates
in both candidate data and exact responses, and solves the matching covariance
principal submatrix. Corrected raw narrow-band Hann/mixing-0.1 primary variance
is 2.477, 2.875, and 2.698 for supports 11, 31, and 47; 11 pixels is the best of
the three. Consistent radial normalization lowers the 11-pixel result to 2.092,
while radius 12 remains high at 3.153. A geometry-only planet-clear subset is
worse at radius 12, so the known planet does not account for that mismatch.
Only the corrected masked values enter subsequent Stage-B decisions.

The subsequent
[post-mean patch-RMS control](../klip-stage-b-patch-rms-20260925/README.md)
changes neither coordinate system materially. In raw narrow-band coordinates,
rectangular/mixing-0.3 changes from variance 2.867 to 2.865 and
Hann/mixing-0.1 from 2.477 to 2.482. After strict radial standardization, those
changes are 2.412 to 2.403 and 2.092 to 2.089. Paired parent and patch-RMS
scores correlate above 0.9997, while radius-12 radial Hann slightly worsens from
3.153 to 3.168. Per-patch RMS is therefore removed from the shortlist. The
direct 11-pixel diagonal/PCA grid is the next discriminator.

The completed
[direct diagonal/PCA screen](../klip-stage-b-direct-covariance-20260925/README.md)
also fails to displace the PSD models. Direct diagonal variance is 6.786 raw and
5.636 after radial standardization. The best direct model, every estimable PCA
mode with floor 1.0, reaches 3.439 raw and 2.886 radial, compared with 2.477 and
2.092 for Hann/mixing-0.1. Its radial split-weight cosine is only 0.730 and its
opposite-half variance is 1.820. A wider band improves those diagnostics but
leaves candidate variance at 2.909. Direct covariance is removed from the
shortlist, and the larger-response radial-mean control is the remaining Stage-B
arm.

The completed
[support-independent radial-mean control](../klip-stage-b-radial-mean-20260925/README.md)
closes the final Stage-B arm. One one-pixel annular profile from the
signal-free image is applied to all three response supports while every
covariance and response quantity remains fixed. For narrow-band
Hann/mixing-0.1, primary variance changes from 2.477 to 2.446 at 11 pixels,
2.875 to 2.831 at 31 pixels, and 2.698 to 2.654 at 47 pixels. Score
correlations exceed 0.9998, and radius 12 remains strongly undercalibrated.
The mean model is therefore neutral and no larger response is promoted.

Stage B is closed. Carry radial-standardized 11-pixel Hann/mixing-0.1 as the
leading data-selected arm and raw 11-pixel rectangular/mixing-0.3 as the
mandatory P4-prior arm into development. Do not carry patch-RMS, direct
diagonal/PCA covariance, or radial-mean subtraction. The frozen Stage-C
development injections are next.

For every fit report:

- training and held-out amplitude variance divided by conditional prediction;
- held-out centered mean-square error;
- fraction within one and two conditional standard deviations;
- split-weight cosine and amplitude-projection agreement;
- covariance spectrum, condition number, retained modes, and fitted samples;
- native-pixel interpolation gain; and
- signed held-out score quantiles and maximum gaps.

Promotion is based on held-out prediction and split stability, not the known
planet. Carry at most two data-selected covariance families into development.
Regardless of that ranking, also carry the P4 priors of raw rectangular PSD,
0.3 isotropic mixing, radial half-width 20, with mean-variance clipping and
half-mean hard truncation. This prevents KLIP results from depending entirely
on another selection made from the same baseline.

### Stage C: development injections

Analyze the 36 development sites at all three source levels. Covariance is fit
on each positive image with the complete source/search exclusion. In parallel,
fit weights once on the signal-free baseline and replay them unchanged on each
positive. This frozen-weight arm measures source contamination of the fitted
mean and covariance; it is a diagnostic rather than a deployable science mode.
Before any development positive is opened, calculate the maximum-null threshold
for every development method from its 20 fixed calibration locations. Those
thresholds are not updated from injection recovery.

For the promoted PSD model, test the complete precision grid used in P4:

- the full inverse;
- eigenvalue clipping at 0.5, 0.75, and 1.0 times mean covariance variance;
- hard truncation below 0.5, 0.75, and 1.0 times mean covariance variance.

Renormalize every weight to unit response after regularization. Record the
number of modified or retained modes and the exact fitted-covariance SNR
efficiency. Evaluate the exact response first. Apply the sparse response to
identity and the final covariance shortlist to detect response-by-covariance
interactions without repeating the whole estimator grid.

The development report includes, by radius, level, and mode:

- recovery at the frozen maximum-null threshold;
- mean and median five-pixel maximum SNR;
- center-pixel SNR and localization error;
- raw amplitude bias and scatter;
- positive-minus-baseline throughput;
- calibration-null threshold and tail diagnostics; and
- refit-minus-frozen changes.

Freeze at most two covariance candidates after this stage. Selection must use
the complete table and may not be based on the actual planet or on the best KL
mode. One candidate may optimize stable faint-source recovery and the other
mean SNR, mirroring the clipped-versus-truncated P4 outcome.

### Stage D: threshold audit and policy freeze

For every frozen method, radius, and mode, copy its unchanged development
threshold from the 20 preassigned calibration nulls into the policy receipt
before opening validation positives. Record the full signed null distribution,
largest and second-largest scores, and their gap. Use deterministic paired
resampling and angular-block deletion to report threshold and development
recovery sensitivity.

The maximum of 20 correlated locations is retained for direct comparison with
P4, but it is not interpreted as a formal false-alarm probability. Report the
effective location geometry and do not call the calibration samples
independent. The six held-out null locations per radius remain unopened until
the policy receipt is written.

The policy receipt contains exact method parameters, response choice, mode-200
primary status, candidate locations, contrasts, thresholds, software hashes,
and fingerprints of every input. Validation must refuse to start if any receipt
input has changed.

After the receipt is durable, expose the held-out null locations and report
their exceedances under the unchanged thresholds, including sensitivity to the
paired calibration resamples. A poor held-out result can reject a frozen
candidate but cannot be used to alter it or nominate another method from the
development grid.

### Stage E: fresh validation

Run and analyze the 36 validation sites only after the policy receipt exists.
No covariance, response, smoothing width, precision cutoff, threshold, search
aperture, or source level changes are allowed.

The primary comparison is paired against Gaussian FWHM 3.6 at mode 200. A
covariance method is accepted as a KLIP detection improvement only if it:

1. has no more held-out-null exceedances than Gaussian;
2. has at least Gaussian recovery at every source level and a strictly larger
   total recovery count;
3. has a positive paired mean-SNR change whose 95-percent radius-stratified
   site-bootstrap lower bound is above zero at either SNR 3 or SNR 5, without
   a significant loss at another level;
4. retains positive-minus-baseline throughput consistent with the finite
   response audit; and
5. does not obtain its advantage by rejecting locations valid for Gaussian.

Report the other seven modes separately as a generalization check. A result
confined to one KL mode is recorded as such rather than generalized to KLIP.

### Stage F: known planet and closure

After validation is immutable, apply all frozen finalists and permanent
references to the original planet-bearing cube using `working/analyze.conf`.
Report nearest-pixel and aperture-maximum SNR, peak position, amplitude, and
support for every mode. This endpoint checks consistency with earlier KLIP
work; it cannot select or rescue a method that failed validation.

The closure report places Gaussian, identity-response, covariance candidates,
and response-smoothed controls side by side for injections, held-out nulls,
and the planet. It must state whether the original matched-filter deficit came
from covariance mismatch, response structure that remains harmful after
whitening, or neither tested model being adequate.

## Implementation and execution checkpoints

1. **Inventory and compatibility checker.** Freeze and verify the archived
   exact/sparse products, current binaries, configuration, PSF, input list, and
   baseline equivalence.
2. **KLIP template adapter.** Load schema-2 47-pixel `PIXEL_EXACT` response
   stamps for all eight modes, expose fixed central 11/31/47 supports,
   transform templates consistently with rotated training patches, and
   reproduce `hciAnalyze` identity filtering exactly.
3. **Noise-only suite.** Freeze support-aware half-overlap sampling, lag grids,
   and effective-rank diagnostics before fitting. Reuse the tested P4
   covariance algebra behind a KLIP mode/template interface; write strict JSON
   diagnostics and a promotion report.
4. **Immutable injection preparer.** Partition sites, calibrate mode-200
   Gaussian SNR 3/5/7 contrasts, write commands, and freeze all fingerprints.
5. **Resumable development runner.** Launch the 108 development reductions,
   analyze all modes, run frozen-covariance replay, and write the precision and
   response comparisons.
6. **Policy freeze.** Write the selected methods and calibration thresholds to
   a separate immutable receipt.
7. **Resumable validation runner.** Launch the 108 unopened validation
   reductions, analyze them without tuning, and expose held-out nulls.
8. **Planet and final report.** Run the descriptive planet endpoint, verify
   independent annular SNR, and commit compact tables and figures.

The prototype may remain in Python because production `hciAnalyze` currently
supports identity, diagonal, and same-radius PCA but not the pooled PSD or
precision-regularization grid. Promote only a validated final estimator into
C++ and add its configuration surface after the scientific comparison.

Checkpoint 4 is implemented by the
[immutable Stage-C preparer](../klip-stage-c-development-setup-20260925/README.md).
It re-audits the full all-mode candidate pool, chooses the smallest fixed
training band that supports all 38 roles at each radius, partitions the sites
by deterministic angular maximin selection, and freezes site-specific
Gaussian-SNR 3/5/7 contrasts for both development and validation. It also
writes the 108 development commands and freezes the complete 18-method
analysis receipt. The prepared ROC geometry and contrasts were reviewed before
the resumable runner was enabled; validation positive products remain unopened.

Checkpoint 5 is complete. The
[Stage-C development result](../klip-stage-c-development-20260926/README.md)
preserves all 108 analyses and selects a compact validation shortlist without
opening held-out nulls or validation images. The development result favors the
identity response filters over covariance weighting, retains raw rectangular
as the mandatory prior, and selects radial-Hann truncation at 0.75 as the one
radial covariance policy for the policy-freeze step.

Checkpoints 6 and 7 are implemented by the
[Stage-D policy and Stage-E validation setup](../klip-stage-d-e-validation-setup-20260926/README.md).
The policy runner copies the unchanged selected-method thresholds, full signed
calibration distributions, deterministic paired bootstrap and angular-block
sensitivity, exact validation tasks, software, and all transitive input
fingerprints into a receipt before a held-out score can be read. The validation
runner requires that receipt, freezes baseline-only weights for the unopened
sites, exposes the preassigned held-out nulls, runs the 108 validation
reductions, and evaluates the preregistered mode-200 gate without method or
mode maximization. A read-only ROC preflight found complete generic annular
coverage for every validation and held-out search pixel. A guarded repair
corrected a driver-only two-versus-three-value unpacking error before any
Stage-E product opened. The resumed run completed all 48 model units, 36
held-out nulls, and 108 validation reductions. Analysis then stopped in the
response-fidelity diagnostic because it imported the original masked-support
helper instead of the recorded Stage-C repair. The filter results are
unaffected and no final validation summary exists. A guarded analysis-only
repair retains all upstream products, archives all 108 partial analysis
directories, and reruns analysis with the repaired full-space convention.

## Required provenance and failure behavior

- Every stage writes `protocol.json`, `manifest.json`, `state.json`, and a
  completion receipt with SHA-256 fingerprints.
- JSON output is strict: nonfinite values are encoded as `null` with a separate
  status or cause field.
- A failed or partial task is archived before retry; completed task directories
  are verified and reused.
- Baseline, development, and validation products remain separate directories.
- Candidate support, annular validity, template availability, and covariance
  support are verified on common masks before any summary is calculated.
- The independent annular oracle must agree with `hciAnalyze` for every
  permanent reference and final candidate at the implemented `rtol=2e-6`,
  `atol=2e-6`; inner-edge substitutions are separately enumerated.
- Raw KLIP inputs, exact response products, positive reductions, and policy
  receipts are never modified in place.

## Interpretation limits

All locations share one observing sequence and correlated residual field.
Development and validation positions prevent direct score reuse but do not
constitute independent observing epochs. The radius-stratified site bootstrap
quantifies spatial variation within this data set, not survey-to-survey
performance. A production claim ultimately needs another target, epoch, or
observing sequence.

The program tests covariance of final KLIP images. It does not model temporal
covariance before derotation/combination and does not combine KL modes as
independent measurements. Those are later extensions if the final-image
prototype demonstrates a reproducible gain.
