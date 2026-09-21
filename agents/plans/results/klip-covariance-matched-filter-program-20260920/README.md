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

The expensive response calculation is already complete at
`working/roc/klip_signal_free_pixel_response_20260914T232006Z`. It supplies:

- an optimized KLIP negative companion at separation 12.3877505 pixels,
  position angle 260.643152 degrees, and contrast 0.0045743625;
- a signal-free eight-plane baseline at KL modes
  125, 150, 175, 200, 225, 250, 300, and 350;
- exact 11-by-11 `refitDifference` response and validity stamps at all 11,192
  native pixels in the 6-to-60-pixel search annulus;
- the production sparse radial response for a response-model control; and
- the original planet-bearing science cube for the final descriptive endpoint.

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
- Use the exact signal-free per-pixel response as the primary template.
- Use the sparse candidate-avoiding radial response only as a response-model
  control.
- Keep the 11-by-11 response stamp and audit edge energy before inference.
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

The same preflight confirms that the fixed partition is feasible. Requiring
all five search centers to have complete exact responses and to lie outside the
seven-pixel fitted-planet exclusion leaves 42, 44, 60, 91, 113, and 160 centers
in the six primary bins, respectively. Site selection uses only this geometry
and validity metadata, never baseline or positive scores.

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
does not repeat the 16.6-hour exact-response calculation.

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

### Stage B: noise-only covariance screening

Use only the signal-free baseline and deterministic angular block splits. Fit
on one set of patches and project weights onto disjoint held-out patches. Swap
the blocks and repeat. Evaluate every mode and radius, with mode 200 and radii
7.5, 10, and 12 controlling promotion.

Test these axes without positive injections:

| Axis | Fixed grid |
| --- | --- |
| Training radial half-width | 0, 5, 10, and 20 pixels; five-pixel ring spacing |
| Patch coordinates | Raw; pixelwise radial-variance standardized |
| Empirical PCA | Retain 0, 3, 8, and every estimable mode; floor fractions 0.1, 0.3, and 1.0 |
| Diagonal covariance | Individual fitted pixel variances with the same fitted mean |
| PSD window | Rectangular and separable Hann on 11-by-11 patches |
| PSD isotropic mixing | 0.1 and 0.3 of the unwindowed mean pixel variance |
| Patch-amplitude control | Post-ensemble-mean per-patch RMS normalization for the selected PSD geometries |
| Fitted mean | On and off for the identity and selected PSD arms |

The PSD estimator uses the P4 contract: a 21-by-21 zero-padded FFT, all finite
lags from -10 through +10, trace rescaling to the unwindowed sample variance,
and untapered candidate data and templates. Radial bands are clipped at the
6-to-60-pixel supported range, and the contribution of each center ring is
recorded.

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
2. **KLIP template adapter.** Load schema-2 `PIXEL_EXACT` response stamps for
   all eight modes, transform templates consistently with rotated training
   patches, and reproduce `hciAnalyze` identity filtering exactly.
3. **Noise-only suite.** Reuse the tested P4 samplers and covariance algebra
   behind a KLIP mode/template interface; write strict JSON diagnostics and a
   promotion report.
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
  permanent reference and final candidate at `rtol=1e-6`, `atol=1e-7`.
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
