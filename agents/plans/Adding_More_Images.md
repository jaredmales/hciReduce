# Adding More Images

I'd like to extend the detector-frame P4 alogorithm to include additional images.  Right now we predict the pixel value based only on the pixels in the same image.  The new extension will include pixels from an image before and an image after.  For each target image, we will base which images to include on the psfRradius, such that we include the first image which has rotated enough to move by one psfRadius at the separation of the inner edge of the SR (not the OR).  Images at the beginning and end that do not have both before and after images available will be excluded from the target images.

Formulate a plan and fill it in below.  Do not alter this prompt.

## Update After Initial Tests:

Ok let's degrade gracefully rather than fail.  Switch to basing rotation amount off the mean-radius, not hte inner radius, of the SR.  Then if the psfRadius spec can't be met for an SR annulus, proceed with the maximum psfRadius that can be used for that SR.  If there is no valid psfRadius (which can happen for the close-in SRs) proceed with numberImages=0 for that SR annulus.

# Plan

1. Define `p4.numberImages` as a nonnegative integer P4 configuration option, defaulting to `0`.  Preserve the
   current detector- and rotated-frame behavior bit-for-bit when it is `0`; values greater than zero enable
   multi-image detector-frame predictors.  Reject an invalid or negative value during configuration validation and
   document the option, its default, and its compatibility contract in the P4 configuration reference and user guide.

2. After the target cube and derotation angles are available, build a deterministic temporal-neighbor selection for
   each central target image.  For a search annulus whose mean radius is `r`, require a field-rotation displacement
   of at least one physical `p4.psfRadius`, i.e. `abs(angleDiff) * r >= p4.psfRadius` (with angles in radians), using
   the existing wrapped angular-difference utility.
   Scan independently backward and forward in input-image order; on each side retain the nearest `p4.numberImages`
   frames satisfying that threshold.  Use the absolute angle difference so the two temporal directions work whether
   the sequence's parallactic angle rises or falls.  Do not wrap around either end of the cube.

3. Treat a central image as usable in an annulus only when it has all requested qualifying earlier and later
   neighbors.  Omit other central images from the local PCA samples and leave their residual/validity output pixels
   invalid. If the requested radius leaves no structurally valid PCA sample set, lower it to the greatest positive
   radius that does. If no positive radius works (including close-in annuli with no meaningful rotation), fall back to
   `numberImages=0` for that annulus. Because the threshold is evaluated at each annulus's mean radius, the usable
   target set and effective temporal radius may differ by annulus.

4. Extend the detector-frame regression assembly in `P4Reduction::regions()` so that each usable central-image row
   contains its original `K` OR-pixel predictors plus `K` predictors for every selected earlier and later image, for
   a total of `K * (1 + 2 * p4.numberImages)` columns.  Keep the response as the central image's search-pixel value.
   Feed the compacted usable-target rows to the existing PCA path, then scatter residuals only to their original
   central-image planes.  Recompute structural degrees of freedom and realized mode limits using the compacted row
   count and expanded predictor count.

5. Keep rotated-frame behavior explicitly unchanged for this increment.  Reject the unsupported configuration
   `p4.regressionFrame=rotated` with `p4.numberImages > 0` during validation; document that explicit contract rather
   than silently ignoring the requested setting.

6. Add provenance and diagnostics sufficient to reproduce the selection: record the configured and effective
   per-annulus temporal image counts and radii in output headers and, when diagnostics are enabled, write
   per-annulus target-selection/neighbor-index data. Update region statistics and progress text to distinguish usable
   target rows from the full cube length and to report the expanded predictor count.

7. Add Catch2 coverage at the P4-reduction and application levels for: the default `0` compatibility path; option
   registration, loading, and invalid values; increasing and decreasing angle sequences; exact threshold inclusion;
   nearest qualifying selection; `N > 1`; endpoint truncation without wraparound; per-annulus mean-radius behavior;
   reduced effective-radius and same-image fallback paths; residual placement/invalidity for dropped targets;
   expanded predictor count and mode limits; output provenance; and rejection of the unsupported rotated-frame
   combination. Update user-facing P4 documentation and example configuration with a nonzero usage example.
