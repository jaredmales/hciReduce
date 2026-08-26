# ToDo

This is a list of to-do items maintained by developers.  As we move through each we'll create a separate plan document for the complicated ones.  Agents should not edit this file unless directed, except to mark [x] for complete and add the commit(s) that correspond(s) to that work

## Doc Update:
[x] under user-guide, clicking on "Common Configuration and Workflow" closes the TOC drop down for User's Guide.  I think that means it is marked as the main page or equivalent for User's Guide.  Fix this behavior.  If needed have the main page just be title or another TOC.

[x]  Add a placeholder Introduction page before "Common Configuration...".  It can be empty for now.

Completed in `c06172f` (`Improve outputs, P4 headers, and user guide navigation`).

## FITS Updates:
[x] change FITS keywords to exploit HIERARCH (happens automatically) to make them more expressive and more closely match the config settings.  E.g. P4DRIN should be P4 OR DELTA RADIUS INNER

[x] Add the FITS keywords to the config tables as a column, i.e. to document the mapping from config to header

Completed in `c06172f` (`Improve outputs, P4 headers, and user guide navigation`).

## Output Directories
[x] make sure we automatically create any output directories

Completed in `c06172f` (`Improve outputs, P4 headers, and user guide navigation`).

## P4 Profiling

[x] added timing measurements for various steps as is done in klipReduce.  Consider a more resuable common design for this than the globals in klipReduce.

Implementation plan: `agents/plans/p4_profiling.md`.

Completed in `2407edc` (`Add shared P4 and KLIP timing`) and `178ab0f` (`added showTiming flag`).

## P4 Memory Management

Implementation plan: `agents/plans/p4_memory_management_optimization.md`.

[x] consider whether we can save on needed memory.  Assume that pre-processing is done and written out to disk.  Then we could load only the pixels needed for the annulus, or even a subset, being calculated into memory at a time.  For very large data sets this could be impactful.

[x] provide automatic reduction in number of threads given available resources.  I.e. if we can support all cores for the current annulus, we drop the number of threads.

Completed in `2590b15` (`Avoid P4 sampling image copies`), `284111f` (`Add bounded P4 memory management`), and
`33e3af9` (`Aggregate P4 finalization progress`). Remaining optional input-cube streaming and residual-product
streaming are recorded as deferred work in `agents/plans/p4_memory_management_optimization.md`.

## KLIP Memory Management

Implementation plan: `agents/plans/klip_memory_management_optimization.md`.

[] Reduce KLIP memory use for large target/reference sequences while preserving its numerical and output contracts.
The planned order is compact residual storage and one-mode finalization, sparse reference-selection bookkeeping, then
memory-aware worker limits and selected-covariance optimization.  The investigation and plan were completed in
`f539e95` (`Close P4 memory work and plan KLIP memory`); implementation is intentionally deferred.

## P4 PSF calculation and Post-Processing

Implementation plan: `agents/plans/p4_psf_calculation_post_processing.md`.

[x] add PSF calculation in P4.  For each output pixel, we calculate the predicted output PSF.  This is done by projecting the modal basis for each input pixel onto the PSF, then rotating and stacking.  To be efficient we'll need to calculate the input-pixel PSFs as we go, store them in a compact form (i.e only a small postage stamp of say 10x10 pixels), possibly writing them to disk.  These will be cubes for the modes. Once processing is done we'll create a PSF for each output pixel by selecting and rotating each of the input-pixel PSFs and coadding the same way.

[x] Once we have PSF calculation, we'll then filter the output image with that spatially variable PSF.  That is, we'll replace each output pixel by the value of that pixel after convolution with the PSF that goes with that pixel.

Completed in `cc1c3c2` (`Add P4 PSF model calculation`), `b1d5bd4` (`Tighten P4 local PSF support`),
`42bdeda` (`Add spatially variable P4 PSF filtering`), and `2795905` (`Align P4 filter products with final images`).

[x] Extend analytic PSF calculation to `p4.numberImages>0` with bounded coefficient/stamp storage.  The implemented
factorization retains one same-image stamp plus the temporal PSF-disk coefficients, reconstructs each target with its
exact selection and detector position against the full template, and does not retain `pixels * frames * modes * stamp
pixels`.

[] Add a one-shot complete-field post-preprocessing injection path for the fitted negative source.  Use the resulting
signal-free reduction to regenerate the spatially variable PSF and filtered product, then determine whether the
response-normalization trough and filtered ring at the true companion separation remain.

[] Complete PSF-filter validation with isolated synthetic sources and application-level real-FITS coverage.  Check
response localization, detection ranking, validity at annulus/image boundaries, and preservation of ordinary
unfiltered results when PSF calculation/filtering is disabled.

## Evaluate extend of double vs float in P4.

[] In P4 double precision is carried through to more than just matrix decomposition.  We should evaluate if there are performance gains for switching to float for more of the algorithm.

## P4 Pixel Local Processing 

Implementation plans: `agents/plans/p4_pixel_local_processing.md` and
`agents/plans/p4_negative_companion_optimization.md`.

Fake planet injection in P4 will exploit the pixel locality of P4 processing.  That is, only a subset of pixels need to be processed to produce a given output pixel.  This subset is the fractional pixel in each input image that rotates to the position of the ouput pixel, plus its interpolation support.  

[x] Define the configuration for this:
   -- for a given fake planet injection, specify the size of the region to process.  Possibly this is handled by the specified PSF itself (i.e. if its 10x10 then we know it's 10x10).  But be able to override this so a large PSF file can be specified and automatically cut down
   -- that might be all we need in addition to existing fake injections

[x] Implement the book-keeping to calculate which input pixels need to be processed for a given PSF.
   -- This should merge the subsets for all output pixels needed for the specified PSF, i.e. if it's 10x10 then we have 100 ouput pixels.  That should only be 10x10 input pixels in each image plus the interpolation support (roughly 2 pixels around the outside).  So the result should be something like 196 pixels x number of images.
   -- optionally include processing for the output pixels that are affected by the presence of the planet.  That is they would have the planet somewhere in their input so it would slightly change their basis.  Default is off b/c we are assuming that the PSF mask minimizes this impact in P4.

[x] For each fake planet, process only the defined pixels.  We nominally only need to save the resulting output pixels as mode-cube.  This can be just the pixels with FITS headers specifying their actual location in the image for later post-processing.  Optionally a new full output cube can be produced replacing the pixels.

Completed in `997d986` (`Add P4 pixel-local processing`).

[x] Add bounded negative-companion contrast optimization using repeated in-memory pixel-local evaluations, with
dense-basin validation and diagnostic/best-fit products.

[x] Extend negative-companion optimization to joint Cartesian position and contrast while preserving the
`0.5*(N-1)` center convention and a fixed sky merit aperture.

[x] Add contiguous delete-one-time-block jackknife covariance and make joint refits search around the converged
complete-sequence solution while retaining the original configured coordinate origin in products.

Completed in `2d6d7d7` (`Add P4 negative companion optimization`), `e05c5ae` (`Add full P4 negative companion
optimization`), `7325ae0` (`Fix P4 joint optimizer aperture support`), `f548592` (`Add P4 jackknife uncertainty and
guide pages`), and `96dc547` (`fixed jackknife failure mode`).  A 621-frame, eight-block real-data run completed all
refits and produced the fitted result and covariance products.

[] add interleaved mode as an alternative to contiguous blocks in jackknife.

[] Validate jackknife uncertainty with negative-source injection/recovery trials at nearby empty locations and
comparable separations.  Compare several reasonable block counts, and compare contiguous with interleaved blocks once
implemented, before adopting either result as a calibrated science uncertainty.

## Better Temporal Prediction

[] Try just adding 1 pixel instead of the whole PSF mask

[] Try adding the mean or median of the pixels

[] when adding numberImages, the extension at the beginning and end by adding 2 later images (for images at the beginning) or 2 earlier images (for images at the end) should be modified.  We should calculate 3 different sets of coefficients.  This is to optimize to exploit any temporal predictability.  As it is, with a single set of coefficients for all three subsets we're basically only calculating the mean.

[] alternatively, we could use PCAT (Long et al, 2023,https://arxiv.org/abs/2303.05559).The flaw in this is that the data have gaps and may not be regularly sampled, meaning that we are not fully exploiting temporal correlations
  - calculate a KL transform of the time-series with a rotation gap using all the other pixels in the OR (or maybe a temoral-OR)
  - project that onto the time-series of the pixels in the masked region
  - use the predicted values of the masked pixels in the main P4 calculation instead of the numberImages pixels above.

[] alternatively, we could use linear prediction.  This could be done to take into account the proper time axis:
  - form an estimate of the auto-correlation (AC) or PSD of the intensity by averaging the AC or PSD of the pixels in the OR. we want to use something like Lomb-Scargle here if the data are not regularly sampled over the whole time-series, which will be common for real observations.  Whether it's AC or PSD is an open question
  - Use the AC/PSD to calculate optimum LP coefficients
  - Use the LP to predict the values in the PSF-masked patch using only prior/later pixels that are outside the rotation exclusion.  
     - This will require multi-step prediction, over non-regular sampling.  We'll need strategies to handle this.
     - There may be better conceptualizations that LP for interpolation vs. extrapolation that exploit the same statistical predictability.
  - Alternatively: use the PSD to calculate a KL transform of the Legendre polynomials and use this as the basis

## Exclusion of target images

Right now we do an inversion for all the images at once.  This means that when we are predictiog x_t[n] the pixels' in OR(x_t) are included in the basis, as is the true value of x_t[n] included in the b vector

[] implement excluding the target image from the basis set
   - ideally we should exclude x_t from it's own basis set
   - if we do this we may also consider having a configurable temporal exclusion region as we do with numberImages > 0.
   
[] this will be much slower, probably on the order of T times slower. So we should look into 
      -- using the downdate algorithm of Long and Males 2021 (https://arxiv.org/abs/2101.11634).  This will require restructuring the algo to use SVD of the data instead of the CV.
      -- batch GPU processing.  I.e. load the CV into GPU memory once then do the decomposition T times, once for each t removed.   

## Other solvers

[] SVD (which is equivalent to eigendecompositon for a covariance matrix)

[] In particular we should look at SVD + Tikhonov regularization.  

## Pixel grouping

[] evaluate grouping pixels such that the solution is based on, say, a 2x2 square of pixels instead of a single pixel.  That solution is used for each pixel individually
  - This may reduce calculation time
  - This may improve the fit given more data

## Fake Planet Cleanup

[x] enable fake injection with preProcess.skip=true

[x] extend the existing "planet" configuration section, parallel to "fake", so the known location and contrast of a planet can be propagated through to hciAnalyze, etc.

[x] add a `fake.subtractPlanet=false/true` flag. When true, inject the negative of the configured planet in addition to
the configured fake planets, without adding it to the fake vectors or `FAKE*` metadata.

## Config Cleanup

[x] make it so boolean options work with both.  This can be done based on isSet.
   [x] --flag  (implicitly true)
   [x] --flag=true/false (explicitly set)
