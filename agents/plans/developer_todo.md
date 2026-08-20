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

[] consider whether we can save on needed memory.  Assume that pre-processing is done and written out to disk.  Then we could load only the pixels needed for the annulus, or even a subset, being calculated into memory at a time.  For very large data sets this could be impactful.

[] provide automatic reduction in number of threads given available resources.  I.e. if we can support all cores for the current annulus, we drop the number of threads.

## P4 PSF calculation and Post-Processing
[] add PSF calculation in P4.  For each output pixel, we calculate the predicted output PSF.  This is done by projecting the modal basis for each input pixel onto the PSF, then rotating and stacking.  To be efficient we'll need to calculate the input-pixel PSFs as we go, store them in a compact form (i.e only a small postage stamp of say 10x10 pixels), possibly writing them to disk.  These will be cubes for the modes. Once processing is done we'll create a PSF for each output pixel by selecting and rotating each of the input-pixel PSFs and coadding the same way.

[] Once we have PSF calculation, we'll then filter the output image with that spatially variable PSF.  That is, we'll replace each output pixel by the value of that pixel after convolution with the PSF that goes with that pixel.

## P4 Pixel Local Processing 

Fake planet injection in P4 will exploit the pixel locality of P4 processing.  That is, only a subset of pixels need to be processed to produce a given output pixel.  This subset is the fractional pixel in each input image that rotates to the position of the ouput pixel, plus its interpolation support.  

[] Define the configuration for this:
   -- for a given fake planet injection, specify the size of the region to process.  Possibly this is handled by the specified PSF itself (i.e. if its 10x10 then we know it's 10x10).  But be able to override this so a large PSF file can be specified and automatically cut down
   -- that might be all we need in addition to existing fake injections

[] Implement the book-keeping to calculate which input pixels need to be processed for a given PSF.
   -- This should merge the subsets for all output pixels needed for the specified PSF, i.e. if it's 10x10 then we have 100 ouput pixels.  That should only be 10x10 input pixels in each image plus the interpolation support (roughly 2 pixels around the outside).  So the result should be something like 196 pixels x number of images.
   -- optionally include processing for the output pixels that are affected by the presence of the planet.  That is they would have the planet somewhere in their input so it would slightly change their basis.  Default is off b/c we are assuming that the PSF mask minimizes this impact in P4.

[] For each fake planet, process only the defined pixels.  We nominally only need to save the resulting output pixels as mode-cube.  This can be just the pixels with FITS headers specifying their actual location in the image for later post-processing.  Optionally a new full output cube can be produced replacing the pixels.


## Better Temporal Prediction

[] when adding numberImages, the extension at the beginning and end by adding 2 later images (for images at the beginning) or 2 earlier images (for images at the end) should be modified.  We should calculate 3 different sets of coefficients.  This is to optimize to exploit any temporal predictability.  As it is, with a single set of coefficients for all three subsets we're basically only calculating the mean.
