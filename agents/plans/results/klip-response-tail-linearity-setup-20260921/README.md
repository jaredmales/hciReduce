# KLIP response-tail contrast-linearity setup

## Question

The response-stamp convergence experiment found that 47 pixels is the first
footprint with complete selected-site support that clears the original edge
thresholds in every KL mode. At the three inner radii, however, the 31-to-47
shell contains 4.5--7.9% median squared response energy and as much as 11.2%
at an individual site. Before generating a full 47-pixel response field, this
experiment tests whether that outer structure is a contrast-linear KLIP
response or a finite-difference floor.

## Fixed design

The completed response-stamp convergence run is the immutable parent. At radii
7.5, 10, and 12 pixels, the experiment takes alternating even-index sites
0, 2, 4, 6, 8, and 10 from the parent's geometry-only angular selection. This
gives six distributed sites per radius without looking at response values.

The parent supplies the positive/negative central difference at contrast scale
1. Each of the 18 sites receives four new reductions: positive and negative
perturbations at scales 0.5 and 2. Thus 72 new reductions produce three
derivative estimates at every site. The known planet is subtracted identically
in every reduction, and all eight KL modes are retained.

Comparisons use a 63-pixel extraction split into square regions:

- the 11-pixel core;
- the shell from 11 through 31 pixels;
- the controlling shell outside 31 and through 47 pixels; and
- the residual shell outside 47 and through 63 pixels.

For each region, contrast, site, radius, and mode, the analysis records common
support, cosine, projection onto the scale-1 derivative, best-scaled residual,
relative difference, and energy ratio. It also records 47-pixel border energy
at all three contrasts, the positive/negative midpoint departure from the
signal-free baseline, and a zero-contrast Richardson estimate

$$
D(0) \simeq \frac{4D(\epsilon/2)-D(\epsilon)}{3}.
$$

Mode 200 controls promotion. At each radius and both new scales, the median
31-to-47-shell cosine must be at least 0.95 and its median projection must be
within 0.10 of unity. The Richardson 47-pixel response must also retain
complete support, at most one percent median border energy, and at most five
percent at every individual site. The other seven modes are mandatory
diagnostics.

If all gates pass, 47 pixels is promoted for one full-field response campaign.
If a gate fails, the report identifies whether the outer energy is contrast
dependent, incoherent, or still reaches the 47-pixel boundary.

At the measured convergence-run rate, the 72 new reductions should take about
6--8 minutes on the fixed 16-CPU ROC allocation.

## ROC commands

```bash
git pull

root=working/roc/klip_response_tail_linearity_20260921

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_response_tail_linearity.py check

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_response_tail_linearity.py prepare "$root"

taskset -c 12-27 python3 \
  "$root/software/run_klip_response_tail_linearity.py" run "$root"
```

The run is resumable. It verifies the complete parent, all selected scale-1
response products, frozen inputs, executable, CPU contract, FITS headers, and
completed products before reuse. On completion it writes `results.json`,
`results.md`, and a fingerprinted receipt.
