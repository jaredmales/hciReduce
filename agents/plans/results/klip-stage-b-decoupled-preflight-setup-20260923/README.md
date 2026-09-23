# KLIP Stage-B decoupled-footprint preflight setup

## Question

The coupled-footprint audit found no usable 47-pixel covariance-training stamps
and inadequate 31-pixel coverage. This follow-up tests the geometry required by
a stationary PSD model: estimate the noise spectrum from 11-by-11 Welch patches
but apply it behind 11-, 31-, or 47-pixel response templates.

This remains a score-blind geometry calculation. It does not fit a PSD or
calculate a matched-filter value.

## Fixed calculation

The response support controls candidate data validity and the source-exclusion
footprint. The training support remains 11 pixels in every arm, with five-pixel
angular and radial center spacing. Thus each candidate excludes the union of
its five 11-, 31-, or 47-pixel response footprints, while accepted training
patches themselves remain 11-by-11.

The calculation reuses the deterministic sites from the completed footprint
preflight and audits all five search pixels. It tests radial half-widths 0, 5,
10, 20, 40, and 60 pixels. For every radius and response support it reports:

- accepted 11-pixel training patches;
- counts wholly contained in each disjoint detector half;
- the centered covariance-rank ceiling; and
- the narrowest band with at least eight patches in both halves for every
  selected query.

The result decides whether an 11-pixel Welch PSD can be estimated without
using pixels touched by a larger candidate response. A pass does not assume
that correlations beyond the 11-pixel lag range vanish; that spectral
extension and its controls will be frozen separately in the noise-only suite.

## ROC commands

```bash
git pull

root=working/roc/klip_stage_b_decoupled_preflight_20260923

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_b_decoupled_preflight.py check

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_b_decoupled_preflight.py prepare "$root"

taskset -c 12-27 python3 \
  "$root/software/run_klip_stage_b_decoupled_preflight.py" run "$root"
```

The runner verifies the complete response and parent-preflight lineage before
analysis. It writes `results.json`, `results.md`, and `complete.json`.
