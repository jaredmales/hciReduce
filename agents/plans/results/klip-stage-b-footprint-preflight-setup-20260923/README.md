# KLIP Stage-B footprint preflight setup

## Question

The accepted KLIP response now contains 47-by-47 pixels, or 2,209 template
components. Before estimating covariance, this preflight determines how much
response is lost by central 11- and 31-pixel crops and whether each footprint
has enough spatially distinct covariance-training patches at the six planned
radii.

This is a geometry and support audit. It does not calculate matched-filter
scores, covariance spectra, or source recovery, so its outcome cannot favor a
method based on the baseline residual values.

## Fixed analysis

For every KL mode and radii 6, 7.5, 10, 12, 16, 20, and 24 pixels, report:

- complete response counts for central 11-, 31-, and 47-pixel supports;
- captured response energy relative to the complete 47-pixel stamp;
- border energy at each support;
- common eight-mode locations supporting the center and four cardinal search
  pixels; and
- the count and centered-rank ceiling of covariance-training patches.

Candidate geometry excludes the fitted planet by seven pixels and requires the
candidate data and response to be complete in every mode. Up to 12 candidate
centers per radius are selected by deterministic angular spacing. Training is
then audited for all five search pixels at every selected center.

For support width $s$, training stamps use radial/tangential bilinear
extraction and half-overlap center spacing $(s-1)/2$: 5, 15, and 23 pixels.
Every fit excludes the known planet and the union of all five candidate
footprints. The tested radial half-widths are:

| Support | Radial half-widths |
| ---: | --- |
| 11 | 0, 5, 10, 20, 40, 60 pixels |
| 31 | 0, 15, 30, 45, 60 pixels |
| 47 | 0, 23, 46, 69 pixels |

The last width reaches every permitted center ring after clipping training
centers to `6 <= r < 60`. The audit also partitions accepted stamps into two
disjoint detector half-planes. It reports the narrowest band with at least
eight patches in each half for every selected query. Patch counts remain
overlapping-sample counts; the report does not interpret them as independent
noise realizations.

## ROC commands

```bash
git pull

root=working/roc/klip_stage_b_footprint_preflight_20260923

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_b_footprint_preflight.py check

taskset -c 12-27 python3 \
  agents/plans/scripts/run_klip_stage_b_footprint_preflight.py prepare "$root"

taskset -c 12-27 python3 \
  "$root/software/run_klip_stage_b_footprint_preflight.py" run "$root"
```

The run verifies and hashes the complete 47-pixel response lineage before
analysis. It writes `results.json`, `results.md`, and `complete.json`.
