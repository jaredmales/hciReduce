# KLIP response-stamp convergence result

## Question and design

Stage A found that the archived 11-by-11 exact response leaves 1.35--3.63%
median mode-200 squared energy on its border at every primary radius. This
experiment measured the response on larger footprints before committing to a
new full-field response campaign.

The completed canonical Stage-A run was the immutable parent. Twelve sites at
each of radii 7.5, 10, 12, 16, 20, and 24 pixels were selected using geometry
alone: candidates were sorted by angle and the midpoint of each of 12
equal-count blocks was retained. At each of the 72 sites, one positive and one
equal negative perturbation was applied after subtracting the fitted known
planet. The perturbation half-contrast was 0.0045743624964148452, matching the
archived paired response calculation.

The 144 KLIP reductions supplied one full-image central difference per site.
The same differences were extracted at 11, 15, 19, 23, 31, 39, 47, 55, and 63
pixels, so footprint comparisons contain no reduction-to-reduction sampling
difference. All eight KL modes were analyzed. Mode 200 controlled the
preregistered decision among sizes through 31 pixels; larger sizes diagnosed
nonlocal response tails.

## Integrity and decision gates

| Check | Result | Diagnostic |
| --- | --- | --- |
| External 11-pixel response versus archive | Pass | Minimum cosine 0.999999; maximum projection error $4.40\times10^{-5}$ |
| Candidate size through 31 pixels | None | No size clears both edge thresholds at every primary mode-200 radius |
| First diagnostic size clearing the edge thresholds | 47 pixels | Complete support and a pass at all six radii in all eight modes |

The archived replay demonstrates that this experiment measures the same local
response as the exact-response archive. The preregistered footprint rule
requires complete support, at most one percent median squared energy on the
outer border in every radial bin, and no selected site above five percent.

The result also confirms why the larger diagnostic extractions were needed.
Border energy is not monotonic with footprint size: a 39-pixel stamp reaches a
structured response ring at radius 10 and fails there. The 47-pixel stamp
contains that ring and passes. Mode 300 has an isolated 19-pixel pass, but a
common response footprint cannot be selected from that local minimum; mode 200
is the fixed primary endpoint and 47 pixels is the first size that passes in
every mode.

## Mode-200 footprint

Energy values are median fractions of the energy in the 63-pixel diagnostic
extraction. Border values are median / maximum across the 12 sites.

| Radius | Energy in 31 | Energy in 39 | Energy in 47 | Border at 39 | Border at 47 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 7.5 | 0.9300 | 0.9645 | 0.9785 | 0.0049 / 0.0078 | 0.0036 / 0.0047 |
| 10 | 0.9043 | 0.9577 | 0.9852 | 0.0156 / 0.0346 | 0.0023 / 0.0041 |
| 12 | 0.9340 | 0.9468 | 0.9822 | 0.0064 / 0.0123 | 0.0053 / 0.0114 |
| 16 | 0.9757 | 0.9787 | 0.9825 | 0.0008 / 0.0011 | 0.0010 / 0.0057 |
| 20 | 0.9915 | 0.9937 | 0.9955 | 0.0005 / 0.0007 | 0.0005 / 0.0005 |
| 24 | 0.9945 | 0.9963 | 0.9975 | 0.0003 / 0.0004 | 0.0003 / 0.0005 |

The 31-to-47-pixel square shell contains median 4.79%, 7.87%, and 4.53% of
the 63-pixel response energy at radii 7.5, 10, and 12, respectively. Individual
sites reach 10.06%, 11.18%, and 6.79%. These tails are large enough to affect a
matched filter if they are physical. The remaining energy outside 47 pixels is
much smaller: median 2.15%, 1.48%, and 1.78% at the same radii.

## All-mode 47-pixel check

| Mode | Worst 39-pixel median border | Worst 47-pixel median border | Worst 47-pixel individual border | Radius of 47-pixel worst median |
| ---: | ---: | ---: | ---: | ---: |
| 125 | 0.0162 | 0.0054 | 0.0116 | 12 |
| 150 | 0.0165 | 0.0052 | 0.0114 | 12 |
| 175 | 0.0155 | 0.0052 | 0.0114 | 12 |
| 200 | 0.0156 | 0.0053 | 0.0114 | 12 |
| 225 | 0.0160 | 0.0053 | 0.0115 | 12 |
| 250 | 0.0163 | 0.0056 | 0.0112 | 12 |
| 300 | 0.0160 | 0.0057 | 0.0113 | 12 |
| 350 | 0.0163 | 0.0051 | 0.0107 | 12 |

All 47-pixel selected-site stamps have complete finite support. The table's
39-pixel worst medians occur at radius 10; the 47-pixel worst medians occur at
radius 12.

## Interpretation and next checkpoint

No local footprint through 31 pixels resolves the Stage-A edge trigger. A
47-pixel footprint does, but it is about 13 $\lambda/D$ wide and was designated
as a nonlocal-tail diagnostic rather than an automatic promotion candidate.
The next checkpoint therefore repeats six angularly distributed sites at each
of radii 7.5, 10, and 12 with half and twice the original perturbation. It tests
whether the 31-to-47-pixel shell has stable cosine and amplitude across
contrast. This short check separates a repeatable KLIP response tail from a
finite-difference floor before the roughly 16.6-hour full-field response is
recomputed at 47 pixels. See the
[contrast-linearity setup](../klip-response-tail-linearity-setup-20260921/README.md).

## ROC receipt

The canonical run is
`working/roc/klip_response_stamp_convergence_20260921`. Its complete receipt
records:

- `results.json`: `a3d631ab07c1806c6337ec88e33d5e35d78bc68808673922607e7162686e1178`;
- `results.md`: `9e5b0c2a01ea6cd4f98be52374321a1cbd1578cdbba6f81f58d71bf259fb6af1`.

The receipt status is complete and fingerprints all 144 reduction receipts.
