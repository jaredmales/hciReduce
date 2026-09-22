# KLIP response-tail contrast-linearity result

## Question and design

The response-stamp convergence experiment found that 47 pixels is the first
footprint with complete selected-site support that clears the original edge
thresholds in every KL mode. At the three inner radii, the 31-to-47 square
shell contains enough energy to affect a matched filter. This experiment tested
whether that outer structure is a contrast-linear KLIP response or a
finite-difference floor before starting a full 47-pixel response campaign.

The completed response-stamp convergence run was the immutable parent. At
radii 7.5, 10, and 12 pixels, the experiment selected alternating even-index
sites 0, 2, 4, 6, 8, and 10 from the parent's geometry-only angular selection.
The parent supplied the central difference at contrast scale 1. Positive and
negative perturbations at scales 0.5 and 2 supplied 72 new reductions. All
three derivatives used the same known-planet subtraction and all eight KL
modes.

The analysis split each 63-pixel response into the 11-pixel core, the 11-to-31
shell, the controlling 31-to-47 shell, and the residual 47-to-63 shell. It also
formed the zero-contrast Richardson estimate

$$
D(0) \simeq \frac{4D(\epsilon/2)-D(\epsilon)}{3}.
$$

## Promotion gates

| Check | Result |
| --- | --- |
| Mode-200 31-to-47 shell cosine | Pass |
| Mode-200 31-to-47 shell projection | Pass |
| Richardson-extrapolated 47-pixel edge | Pass |
| Promote 47-pixel response | Yes |

The fixed tail gate required median cosine at least 0.95 and median projection
within 0.10 of unity at both new contrast scales and every inner radius. The
Richardson response also had to retain complete support, at most one percent
median border energy, and at most five percent at every site.

## Mode-200 tail

Values are medians across six angularly distributed sites. Projection is the
candidate derivative onto the scale-1 derivative. Tail energy is the scale-1
31-to-47 shell as a fraction of 63-pixel response energy.

| Radius | Scale | Tail cosine | Tail projection | Scaled residual | Tail energy fraction |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 7.5 | 0.5 | 0.9676 | 1.0189 | 0.2526 | 0.0479 |
| 7.5 | 2 | 0.9876 | 0.9463 | 0.1567 | 0.0479 |
| 10 | 0.5 | 0.9955 | 1.0096 | 0.0951 | 0.0803 |
| 10 | 2 | 0.9969 | 0.9715 | 0.0777 | 0.0803 |
| 12 | 0.5 | 0.9948 | 1.0098 | 0.1018 | 0.0462 |
| 12 | 2 | 0.9982 | 0.9632 | 0.0603 | 0.0462 |

The weakest agreement is the half-contrast radius-7.5 shell, but it remains
above the preregistered cosine threshold and close to unit projection. The
larger scaled residual there reflects shape variation within a shell containing
only 4.79% of the total response energy; it does not indicate an amplitude
collapse.

## All-mode diagnostics

Each row takes the worst median over radii 7.5, 10, and 12 and contrast scales
0.5 and 2. Edge values use the Richardson 47-pixel response.

| Mode | Minimum tail cosine | Maximum projection error | Maximum median border energy | Maximum individual border energy |
| ---: | ---: | ---: | ---: | ---: |
| 125 | 0.9766 | 0.0364 | 0.0057 | 0.0108 |
| 150 | 0.9793 | 0.0358 | 0.0054 | 0.0097 |
| 175 | 0.9814 | 0.0350 | 0.0052 | 0.0096 |
| 200 | 0.9676 | 0.0537 | 0.0053 | 0.0145 |
| 225 | 0.9781 | 0.0371 | 0.0053 | 0.0126 |
| 250 | 0.9852 | 0.0461 | 0.0059 | 0.0105 |
| 300 | 0.9781 | 0.0367 | 0.0053 | 0.0103 |
| 350 | 0.9673 | 0.0470 | 0.0052 | 0.0100 |

Every Richardson stamp has complete support. The outer response therefore
persists across perturbation amplitude and KL mode, rather than behaving like a
single-contrast numerical floor.

## Initial-run metadata repair

The first run stopped after two completed tasks because the six-significant-
digit `FAKECONT` header rounded the doubled contrast to `0.00914873`, only
$7.17\times10^{-12}$ beyond the original absolute-tolerance boundary. The
configured reduction and product were correct. The repair changed validation
to a one-part-per-million relative comparison, fingerprinted the old and new
runners, validated and receipted the existing product, and resumed without
rerunning it. One recovered task consequently lacks a wall-time entry; the
other 71 reductions used 5.92 summed task-minutes with a 4.99-second median.

## Interpretation and next checkpoint

The 47-pixel footprint is promoted. The next
[exact-response campaign](../klip-response-47-setup-20260921/README.md) measures
all 11,192 integer search locations with native 47-by-47 response stamps. Its
acceptance checks require a bitwise-identical signal-free baseline, unchanged
coordinates, agreement of every central 11-pixel response with the archived
exact field, and agreement of complete 47-pixel responses with the independent
fixed-site derivatives from this experiment.

## ROC receipt

The canonical run is
`working/roc/klip_response_tail_linearity_20260921`. Its completion receipt
records:

- `results.json`: `83e0a74c86e7be67d97339757d3b19acfdf957c15d3f5c22bb6e69b2748a5987`;
- `results.md`: `cfa323a842aab1ca10af0c9cd99ce82273ae2e2af845a07c07a0348631bdf0e5`.

The receipt status is complete and fingerprints all 72 reduction receipts.
