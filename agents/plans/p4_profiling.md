# P4 Profiling Plan

## Goal

Provide trustworthy, comparable P4 timing measurements for the reduction lifecycle and the per-annulus local-regression
work. Replace KLIP's public timing globals with a reusable, caller-owned timing facility that P4 and KLIP can both use.
The facility must distinguish elapsed wall time from aggregate worker CPU time and must not affect science values, the
chosen PCA modes, or the scheduling policy.

## Current state

- `HCIobservation` already records wall-clock boundaries for loading, fake injection, coaddition, preprocessing,
  derotation, combination, and final output.
- `KLIPreduction` adds public `t_worker_*`, `t_eigenv`, `t_klim`, and `t_psf` globals plus `dump_times()`. The latter
  three values accumulate work across OpenMP workers but are not represented as a reusable structured result.
- P4 records only `m_geometrySeconds` and `m_regressionSeconds`; enabled diagnostics write these two values in
  `p4Timing.fits`. This is insufficient to attribute time between sampling, PCA setup/solve, projection, and output.

## Design

1. Introduce a small shared timing type owned by a reduction instance, rather than process-global variables. It should
   hold named elapsed wall-clock intervals and named aggregate worker intervals, support reset/snapshot/report, and
   document that worker totals can exceed elapsed time with parallel execution.
2. Keep lifecycle measurements in the existing `HCIobservation` timing boundaries, but expose them through the shared
   report rather than having KLIP print directly from public member data.
3. Define algorithm-stage names consistently:
   - geometry and predictor/target sampling;
   - normal-equation/Gram construction;
   - eigensolve and rank selection;
   - modal projection and residual application;
   - diagnostic writing, derotation, combination, and final-product writing.
4. Use monotonic `mx::sys::get_curr_time()` boundaries for elapsed stages. Inside OpenMP loops, accumulate independent
   per-thread stage totals with explicit OpenMP reductions or thread-local result objects merged after the loop; never
   mutate a shared `double` from workers.
5. Preserve the existing P4 diagnostic contract by extending `p4Timing.fits` with named rows/metadata (or an
   explicitly versioned table) rather than silently changing the meaning of its two existing values. Include elapsed
   and aggregate-worker units in seconds.
6. Make human-readable timing output opt-in or tied to the existing diagnostics/progress policy, so normal reductions
   do not gain unsolicited output. The structured snapshot remains available to applications and tests.

## Implementation sequence

1. Inventory every existing timing boundary in `HCIobservation`, `ADIobservation`, KLIP, and P4; establish the exact
   ownership/lifetime of the common timing object and migrate the KLIP counters without changing reported values.
2. Add P4 boundaries around geometry, per-pixel sampling, P4PCA Gram construction, eigensolve/rank selection, and
   projection/residual application. Ensure exception paths still close or clearly mark incomplete intervals.
3. Add a versioned P4 timing diagnostic and update the P4 configuration/user documentation to describe its fields.
4. Add deterministic unit tests using a clock seam or injected timing recorder. Tests must verify reset between
   reductions, correct aggregation of multiple worker records, no data races under OpenMP, and diagnostic schema/value
   placement. Do not assert real elapsed durations.
5. Run a representative AF Lep configuration across controlled OpenMP and BLAS thread counts. Report per-stage elapsed
   and worker totals, CPU utilization, peak RSS, and variance over repeated runs before choosing any optimization or
   GPU work.

## Acceptance criteria

- P4 and KLIP use the same documented timing representation; neither exposes mutable public timing globals.
- Reported wall and worker times remain finite, reset per reduction, and are race-free with more than one OpenMP thread.
- P4 diagnostics identify every recorded timing value and retain compatibility or an explicit schema version.
- Focused timing tests and existing P4/KLIP regression tests pass without changing reduction output.
- A reproducible benchmark record identifies the dominant P4 stages before further performance work is authorized.

## Risks and decisions deferred

- The measurement clock must remain suitable for elapsed time; CPU time is an aggregate accounting metric, not a
  replacement for wall time.
- Fine-grained timers can distort very small local fits. Establish stage granularity from baseline profiling and measure
  instrumentation overhead before retaining per-fit timers.
- Do not add CUDA or alter OpenMP/BLAS thread limits as part of timing instrumentation. Those are follow-on decisions
  driven by the benchmark record.
