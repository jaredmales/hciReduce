Tools for reducing high contrast imaging observations.

Use `python3 scripts/contrast_curve.py --help` to run repeated KLIP or P4 planet
injections and plot contrast versus separation. The [contrast-curve guide](doc/contrast_curve.dox)
describes configuration overrides, radius-dependent contrasts, and resuming runs.

# Building

```sh
cmake -S . -B _build
cmake --build _build -j
```

## Tests and documentation

Build tests on demand with `cmake --build _build --target tests`, or configure
them for CTest with `-DHCIREDUCE_BUILD_TESTS=ON`. Generate the Doxygen site with
`cmake --build _build --target docs`.

The `coverage` target generates an instrumented test build and embeds its lcov
report in the Doxygen site when `lcov` and `genhtml` are installed.
