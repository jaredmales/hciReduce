#!/usr/bin/env python3
"""Verify extension of an 11-pixel Welch PSD to larger KLIP templates."""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import cho_solve
from scipy.sparse.linalg import LinearOperator, cg

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_p4_step5_welch_psd as p4_psd  # noqa: E402
import run_klip_covariance_stage_a as stage  # noqa: E402


TRAINING_SIZE = 11
RESPONSE_SIZES = (11, 31, 47)


def estimation_window(name: str) -> np.ndarray:
    """Return one fixed 11-pixel PSD estimation window."""
    stage.require(name in ("rectangular", "hann"), "unknown PSD estimation window")
    if name == "rectangular":
        return np.ones((TRAINING_SIZE, TRAINING_SIZE), dtype=np.float64)
    return np.outer(np.hanning(TRAINING_SIZE), np.hanning(TRAINING_SIZE))


def fit_extended_psd(samples: np.ndarray, response_size: int,
                     window_name: str, mixing: float) -> dict[str, object]:
    """Fit an 11-pixel periodogram on the linear-lag grid of a response support."""
    samples = np.asarray(samples, dtype=np.float64)
    stage.require(samples.ndim == 2 and samples.shape[1] == TRAINING_SIZE * TRAINING_SIZE and
                  len(samples) >= 8 and np.all(np.isfinite(samples)), "invalid PSD training patches")
    stage.require(response_size in RESPONSE_SIZES and 0 < mixing <= 1, "invalid PSD extension policy")
    mean = np.mean(samples, axis=0)
    centered = samples - mean
    target = float(np.sum(np.square(centered)) /
                   ((len(samples) - 1) * TRAINING_SIZE * TRAINING_SIZE))
    stage.require(np.isfinite(target) and target > 0, "invalid unwindowed training variance")
    window = estimation_window(window_name)
    window_energy = float(np.sum(np.square(window)))
    fft_size = 2 * response_size - 1
    transformed = np.fft.fft2(centered.reshape((-1, TRAINING_SIZE, TRAINING_SIZE)) * window,
                              s=(fft_size, fft_size), axes=(-2, -1))
    raw_power = np.sum(np.square(np.abs(transformed)), axis=0) / ((len(samples) - 1) * window_energy)
    raw_zero_lag = float(np.mean(raw_power))
    stage.require(np.isfinite(raw_zero_lag) and raw_zero_lag > 0, "invalid windowed PSD variance")
    rescaling = target / raw_zero_lag
    power = (1 - mixing) * raw_power * rescaling + mixing * target
    stage.require(np.all(np.isfinite(power)) and np.min(power) >= mixing * target * (1 - 1e-12),
                  "extended PSD lost its positive spectral floor")
    lag = np.fft.ifft2(power)
    stage.require(np.max(np.abs(lag.imag)) <= 1e-12 * target, "extended PSD produced a complex lag kernel")
    return {"mean": mean, "target_variance": target, "window": window_name,
            "window_energy": window_energy, "mixing": mixing, "fft_size": fft_size,
            "raw_windowed_zero_lag": raw_zero_lag, "psd_rescaling": rescaling,
            "power": power, "lag": lag.real}


def dense_covariance(model: dict[str, object], response_size: int) -> np.ndarray:
    """Construct the finite block-Toeplitz covariance for verification."""
    stage.require(int(model["fft_size"]) == 2 * response_size - 1, "PSD grid and response support differ")
    yy, xx = np.indices((response_size, response_size))
    delta_y = yy.ravel()[:, None] - yy.ravel()[None, :]
    delta_x = xx.ravel()[:, None] - xx.ravel()[None, :]
    size = int(model["fft_size"])
    covariance = np.asarray(model["lag"])[delta_y % size, delta_x % size]
    return 0.5 * (covariance + covariance.T)


def covariance_operator(model: dict[str, object], response_size: int) -> LinearOperator:
    """Return an FFT convolution operator for the finite covariance."""
    size = int(model["fft_size"])
    stage.require(size == 2 * response_size - 1, "PSD grid and response support differ")
    power = np.asarray(model["power"], dtype=np.float64)
    pixels = response_size * response_size

    def matvec(vector: np.ndarray) -> np.ndarray:
        padded = np.zeros((size, size), dtype=np.float64)
        padded[:response_size, :response_size] = np.asarray(vector).reshape((response_size, response_size))
        filtered = np.fft.ifft2(np.fft.fft2(padded) * power).real
        return filtered[:response_size, :response_size].ravel()

    return LinearOperator((pixels, pixels), matvec=matvec, rmatvec=matvec, dtype=np.float64)


def solve_template(model: dict[str, object], template: np.ndarray,
                   relative_tolerance: float = 1e-10) -> dict[str, object]:
    """Solve the finite PSD covariance and normalize weights to unit response."""
    template = np.asarray(template, dtype=np.float64)
    stage.require(template.ndim == 2 and template.shape[0] == template.shape[1] and
                  template.shape[0] in RESPONSE_SIZES and np.all(np.isfinite(template)),
                  "invalid response template")
    support = template.shape[0]
    vector = template.ravel()
    operator = covariance_operator(model, support)
    target = float(model["target_variance"])
    preconditioner = LinearOperator(operator.shape, matvec=lambda value: np.asarray(value) / target,
                                    rmatvec=lambda value: np.asarray(value) / target, dtype=np.float64)
    iterations = 0

    def count_iteration(_: np.ndarray) -> None:
        nonlocal iterations
        iterations += 1

    inverse_template, information = cg(operator, vector, M=preconditioner, rtol=relative_tolerance,
                                       atol=0, maxiter=4 * vector.size, callback=count_iteration)
    residual = operator @ inverse_template - vector
    relative_residual = float(np.linalg.norm(residual) / np.linalg.norm(vector))
    stage.require(information == 0 and relative_residual <= 5 * relative_tolerance,
                  "finite PSD solve did not converge")
    energy = float(vector @ inverse_template)
    stage.require(np.isfinite(energy) and energy > 0, "finite PSD solve has nonpositive template energy")
    return {"weight": inverse_template / energy, "inverse_template": inverse_template,
            "energy": energy, "sigma": 1 / math.sqrt(energy),
            "iterations": iterations, "relative_residual": relative_residual}


def check() -> None:
    """Verify legacy equivalence, positivity, lag support, and iterative solves."""
    generator = np.random.default_rng(67381)
    samples = generator.normal(size=(64, TRAINING_SIZE * TRAINING_SIZE))
    samples += 0.35 * np.roll(samples, 1, axis=1)
    template_11 = generator.normal(size=(TRAINING_SIZE, TRAINING_SIZE))
    checks = 0
    for window_name in ("rectangular", "hann"):
        for mixing in (0.1, 0.3, 1.0):
            model = fit_extended_psd(samples, 11, window_name, mixing)
            legacy = p4_psd.fit_psd(samples, window_name, mixing)
            stage.require(legacy is not None and np.array_equal(model["mean"], legacy["mean"]) and
                          np.isclose(model["target_variance"], legacy["target_variance"],
                                     rtol=1e-15, atol=0) and
                          np.allclose(model["power"], legacy["power"], rtol=2e-15, atol=1e-15),
                          "extended 11-pixel PSD differs from the tested P4 estimator")
            covariance = dense_covariance(model, 11)
            stage.require(np.allclose(covariance, legacy["covariance"], rtol=2e-14,
                                      atol=2e-14 * model["target_variance"]),
                          "extended 11-pixel covariance differs from the tested P4 estimator")
            solved = solve_template(model, template_11)
            legacy_inverse = cho_solve(legacy["factorization"], template_11.ravel(), check_finite=False)
            legacy_energy = float(template_11.ravel() @ legacy_inverse)
            stage.require(np.allclose(solved["weight"], legacy_inverse / legacy_energy,
                                      rtol=2e-9, atol=2e-11),
                          "iterative 11-pixel weights differ from the dense tested solution")
            checks += 1

    diagnostics = {}
    for response_size in RESPONSE_SIZES:
        model = fit_extended_psd(samples, response_size, "rectangular", 0.3)
        lag = np.asarray(model["lag"])
        indices = np.arange(int(model["fft_size"]))
        signed = np.where(indices <= int(model["fft_size"]) // 2, indices,
                          indices - int(model["fft_size"]))
        outside = ((np.abs(signed[:, None]) > TRAINING_SIZE - 1) |
                   (np.abs(signed[None, :]) > TRAINING_SIZE - 1))
        if np.any(outside):
            stage.require(np.max(np.abs(lag[outside])) <= 2e-14 * model["target_variance"],
                          "11-pixel periodogram created unsupported long-lag covariance")
        template = generator.normal(size=(response_size, response_size))
        operator = covariance_operator(model, response_size)
        first = generator.normal(size=template.size)
        second = generator.normal(size=template.size)
        stage.require(np.isclose(first @ (operator @ second), second @ (operator @ first),
                                 rtol=2e-13, atol=2e-11 * model["target_variance"]),
                      "finite PSD operator is not symmetric")
        stage.require(first @ (operator @ first) > 0, "finite PSD operator is not positive")
        solved = solve_template(model, template)
        isotropic = fit_extended_psd(samples, response_size, "hann", 1.0)
        isotropic_solve = solve_template(isotropic, template)
        expected_weight = template.ravel() / float(np.sum(np.square(template)))
        expected_sigma = math.sqrt(float(isotropic["target_variance"]) /
                                   float(np.sum(np.square(template))))
        stage.require(np.allclose(isotropic_solve["weight"], expected_weight,
                                  rtol=2e-11, atol=2e-13) and
                      np.isclose(isotropic_solve["sigma"], expected_sigma, rtol=2e-11),
                      "extended isotropic endpoint differs from identity weighting")
        diagnostics[str(response_size)] = {"iterations": solved["iterations"],
                                           "relative_residual": solved["relative_residual"]}
        checks += 4
    print(f"KLIP PSD extension checks passed: {checks} contracts; {diagnostics}", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the numerical-contract command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("action", choices=("check",))
    return result


def main() -> None:
    """Run the requested numerical contract check."""
    parser().parse_args()
    check()


if __name__ == "__main__":
    main()
