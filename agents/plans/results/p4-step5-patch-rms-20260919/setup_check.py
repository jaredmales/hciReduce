#!/usr/bin/env python3
"""Check post-mean RMS powers against direct linear lags and physical filter responses."""
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO/'agents/plans/scripts'))
import compare_p4_step5_patch_rms as test
from run_p4_step5_full_injections import fingerprint, write_json


def direct_covariance(samples, window_name, mixing):
    """Sum every linear pixel lag without Fourier transforms or estimator helpers."""
    residual = samples-samples.mean(axis=0)
    residual /= np.sqrt(np.mean(residual**2, axis=1))[:, None]
    window = np.ones((11, 11)) if window_name == 'rectangular' else np.outer(np.hanning(11), np.hanning(11))
    patches = residual.reshape((-1, 11, 11))*window
    lags = np.empty((21, 21))
    for dy in range(-10, 11):
        for dx in range(-10, 11):
            a = patches[:, max(dy, 0):min(11, 11+dy), max(dx, 0):min(11, 11+dx)]
            b = patches[:, max(-dy, 0):min(11, 11-dy), max(-dx, 0):min(11, 11-dx)]
            lags[dy+10, dx+10] = np.sum(a*b)/((len(samples)-1)*np.sum(window**2))
    target = np.sum((samples-samples.mean(axis=0))**2)/((len(samples)-1)*121)
    lags *= target/lags[10, 10]
    yy, xx = np.indices((11, 11))
    covariance = lags[yy.ravel()[:, None]-yy.ravel()[None, :]+10, xx.ravel()[:, None]-xx.ravel()[None, :]+10]
    return (1-mixing)*covariance+mixing*target*np.eye(121)


def main():
    """Verify real-data covariance, protected source response, and synthetic limits."""
    root = REPO/'working/roc/p4_psd_full_20260918'
    science = fits.getdata(root/'reductions/baseline/finim.fits').squeeze().astype(float)
    protocol = test.full.read(root/'protocol.json')
    templates = test.full.load_templates(root/'payload/response')
    counts = {'direct_lag_covariances': 0, 'generic_filter_solves': 0, 'protected_training_matrices': 0,
              'unit_response_increments': 0, 'synthetic_identities': 0}
    max_difference = 0.
    for site in protocol['sites']:
        x, y = site['row'], site['column']
        template = templates[(x, y)]
        forbidden = test.full.holdout_mask(science.shape, protocol, site)
        matrices = test.previous.narrow_samples(science, (x, y), forbidden)
        before, diagnostics = test.filter_pixel(science, (x, y), template, forbidden)
        for prefix, window, mixing, width, _ in test.FAMILIES:
            samples = matrices[0] if width == 0 else np.vstack(list(matrices.values()))
            fitted = test.fit_patch_rms(samples, window, mixing)
            direct = direct_covariance(samples, window, mixing)
            assert np.allclose(direct, fitted['covariance'], rtol=1e-10, atol=1e-12*fitted['target_variance'])
            max_difference = max(max_difference, float(np.max(abs(direct-fitted['covariance']))/fitted['target_variance']))
            q = np.linalg.solve(direct, template)
            w = q/(template @ q)
            data = science[y-5:y+6, x-5:x+6].ravel()
            assert np.isclose(w @ data, before[prefix+'_patch_rms_no_mean'], rtol=1e-10, atol=1e-13)
            assert np.isclose(w @ (data-samples.mean(axis=0)), before[prefix+'_patch_rms_mean'], rtol=1e-10, atol=1e-13)
            counts['generic_filter_solves'] += 1
            counts['direct_lag_covariances'] += 1
        changed = science.copy()
        changed[y-5:y+6, x-5:x+6] += 2e-5*template.reshape(11, 11)
        after_matrices = test.previous.narrow_samples(changed, (x, y), forbidden)
        for offset in matrices:
            assert np.array_equal(matrices[offset], after_matrices[offset])
            counts['protected_training_matrices'] += 1
        after, _ = test.filter_pixel(changed, (x, y), template, forbidden)
        for name in before:
            assert np.isclose(after[name]-before[name], 2e-5, rtol=1e-10, atol=1e-13)
            counts['unit_response_increments'] += 1

    rng = np.random.default_rng(28317)
    base = rng.normal(size=(8, 121))
    base /= np.sqrt(np.mean(base**2, axis=1))[:, None]
    paired = np.vstack((base, -base))
    scales = np.geomspace(.1, 10, len(base))[:, None]
    unequal = np.vstack((base*scales, -base*scales))
    for window in ('rectangular', 'hann'):
        equal = test.fit_patch_rms(paired, window, .3)
        old = test.full.psd.fit_psd(paired, window, .3)
        assert np.allclose(equal['covariance'], old['covariance'], rtol=1e-10, atol=1e-12)
        varied = test.fit_patch_rms(unequal, window, .3)
        assert np.allclose(equal['covariance']/equal['target_variance'], varied['covariance']/varied['target_variance'], rtol=1e-10, atol=1e-12)
        assert not np.allclose(test.full.psd.fit_psd(unequal, window, .3)['covariance']/varied['target_variance'], equal['covariance']/equal['target_variance'], rtol=1e-3, atol=1e-3)
        endpoint = test.fit_patch_rms(unequal, window, 1)
        assert np.allclose(endpoint['covariance'], np.eye(121)*endpoint['target_variance'], rtol=1e-12, atol=1e-12)
        scaled = test.fit_patch_rms(3*unequal, window, .3)
        assert np.allclose(scaled['covariance'], 9*varied['covariance'], rtol=1e-10, atol=1e-12)
        assert test.fit_patch_rms(paired[:7], window, .3) is None
        assert test.fit_patch_rms(np.zeros((8, 121)), window, .3) is None
        assert test.fit_patch_rms(np.vstack((paired, np.zeros((1, 121)))), window, .3) is None
        counts['synthetic_identities'] += 8
    write_json(Path(__file__).with_name('setup_checks.json'), {'all_passed': True, 'checks': counts,
        'maximum_covariance_difference_over_variance': max_difference, 'new_reductions': 0,
        'script': fingerprint(Path(test.__file__).resolve()), 'source': fingerprint(root/'reductions/baseline/finim.fits')})
    print(counts)


if __name__ == "__main__":
    main()
