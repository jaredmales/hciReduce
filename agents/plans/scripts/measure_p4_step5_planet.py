#!/usr/bin/env python3
"""Measure AF Lep b with fixed Step-5 filters and the user's hciAnalyze settings.

Use the completed ROC baseline without another reduction. Export filtered
amplitudes for production annular SNR, retaining conditional scores separately.
Invalid inner-radius covariance fits never fall back to another estimator.
"""
from __future__ import annotations

import argparse
import configparser
from pathlib import Path
import shutil
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
from astropy.io import fits

import run_p4_step5_roc_full as full
from run_p4_step5_full_injections import fingerprint, write_json

LABELS = {'unfiltered': 'Unfiltered image', 'gaussian': 'Gaussian FWHM 3.6', 'identity': 'Identity matched filter',
          'psd_hann_b0_m0.1': 'Same-radius Hann PSD', 'psd_rectangular_b5_m0.3': 'Pooled rectangular PSD',
          'pca_b5_f1': 'Pooled PCA, floor 1', 'isotropic_b5': 'Pooled isotropic'}


def annular_oracle(amplitude: np.ndarray, settings: dict) -> tuple[np.ndarray, list]:
    """Reconstruct production one-pixel annular means/stddevs and the small-sample multiplier."""
    yy, xx = np.indices(amplitude.shape)
    radius = np.hypot(xx-127.5, yy-127.5).astype('f4')
    source_radius = np.hypot(xx-settings['source_x'], yy-settings['source_y']).astype('f4')
    noise = np.isfinite(amplitude) & (source_radius > settings['source_radius']+.5)
    profile = []
    for lower in range(181):
        # The image center is half-integral: no native radius lies exactly on an integer boundary.
        values = amplitude[noise & (radius > lower) & (radius <= lower+1)].astype(float)
        profile.append({'radius': lower+.5, 'pixels': len(values),
            'mean': float(np.mean(values)) if len(values) else None,
            'stddev': float(np.std(values, ddof=1)) if len(values) > 1 else None})
    mean = np.interp(radius, [p['radius'] for p in profile],
                     [p['mean'] if p['mean'] is not None else np.nan for p in profile])
    sigma = np.interp(radius, [p['radius'] for p in profile],
                      [p['stddev'] if p['stddev'] is not None else np.nan for p in profile]).astype('f4')
    with np.errstate(divide='ignore', invalid='ignore'):
        score = ((amplitude-mean)/sigma).astype('f4')
        count = 2*np.pi*np.hypot(xx-127.5, yy-127.5)/settings['lambda_d']-1
        correction = np.where(count > 0, 1/np.sqrt(1+1/count), 0).astype('f4')
        score *= correction
    score[(radius < settings['min_radius']) | (radius > settings['max_radius']) | ~np.isfinite(amplitude)] = np.nan
    return score, profile


def run(args: argparse.Namespace) -> None:
    """Freeze inputs, build fixed filter maps, and measure the known planet with production SNR."""
    root, study, software = args.output.resolve(), args.study.resolve(), args.software.resolve()
    source = study/'reductions/baseline/finim.fits'
    science, header = fits.getdata(source, header=True)
    science = science.squeeze().astype(float)
    protocol = full.read(study/'protocol.json')
    config = configparser.ConfigParser()
    config.read_string('[root]\n'+args.config.read_text())
    full.radial.require(float(config['planet']['sep']) == float(header['PLANETSEP']) and
                        float(config['planet']['PA']) == float(header['PLANETPA']), 'planet coordinates changed')
    settings = {'lambda_d': float(config['root']['lambdaD']), 'source_x': protocol['known_source_circle'][0],
        'source_y': protocol['known_source_circle'][1], 'source_radius': float(config['planet']['R']),
        'aperture_radius': float(config['snr']['apertureR']), 'pixel_buffer': .5,
        'min_radius': float(config['snr']['minRad']), 'max_radius': 60., 'nearest_pixel': [139, 126],
        'models': list(LABELS), 'minimum_training_patches': 8,
        'covariance_holdout': 'original 28 calibration footprints plus every 15x15 Gaussian footprint in the configured planet aperture and the source circle',
        'annular_noise_mask': 'finite filtered pixels outside planet.R plus the production half-pixel buffer; original calibration pixels may enter annular normalization',
        'statistic': 'production annular SNR of filtered amplitude; conditional scores are separate diagnostics',
        'parent_study': str(study)}
    full.radial.require(settings['lambda_d'] == 3.6 and settings['aperture_radius'] == 3 and settings['source_radius'] == 7 and
                        settings['min_radius'] == 6 and 'maxRad' not in config['snr'], 'changed requested guidance')
    software_records = full.read(args.software.parent/'manifest.json')['software']
    full.verify(software_records)
    full.radial.require((software/'hciAnalyze.cpp').read_bytes() == Path('src/apps/hciAnalyze.cpp').read_bytes(),
                        'frozen analysis source differs from production')
    inputs = [fingerprint(source), fingerprint(args.config.resolve()), fingerprint(study/'protocol.json'),
              *[fingerprint(p) for p in sorted((study/'payload/response').glob('*.fits'))],
              fingerprint(study/'references/baseline/gaussian.fits'),
              *[fingerprint(Path(__file__).with_name(name).resolve()) for name in full.SCRIPTS], fingerprint(Path(__file__).resolve())]
    full.verify(full.read(study/'reductions/baseline/complete.json')['products'])
    root.mkdir(parents=True, exist_ok=False)
    (root/'maps').mkdir()
    (root/'analysis').mkdir()
    shutil.copy2(args.config, root/'analyze.conf')
    yy, xx = np.indices(science.shape)
    distance = np.hypot(xx-settings['source_x'], yy-settings['source_y']).astype('f4')
    aperture = distance <= settings['aperture_radius']+.5
    forbidden = distance <= settings['source_radius']+.5
    for trial in protocol['calibration_trials']:
        full.mark_footprint(forbidden, trial['row'], trial['column'])
    for y, x in np.argwhere(aperture):
        forbidden[y-7:y+8, x-7:x+8] = True
    fits.writeto(root/'training_exclusion.fits', forbidden.astype('u1'), header)
    settings['aperture_pixels'] = int(aperture.sum())
    write_json(root/'settings.json', settings)
    write_json(root/'manifest.json', {'inputs': inputs, 'software': software_records,
        'settings': fingerprint(root/'settings.json'), 'config': fingerprint(root/'analyze.conf'),
        'new_reductions': 0, 'post_study_known_planet_diagnostic': True})
    templates = full.load_templates(study/'payload/response')
    names = (*full.COVARIANCE_MODELS, 'identity')
    amplitude = {name: np.full(science.shape, np.nan) for name in names}
    sigma = {name: np.full(science.shape, np.nan) for name in names}
    scores = {name: np.full(science.shape, np.nan) for name in names}
    training = {name: np.zeros(science.shape, dtype='i4') for name in full.COVARIANCE_MODELS}
    pixels, checked = [], 0
    for index, ((x, y), template) in enumerate(sorted(templates.items())):
        data = science[y-5:y+6, x-5:x+6].ravel()
        if data.size != 121 or not np.isfinite(data).all():
            continue
        models = full.candidate_models(science, (x, y), forbidden)
        row = {'x': x, 'y': y, 'models': {}}
        rings = full.radial.geometry(science, (x, y), forbidden) if aperture[y, x] else None
        for name, model in models.items():
            if model is None:
                row['models'][name] = {'valid': False, 'reason': 'fewer than eight accepted training patches',
                    'training_patches': len(rings[0]['halves']) if rings and name.startswith('psd_hann') else None}
                continue
            value = full.radial.filter_stamp(data, template, model, np.ones(121))
            amplitude[name][y, x], sigma[name][y, x], scores[name][y, x] = value['amplitude'], value['sigma'], value['score']
            training[name][y, x] = model['samples']
            row['models'][name] = {'valid': True, 'training_patches': model['samples'],
                                  **{k: value[k] for k in ('amplitude', 'sigma', 'score')}}
            if aperture[y, x]:
                weight = np.linalg.solve(model['covariance'], template)
                energy = float(template @ weight)
                full.radial.require(np.isclose(float(weight @ (data-model['mean'])/energy), value['amplitude'], rtol=1e-9, atol=1e-12)
                                    and np.isclose(1/np.sqrt(energy), value['sigma'], rtol=1e-10), 'generic solve mismatch')
                checked += 1
        energy = float(template @ template)
        amplitude['identity'][y, x], sigma['identity'][y, x] = float(template @ data/energy), 1/np.sqrt(energy)
        scores['identity'][y, x] = amplitude['identity'][y, x]/sigma['identity'][y, x]
        row['models']['identity'] = {'valid': True, 'amplitude': amplitude['identity'][y, x],
                                     'sigma': sigma['identity'][y, x], 'score': scores['identity'][y, x]}
        if aperture[y, x]:
            pixels.append(row)
        if index % 1000 == 0:
            print(f'filtered {index}/{len(templates)} native positions', flush=True)
    full.radial.require(len(pixels) == int(aperture.sum()), 'missing aperture responses')
    for name in names:
        fits.writeto(root/'maps'/(name+'.fits'), amplitude[name].astype('f4'), header)
        fits.writeto(root/'maps'/(name+'_conditional.fits'), scores[name].astype('f4'), header)
        fits.writeto(root/'maps'/(name+'_sigma.fits'), sigma[name].astype('f4'), header)
        if name in training:
            fits.writeto(root/'maps'/(name+'_samples.fits'), training[name], header)
    shutil.copy2(source, root/'maps/unfiltered.fits')
    shutil.copy2(study/'references/baseline/gaussian.fits', root/'maps/gaussian.fits')
    environment = full.binary_environment(study, reduction=False)
    environment['LD_LIBRARY_PATH'] = str(software)
    summaries, profiles, annular_maps = {}, {}, {}
    for name in LABELS:
        image_path = root/'maps'/(name+'.fits')
        command = ['taskset', '-c', '12,13', str(software/'hciAnalyze'), '--config', str(root/'analyze.conf'),
            '--file='+str(image_path), '--filter.psfResponse=', '--filter.lpfGaussFW=0', '--filter.hpfGaussFW=0',
            '--noise.model=identity', '--noise.outputDiagnostics=false', '--noise.only=false', '--diagnostics=true']
        write_json(root/'analysis'/(name+'_command.json'), command)
        with (root/'analysis'/(name+'.log')).open('w') as log:
            subprocess.run(command, env=environment, cwd=root/'analysis', stdout=log, stderr=subprocess.STDOUT, check=True)
        snr_path = image_path.with_name(name+'_snr.fits')
        actual, actual_header = fits.getdata(snr_path, header=True)
        actual = actual.squeeze()
        native = fits.getdata(image_path).squeeze()
        full.radial.require(actual_header['SNRAPER'] == 3 and actual_header['SNRMINR'] == 6 and
                            actual_header['SNRMAXR'] == 60 and actual_header['SNRSMALL'] == 1 and actual_header['SNRMEAN'] == 1,
                            'production SNR settings mismatch')
        expected, profiles[name] = annular_oracle(native, settings)
        valid = aperture & np.isfinite(native)
        full.radial.require(np.isfinite(expected[valid]).all() and np.allclose(actual[valid], expected[valid], rtol=2e-6, atol=2e-6),
                            'independent annular normalization mismatch: '+name)
        error = float(np.max(np.abs(actual[valid]-expected[valid])))
        yx = np.argwhere(valid)
        peak = int(np.argmax(actual[valid])); py, px = yx[peak]
        full.radial.require(f'{float(actual[py, px]):.6g}' == (root/'analysis'/(name+'.log')).read_text().splitlines()[-1].split()[-1],
                            'CLI aperture maximum mismatch')
        actual = np.where(np.isfinite(native), actual, np.nan)
        annular_maps[name] = actual
        center = np.isfinite(native[126, 139])
        summaries[name] = {'label': LABELS[name], 'center_valid': bool(center), 'valid_aperture_pixels': int(valid.sum()),
            'total_aperture_pixels': int(aperture.sum()), 'annular_snr_center': float(actual[126,139]) if center else None,
            'annular_snr_aperture_maximum': float(actual[py,px]), 'peak_x_y': [int(px),int(py)],
            'annular_oracle_max_absolute_difference': error,
            'conditional_center': next(p['models'][name] for p in pixels if (p['x'],p['y']) == (139,126)) if name in names else None}
    # Confirm exported Gaussian and identity amplitudes agree with the production filtering paths.
    replay_errors = {}
    for name, response, fwhm in (('gaussian', '', '3.6'), ('identity', str(study/'payload/response/p4PSF_manifest.fits'), '0')):
        directory = root/'analysis'/('direct_'+name)
        directory.mkdir()
        shutil.copy2(source, directory/'science.fits')
        command = ['taskset','-c','12,13',str(software/'hciAnalyze'),'--config',str(root/'analyze.conf'),
                   '--file='+str(directory/'science.fits'),'--filter.psfResponse='+response,'--filter.lpfGaussFW='+fwhm,
                   '--filter.hpfGaussFW=0','--noise.model=identity']
        with (directory/'run.log').open('w') as log:
            subprocess.run(command, env=environment, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)
        direct = fits.getdata(directory/'science_snr.fits').squeeze()
        replay = fits.getdata(root/'maps'/(name+'_snr.fits')).squeeze()
        full.radial.require(np.allclose(direct, replay, rtol=3e-6, atol=2e-5, equal_nan=True), 'direct production filter mismatch')
        replay_errors[name] = float(np.nanmax(np.abs(direct-replay)))
    write_json(root/'results.json', {'settings': settings, 'models': summaries, 'aperture_pixels': pixels,
        'generic_covariance_solves_verified': checked, 'production_filter_replay_max_absolute_errors': replay_errors,
        'study_reference_aliases': {'gaussian_raw': 'gaussian', 'gaussian_snr': 'gaussian', 'identity_snr': 'identity'},
        'caveats': ['Same-radius Hann is invalid at the planet center; its partial-aperture peak is diagnostic only.',
                    'Conditional scores and annular SNRs use different noise estimates; C=I conditional score has no fitted noise scale.',
                    'The known planet lies inside the radii tested by the 90-injection study; outer-site thresholds are not applied.']})
    write_json(root/'annular_profiles.json', profiles)
    fig, axes = plt.subplots(2, 4, figsize=(13, 7), layout='constrained')
    for axis, (name, values) in zip(axes.flat, annular_maps.items()):
        shown = axis.imshow(values[116:137,129:150], origin='lower', extent=(128.5,149.5,115.5,136.5), vmin=-6, vmax=6, cmap='RdBu_r')
        axis.add_patch(Circle((settings['source_x'],settings['source_y']),3.5,fill=False,color='black',linestyle='--'))
        axis.plot(settings['source_x'],settings['source_y'],'k+',ms=8)
        axis.set_title(LABELS[name]+f'\nSNR peak {summaries[name]["annular_snr_aperture_maximum"]:.2f}; {summaries[name]["valid_aperture_pixels"]}/39 pixels',fontsize=9)
        axis.set(xlabel='Native x',ylabel='Native y')
    axes.flat[-1].axis('off')
    axes.flat[-1].text(0,.85,'AF Lep b, original ROC baseline\n\nConfigured aperture radius: 3 px\nProduction pixel buffer: 0.5 px\nSource exclusion: 7 + 0.5 px\nλ/D: 3.6 px\n\nWhite pixels: unsupported\nHann center: invalid',transform=axes.flat[-1].transAxes,va='top',fontsize=10)
    fig.colorbar(shown,ax=list(axes.flat[:-1]),label='Application annular SNR',shrink=.8)
    fig.savefig(root/'comparison.png',dpi=170)
    plt.close(fig)
    full.verify([*inputs,*software_records])
    write_json(root/'complete.json', {'new_reductions': 0, 'filters': len(LABELS), 'inputs_unchanged': True,
        'products': [fingerprint(p) for p in sorted(root.rglob('*')) if p.is_file()]})
    print(summaries,flush=True)


def main() -> None:
    """Measure the known planet without changing the completed injection study."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--config', type=Path, default=Path('working/analyze.conf'))
    parser.add_argument('--software', type=Path, default=Path('working/roc/p4_noise_step5_gaussian_20260918/software'))
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
