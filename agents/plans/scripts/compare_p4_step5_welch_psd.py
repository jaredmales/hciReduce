#!/usr/bin/env python3
"""Test a fixed-frequency, Welch-style covariance on the saved Step-5 noise splits.

Average windowed patch periodograms on a zero-padded grid, reconstruct a finite
stationary covariance, and solve on untapered candidate coordinates. Compare the
prespecified windows/floors with three-mode and isotropic controls; no reductions
or injection recovery are performed by this development diagnostic.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.linalg import cho_factor, eigvalsh

import compare_p4_step5_radial_pooling as radial
import compare_p4_step5_shrinkage as shrinkage
import diagnose_p4_step5_projected_noise as projection
from compare_p4_step5_variance_floor import compare_saved
from run_p4_step5_full_injections import fingerprint, write_json

VARIANTS = [('rectangular', 0.1), ('rectangular', 0.3), ('hann', 0.1), ('hann', 0.3)]


def fit_psd(samples: np.ndarray, window_name: str, mixing: float) -> dict | None:
    """Estimate tapered power, preserve unwindowed variance, and build a positive finite lag covariance."""
    radial.require(window_name in ('rectangular', 'hann'), 'unknown estimation window')
    radial.require(np.isfinite(mixing) and 0 < mixing <= 1, 'spectral mixing must be in (0, 1]')
    radial.require(samples.ndim == 2 and samples.shape[1] == 121 and np.isfinite(samples).all(), 'expected finite 11x11 patches')
    if len(samples) < 8:
        return None
    mean = samples.mean(axis=0)
    centered = samples - mean
    target = float(np.sum(centered * centered) / ((len(samples)-1) * 121))
    if not np.isfinite(target) or target <= 0:
        return None
    window = np.ones((11,11)) if window_name == 'rectangular' else np.outer(np.hanning(11), np.hanning(11))
    energy = float(np.sum(window * window))
    # Padding to 2L-1 retains every linear lag without aliasing opposite stamp edges.
    transformed = np.fft.fft2(centered.reshape((-1,11,11)) * window, s=(21,21), axes=(-2,-1))
    power = np.sum(np.abs(transformed)**2, axis=0) / ((len(samples)-1) * energy)
    raw_zero_lag = float(power.mean())
    if not np.isfinite(raw_zero_lag) or raw_zero_lag <= 0:
        return None
    rescaling = target / raw_zero_lag
    power = power * rescaling
    regularized = (1-mixing) * power + mixing * target
    autocovariance = np.fft.ifft2(regularized)
    radial.require(np.max(np.abs(autocovariance.imag)) <= 1e-12 * target, 'nonreal autocovariance')
    yy, xx = np.indices((11,11))
    dy, dx = yy.ravel()[:,None]-yy.ravel()[None,:], xx.ravel()[:,None]-xx.ravel()[None,:]
    covariance = autocovariance.real[dy % 21, dx % 21]
    radial.require(np.allclose(covariance, covariance.T, rtol=1e-12, atol=1e-12*target), 'asymmetric lag covariance')
    covariance = (covariance + covariance.T) / 2
    radial.require(np.isclose(np.trace(covariance), 121*target, rtol=1e-12, atol=0), 'PSD covariance changes total variance')
    eigenvalues = eigvalsh(covariance, check_finite=False)
    radial.require(eigenvalues[0] >= mixing*target*(1-1e-10), 'PSD covariance loses its positive floor')
    return {'mean': mean, 'covariance': covariance, 'factorization': cho_factor(covariance, lower=True, check_finite=False),
        'target_variance': target, 'window': window_name, 'window_energy': energy, 'mixing': mixing,
        'raw_windowed_zero_lag': raw_zero_lag, 'psd_rescaling': rescaling,
        'condition_number': float(eigenvalues[-1]/eigenvalues[0]), 'power': regularized}


def policy_name(base: str, window: str, mixing: float) -> str:
    """Name the fixed sampling, estimation-window, and spectral-mixing choices."""
    return f'{base}_{window}_m{mixing:g}'


def plots(summary: dict, output: Path) -> None:
    """Compare calibrated variance, actual noise, and stability with both fixed controls."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 12), layout='constrained')
    panels = [
        ('Training variance / prediction', lambda r: r['training']['median_variance_over_predicted'], True),
        ('Held-out variance / prediction', lambda r: r['same_radius']['median_variance_over_predicted'], True),
        ('Held-out amplitude variance / PCA control', lambda r: r['median_paired_variance_over_pca1'], True),
        ('Held-out amplitude variance / isotropic control', lambda r: r['median_paired_variance_over_isotropic'], True),
        ('Split-weight cosine: same 16 sites', lambda r: r['median_split_weight_cosine'], False),
        ('Median held-out fraction within ±1 sigma', lambda r: r['same_radius']['median_within_one_sigma_fraction'], False)]
    for axis, (title, getter, logarithmic) in zip(axes.flat, panels):
        values = np.array([[getter(summary['pca_controls'][base]), getter(summary['isotropic_controls'][base])] +
                           [getter(summary['policies'][policy_name(base,w,m)]) for w,m in VARIANTS]
                           for base,_,_ in radial.POLICIES])
        rendered = np.log10(values) if logarithmic else values
        extent = max(1, math.ceil(np.max(np.abs(rendered)))) if logarithmic else 1
        lower, upper = (-extent, extent) if logarithmic else (0, 1)
        mesh = axis.imshow(rendered, vmin=lower, vmax=upper, cmap='RdBu_r' if logarithmic else 'Blues', aspect='auto')
        axis.set(title=title, xticks=range(6), xticklabels=('PCA f=1','Iso','Rect .1','Rect .3','Hann .1','Hann .3'),
                 yticks=range(8), yticklabels=[f'{"Norm" if n else "Raw"} ±{w}' for _,w,n in radial.POLICIES])
        axis.tick_params(axis='x', labelsize=9)
        for y in range(8):
            for x in range(6):
                dark = abs(rendered[y,x]) > .6*extent if logarithmic else rendered[y,x] > .6
                axis.text(x,y,f'{values[y,x]:.2g}',ha='center',va='center',color='white' if dark else 'black',fontsize=9)
        fig.colorbar(mesh,ax=axis,shrink=.7,label='log10 ratio' if logarithmic else 'fraction / cosine')
    fig.suptitle('Step 5 fixed-frequency PSD: 21×21 zero padding, finite 11×11 lag covariance\n'
                 '16 reused sites × two directions; same held-out ring; development diagnostic')
    fig.savefig(output/'comparison.png',dpi=170)
    plt.close(fig)


def main() -> None:
    """Freeze the PSD design, reproduce both controls, and evaluate the unchanged angular splits."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shrinkage',type=Path,required=True)
    parser.add_argument('--projected-noise',type=Path,required=True)
    parser.add_argument('--floor-comparison',type=Path,required=True)
    parser.add_argument('--radial-comparison',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    prior=radial.read_json(args.shrinkage/'manifest.json')
    common=prior['common_sites']
    radial.require(len(common)==16,'changed common-site protocol')
    previous_pca={(r['site'],r['policy'],r['training_half']):r for r in radial.read_json(args.projected_noise/'projections.json')}
    previous_iso={(r['site'],r['policy'],r['training_half']):r for r in radial.read_json(args.shrinkage/'projections.json') if r['gamma']==1}
    sampling=radial.read_json(args.radial_comparison/'manifest.json')
    nulls=radial.read_json(args.floor_comparison/'nulls.json')
    trials=[r for r in nulls if r['name'] in common]
    radial.require(len(trials)==16,'missing common trial')
    inputs=[*prior['inputs'],*[fingerprint(args.shrinkage/n) for n in ('manifest.json','projections.json','summary.json')]]
    for root,name in ((args.projected_noise,'projections.json'),(args.radial_comparison,'manifest.json'),(args.floor_comparison,'nulls.json')):
        radial.require(fingerprint(root/name) in inputs,'changed archived control input')
    for record in inputs:
        radial.require(fingerprint(Path(record['path']))==record,'changed original input')
    scripts=[fingerprint(Path(__file__)),*[fingerprint(Path(__file__).with_name(n)) for n in
        ('compare_p4_step5_radial_pooling.py','compare_p4_step5_shrinkage.py','diagnose_p4_step5_projected_noise.py',
         'compare_p4_step5_variance_floor.py','analyze_p4_step5_development.py','run_p4_step5_full_injections.py')]]
    args.output.mkdir(parents=True,exist_ok=False)
    write_json(args.output/'manifest.json',{'purpose':'fixed-frequency Welch-style PSD development comparison',
        'variants':[{'window':w,'spectral_isotropic_mixing':m} for w,m in VARIANTS], 'sampling_policies':radial.POLICIES,
        'common_sites':common,'directions_per_site':2,'fft_shape':[21,21],'covariance_shape':[121,121],
        'mean':'subtract the unwindowed across-patch mean; retain the same fitted mean in all filters',
        'psd':'sum(abs(FFT21(window*(patch-mean)))**2)/((n-1)*sum(window**2)); rescale spectral mean to unwindowed mean pixel variance',
        'regularization':'(1-mixing)*rescaled_power + mixing*unwindowed_mean_pixel_variance',
        'edge_policy':'finite covariance from linear lags -10..10; no 11-pixel periodic wrap and no deconvolution of window autocorrelation',
        'candidate':'untapered data and response; only existing native radial normalization applies',
        'training':prior['training'],'primary_validation':prior['primary_validation'],'normalization':prior['normalization'],
        'controls':'three modes at floor 1 and mean-variance isotropic covariance, same fitted means and samples',
        'inputs':inputs,'scripts':scripts})
    science=fits.getdata(sampling['inputs'][0]['path']).squeeze().astype(float)
    field=Path(sampling['inputs'][1]['path']).parent
    yy,xx=np.indices(science.shape);forbidden=np.zeros(science.shape,bool)
    for x,y,radius in sampling['exclusions']:
        forbidden |= np.hypot(xx-x,yy-y)<=radius
    scale_map,profile=radial.variance_profile(science,forbidden)
    checks=compare_saved(profile,radial.read_json(args.floor_comparison/'variance_profiles.json')['baseline'],'profile')
    coordinates=fits.getdata(field/'p4PSF_coordinates.fits').T
    responses=fits.getdata(field/'p4PSF_model_0000.fits').astype(float)
    validity=fits.getdata(field/'p4PSF_validity_0000.fits').ravel()
    templates={(int(c[0]),int(c[1])):t.ravel() for c,t,v in zip(coordinates,responses,validity) if v==1 and np.isfinite(t).all()}
    records,stability,pca_rows,iso_rows,pca_stability,iso_stability=[],[],[],[],[],[]
    endpoints=0
    for trial in trials:
        x,y=trial['row'],trial['column'];rings=radial.geometry(science,(x,y),forbidden)
        supports=[np.unique(np.concatenate([projection.native_support(r,r['halves']==h) for r in rings.values()])) for h in (0,1)]
        radial.require(not np.intersect1d(*supports).size,'train/validation supports overlap')
        radial.require(not np.any(forbidden.ravel()[np.concatenate(supports)]),'forbidden native pixels used')
        template=templates[(x,y)]
        for base,width,normalized in radial.POLICIES:
            image=science/scale_map if normalized else science
            scale=scale_map[y-5:y+6,x-5:x+6].ravel() if normalized else np.ones(121)
            matrices={o:radial.extract(image,r) for o,r in rings.items()}
            selected=[o for o in rings if abs(o)<=width]
            details={'pca':[],'iso':[]}|{policy_name(base,w,m):[] for w,m in VARIANTS}
            for half in (0,1):
                training=np.vstack([matrices[o][rings[o]['halves']==half] for o in selected])
                controls={}
                for kind,model,saved,destination in (
                    ('pca',radial.fit(training,1),previous_pca[(trial['name'],base+'_f1',half)],pca_rows),
                    ('iso',shrinkage.fit_shrinkage(training,1),previous_iso[(trial['name'],base+'_g1',half)],iso_rows)):
                    radial.require(model is not None,'invalid common control fit')
                    row,detail=shrinkage.measure(training,matrices,rings,selected,model,template,scale,half)
                    checks+=compare_saved(row,{k:saved[k] for k in row},'archived '+kind)
                    row.update({'site':trial['name'],'policy':base,'training_half':half})
                    controls[kind]=row;details[kind].append(detail);destination.append(row)
                for window in ('rectangular','hann'):
                    endpoint=fit_psd(training,window,1)
                    iso_covariance=shrinkage.fit_shrinkage(training,1)['covariance']
                    radial.require(endpoint is not None and np.allclose(endpoint['covariance'],iso_covariance,
                        rtol=1e-12,atol=1e-13*endpoint['target_variance']),'PSD isotropic endpoint differs')
                    endpoints+=1
                for window,mixing in VARIANTS:
                    model=fit_psd(training,window,mixing)
                    radial.require(model is not None,'invalid PSD fit on common support')
                    row,detail=shrinkage.measure(training,matrices,rings,selected,model,template,scale,half)
                    name=policy_name(base,window,mixing)
                    row.update({'site':trial['name'],'policy':name,'training_half':half,
                        'window':window,'mixing':mixing,'condition_number':model['condition_number'],
                        'target_variance':model['target_variance'],'raw_windowed_zero_lag':model['raw_windowed_zero_lag'],
                        'psd_rescaling':model['psd_rescaling']})
                    for kind,label in (('pca','pca1'),('iso','isotropic')):
                        reference=controls[kind]
                        row['paired_sigma_over_'+label]=row['sigma']/reference['sigma']
                        for metric,short in (('variance_over_predicted','variance'),('mse_over_predicted','mse')):
                            row['paired_'+short+'_over_'+label]=row['same_radius'][metric]*row['sigma']**2/(
                                reference['same_radius'][metric]*reference['sigma']**2)
                    records.append(row);details[name].append(detail)
            pca_stability.append({'site':trial['name'],'policy':base,**shrinkage.split_stability(details['pca'])})
            iso_stability.append({'site':trial['name'],'policy':base,**shrinkage.split_stability(details['iso'])})
            for window,mixing in VARIANTS:
                name=policy_name(base,window,mixing)
                stability.append({'site':trial['name'],'policy':name,**shrinkage.split_stability(details[name])})
        print(f'completed {trial["name"]}',flush=True)
    summary={'common_sites':common,'directional_fits_per_policy':32,'pca_controls':{},'isotropic_controls':{},'policies':{}}
    for base,_,_ in radial.POLICIES:
        pca=[r for r in pca_rows if r['policy']==base];iso=[r for r in iso_rows if r['policy']==base]
        paired=[p['same_radius']['variance_over_predicted']*p['sigma']**2/(i['same_radius']['variance_over_predicted']*i['sigma']**2)
                for p,i in zip(pca,iso)]
        pc=shrinkage.aggregate(pca,[r for r in pca_stability if r['policy']==base])
        ic=shrinkage.aggregate(iso,[r for r in iso_stability if r['policy']==base])
        pc.update({'median_paired_variance_over_pca1':1.,'median_paired_variance_over_isotropic':float(np.median(paired))})
        ic.update({'median_paired_variance_over_isotropic':1.,'median_paired_variance_over_pca1':float(np.median(1/np.array(paired)))})
        summary['pca_controls'][base]=pc;summary['isotropic_controls'][base]=ic
        for window,mixing in VARIANTS:
            name=policy_name(base,window,mixing);rows=[r for r in records if r['policy']==name]
            one=shrinkage.aggregate(rows,[r for r in stability if r['policy']==name])
            for metric in ('paired_variance_over_pca1','paired_variance_over_isotropic','paired_mse_over_pca1',
                           'paired_mse_over_isotropic','paired_sigma_over_pca1','paired_sigma_over_isotropic','condition_number','psd_rescaling'):
                one['median_'+metric]=float(np.median([r[metric] for r in rows]))
            one['psd_rescaling_range']=[min(r['psd_rescaling'] for r in rows),max(r['psd_rescaling'] for r in rows)]
            summary['policies'][name]=one
    write_json(args.output/'projections.json',records);write_json(args.output/'stability.json',stability)
    write_json(args.output/'summary.json',summary)
    for record in [*inputs,*scripts]:
        radial.require(fingerprint(Path(record['path']))==record,'input or script changed during PSD experiment')
    verification={'archived_control_scalar_values':checks,'psd_directional_fits':len(records),'control_directional_fits':len(pca_rows)+len(iso_rows),
        'split_comparisons':len(stability),'windowed_isotropic_endpoint_checks':endpoints,'native_support_disjoint_and_excluded':True,
        'spectral_floor_and_trace_preserved':True,'physical_unit_and_standardized_solutions_agree':True,
        'input_and_script_fingerprints_unchanged':True}
    write_json(args.output/'verification.json',verification)
    plots(summary,args.output)
    write_json(args.output/'complete.json',{'baseline_images':1,'new_reductions':0,'new_positive_images':0,
        'policies':32,'control_policies':16,'common_sites':16,'directional_fits':len(records),'verification':verification})


if __name__=='__main__':
    main()
