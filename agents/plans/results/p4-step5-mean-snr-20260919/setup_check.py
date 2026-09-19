"""Verify optimized sampling and mean-toggle algebra against the original study."""
from pathlib import Path
import sys
import numpy as np
from astropy.io import fits
sys.path.insert(0,str(Path('agents/plans/scripts').resolve()))
import compare_p4_step5_mean_snr as comparison
import run_p4_step5_roc_full as full
from run_p4_step5_full_injections import fingerprint,write_json

root=Path('working/roc/p4_psd_full_20260918')
protocol=full.read(root/'protocol.json')
science=fits.getdata(root/'reductions/baseline/finim.fits').squeeze().astype(float)
templates=full.load_templates(root/'payload/response')
checks={'training_matrices':0,'protected_stamps':0,'unit_increments':0,'generic_psd_solves':0}
for site in protocol['sites']:
    mask=full.holdout_mask(science.shape,protocol,site)
    positions=[(site['row']+dx,site['column']+dy) for dx,dy in full.radial.OFFSETS]
    positions += [(t['row'],t['column']) for t in (protocol['calibration_trials'][0],protocol['calibration_trials'][-1])]
    for position in positions:
        narrow=comparison.narrow_samples(science,position,mask)
        rings=full.radial.geometry(science,position,mask)
        for offset,matrix in narrow.items():
            assert np.array_equal(matrix,full.radial.extract(science,rings[offset]))
            checks['training_matrices']+=1
    x,y=site['row'],site['column'];template=templates[(x,y)]
    before,diagnostic=comparison.filter_pixel(science,(x,y),template,mask)
    changed=science.copy();changed[y-5:y+6,x-5:x+6]+=2e-5*template.reshape(11,11)
    matrices=comparison.narrow_samples(changed,(x,y),mask)
    original=comparison.narrow_samples(science,(x,y),mask)
    for offset in matrices:assert np.array_equal(matrices[offset],original[offset])
    after,_=comparison.filter_pixel(changed,(x,y),template,mask)
    for name in before:
        assert np.isclose(after[name]-before[name],2e-5,rtol=1e-10,atol=1e-13)
        checks['unit_increments']+=1
    checks['protected_stamps']+=1
    models=full.candidate_models(science,(x,y),mask)
    data=science[y-5:y+6,x-5:x+6].ravel()
    for prefix,_,_,_,old in comparison.FAMILIES:
        model=models[old];q=np.linalg.solve(model['covariance'],template);weight=q/(template@q)
        assert np.isclose(weight@data,before[prefix+'_psd_no_mean'],rtol=1e-10,atol=1e-13)
        assert np.isclose(weight@(data-model['mean']),before[prefix+'_psd_mean'],rtol=1e-10,atol=1e-13)
        checks['generic_psd_solves']+=1
write_json(Path(__file__).with_name('setup_checks.json'),{'all_passed':True,'checks':checks,
    'script':fingerprint(Path(comparison.__file__).resolve()),'source':fingerprint(root/'reductions/baseline/finim.fits'),
    'new_reductions':0})
print(checks)
