"""Reproduce the local full-study setup checks from the repository root."""
import json,sys,math
from pathlib import Path
import numpy as np
from astropy.io import fits
sys.path.insert(0,str(Path('agents/plans/scripts').resolve()))
import run_p4_step5_roc_full as full
from run_p4_step5_full_injections import fingerprint,write_json
root=Path('working/roc/p4_psd_full_20260918')
p=full.read(root/'protocol.json');g=full.read(root/'geometry.json');staged=full.read(root/'staging.json')
full.verify(staged['records'])
expected=json.loads(Path(__file__).with_name('independent_geometry.json').read_text())
assert [(t['row'],t['column']) for t in p['sites']]==[(t['row'],t['column']) for t in expected]
assert len(p['sites'])==30 and len(p['calibration_trials'])==28 and len(p['models'])==8
assert list(full.LEVELS)==[.5,.75,1.] and len(full.LEVELS)*len(p['sites'])==90
assert math.ceil((28+1)*.95)==28
science=fits.getdata(root/'payload/previous_baseline.fits').squeeze().astype(float)
templates=full.load_templates(root/'payload/response')
yy,xx=np.indices(science.shape)
checks={'independent_site_masks':0,'protected_covariance_fits':0,'candidate_search_stencils':0,'calibration_search_stencils':0,'unit_response_increments':0}
for site,audit in zip(p['sites'],g['sites']):
 independent=np.zeros(science.shape,bool)
 sx,sy,r=p['known_source_circle'];independent|=(xx-sx)**2+(yy-sy)**2<=r*r
 for t in [*p['calibration_trials'],site]:
  x,y=t['row'],t['column']
  for ox,oy in full.radial.OFFSETS:
   for dx in range(-7,8):
    for dy in range(-7,8):independent[y+oy+dy,x+ox+dx]=True
 mask=full.holdout_mask(science.shape,p,site)
 assert np.array_equal(independent,mask)
 checks['independent_site_masks']+=1
 checks['candidate_search_stencils']+=len(audit['candidate'])
 assert all(c['same_radius']>=8 and c['band5']>=c['same_radius'] for c in audit['candidate'])
 for records in audit['calibration'].values():
  checks['calibration_search_stencils']+=len(records)
  assert all(c['same_radius']>=8 and c['band5']>=c['same_radius'] for c in records)
 changed=science.copy();source_mask=np.zeros_like(mask)
 full.mark_footprint(source_mask,site['row'],site['column']);changed[source_mask]=np.nan
 for t in [site,p['calibration_trials'][0],p['calibration_trials'][-1]]:
  position=(t['row'],t['column'])
  before=full.candidate_models(science,position,mask);after=full.candidate_models(changed,position,mask)
  for name in full.COVARIANCE_MODELS:
   assert before[name] is not None and after[name] is not None
   assert np.array_equal(before[name]['mean'],after[name]['mean'])
   assert np.array_equal(before[name]['covariance'],after[name]['covariance'])
   assert before[name]['samples']==after[name]['samples']
   checks['protected_covariance_fits']+=1
# Known previously evaluated site: exact template-only injection checks the new measurement integration.
trial={'name':'unit_known_site','row':150,'column':140}
mask=full.holdout_mask(science.shape,p,trial)
dummy={k:np.zeros_like(science) for k in ('gaussian_raw','gaussian_snr','identity_snr')}
changed=science.copy();t=templates[(150,140)]
changed[135:146,145:156]+=2e-5*t.reshape(11,11)
a=full.measure_trial(science,trial,templates,mask,dummy);b=full.measure_trial(changed,trial,templates,mask,dummy)
for name in (*full.COVARIANCE_MODELS,'identity'):
 assert np.isclose(b[name]['center']['amplitude']-a[name]['center']['amplitude'],2e-5,rtol=1e-10,atol=0)
 checks['unit_response_increments']+=1
# The staged injection PSF retains its original 12-square support and stored normalization.
assert fingerprint(root/'payload/injection_psf.fits')['sha256']==fingerprint(Path('working/roc/p4_noise_step5_development_20260918/injection_psf.fits'))['sha256']
full.verify(staged['records'])
write_json(root/'setup_checks.json',{'checks':checks,'all_passed':True,
 'new_site_values_masked_with_nan_leave_candidate_and_calibration_covariances_unchanged':True,
 'independent_geometry_prototype_matches_all_sites':True,
 'source_psf_and_script_hashes_unchanged':True,
 'checks_do_not_run_a_roc_reduction':True})
print(json.dumps(checks,indent=2))
