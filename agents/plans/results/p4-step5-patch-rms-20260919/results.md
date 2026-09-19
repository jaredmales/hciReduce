# Per-patch normalization after mean subtraction

| Method | 0.5× /30 | 0.75× /30 | 1× /30 | Nulls /30 |
| --- | --- | --- | --- | --- |
| gaussian | 11 | 14 | 16 | 3 |
| identity | 12 | 16 | 19 | 2 |
| hann_psd_no_mean | 12 | 16 | 19 | 1 |
| hann_psd_mean | 12 | 16 | 19 | 1 |
| hann_patch_rms_no_mean | 12 | 16 | 19 | 1 |
| hann_patch_rms_mean | 12 | 16 | 19 | 1 |
| rect_psd_no_mean | 12 | 16 | 21 | 2 |
| rect_psd_mean | 12 | 16 | 21 | 2 |
| rect_patch_rms_no_mean | 12 | 16 | 21 | 2 |
| rect_patch_rms_mean | 12 | 16 | 21 | 2 |
