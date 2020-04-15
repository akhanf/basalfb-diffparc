import nibabel as nib
import numpy as np
mask_nib = nib.load(snakemake.input.mask)
mask_vol = mask_nib.get_fdata()
mask_indices = mask_vol > 0 
masked = mask_vol[mask_indices]

nvoxels = len(masked)
ntargets = len(snakemake.input.connmap_3d)
conn = np.zeros((nvoxels,ntargets))
for i,conn_file in enumerate(snakemake.input.connmap_3d):
    vol = nib.load(conn_file).get_fdata()
    masked =  vol[mask_indices].T
    conn[:,i] = masked
np.savez(snakemake.output.connmap_npz, conn=conn,mask=mask_vol,affine=mask_nib.affine)


