

addpath(genpath('/data/fniess/matlab/'));

slopes=readNPY('mrsi_3D_slopes_map.npy');
pvalue=readNPY('mrsi_3D_pvalue_map.npy');
rvalue=readNPY('mrsi_3D_rvalue_map.npy');
intercept=readNPY('mrsi_3D_intercept_map.npy');
stderr=readNPY('mrsi_3D_stderr_map.npy');
tau=readNPY('mrsi_3D_tau.npy');
tau_sd=readNPY('mrsi_3D_tau_sd.npy');

niftiwrite(slopes,'slopes.nii');
niftiwrite(pvalue,'pvalue.nii');
niftiwrite(rvalue,'rvalue.nii');
niftiwrite(intercept,'intercept.nii');
niftiwrite(stderr,'stderr.nii');
niftiwrite(tau,'tau.nii');
niftiwrite(tau_sd,'tau_sd.nii');


