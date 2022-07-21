%%% GO TO Subject output part 1 folder
load('CombinedCSI.mat')
csi_combined=csi;
image_combined=image_new;
weights_combined=weights;

%%% GO TO Subject output part 2 folder
load('CombinedCSI.mat')

csi_combined(:,:,:,:,8:14)=csi;
image_combined(:,:,:,:,:,8:14)=image_new;
weights_combined(:,:,:,:,:,8:14)=weights;

dwelltime_c=ReadInInfo.Par.Dwelltime*10e-10;
image_new_cc = squeeze(sum(image_combined .* repmat(conj(weights_combined),[1,1,1,1,size(image_combined,5),1]),1));
b0_map = squeeze(median(angle(image_new_cc(:,:,:,5+1:end,:)) - angle(image_new_cc(:,:,:,5:end-1,:)),4)./(2*pi*dwelltime_c)) ;

mask_rep=repmat(mask,[1 1 1 size(b0_map,4)]);
b0_map=b0_map.*mask_rep;

N_fid = size(csi_combined,4);
time_vec = linspace(0,(N_fid-1)*dwelltime_c,N_fid);
time_vec=reshape(time_vec,[1 1 1 size(time_vec,2) 1]);
time_vec_matrix=repmat(time_vec,[size(csi_combined,1) size(csi_combined,2) size(csi_combined,3) 1 size(csi_combined,5)]);

b0_map=reshape(b0_map,[size(b0_map,1) size(b0_map,2) size(b0_map,3) 1 size(b0_map,4)]);
b0_map=repmat(b0_map,[1 1 1 size(csi_combined,4) 1]);

%%% GO TO Subject output combined folder

%if .npy files do not exist yet
%GM=load_untouch_nii('GM_mask.nii');
%GM=GM.img;
%WM=load_untouch_nii('WM_mask.nii');
%WM=WM.img;

GM=readNPY('GM_mask.npy');
WM=readNPY('WM_mask.npy');

GM_spectra=flipud(fliplr(GM));
WM_spectra=flipud(fliplr(WM));

load('FWHM_map_flipped.mat')
load('SNR_map_flipped.mat')

% Stricter QC for better results
QC=(SNR_map_flipped>25) & (FWHM_map_flipped<0.05);

%alread done
%GM_spectra(GM_spectra<0.8)=0;
%GM_spectra(GM_spectra>0.79)=1;
%WM_spectra(WM_spectra<0.8)=0;
%WM_spectra(WM_spectra>0.79)=1;

GM_spectra=repmat(GM_spectra,[1 1 1 14]);
WM_spectra=repmat(WM_spectra,[1 1 1 14]);

GM_spectra=GM_spectra.*QC;
WM_spectra=WM_spectra.*QC;

GM_rep=reshape(GM_spectra,[size(GM_spectra,1) size(GM_spectra,2) size(GM_spectra,3) 1 size(GM_spectra,4)]);
WM_rep=reshape(WM_spectra,[size(WM_spectra,1) size(WM_spectra,2) size(WM_spectra,3) 1 size(WM_spectra,4)]);

GM_rep=repmat(GM_rep,[1 1 1 1176 1]);
WM_rep=repmat(WM_rep,[1 1 1 1176 1]);


FWHM_spectra=reshape(FWHM_map_flipped,[32 32 21 1 14]);
FWHM_spectra=repmat(FWHM_spectra,[1 1 1 1176 1]);

mask_rep=reshape(mask_rep,[32 32 21 1 14]);
mask_rep=repmat(mask_rep,[1 1 1 1176 1]);

exp_factor=exp(2*pi*i*b0_map.*time_vec_matrix-time_vec_matrix.*FWHM_spectra).*mask_rep;
exp_factor_onlyb0=exp(2*pi*i*b0_map.*time_vec_matrix).*mask_rep;

%b0_csi_combined=csi_combined.*exp_factor_onlyb0;
b0_csi_combined=csi_combined.*exp_factor;

GM_csi=b0_csi_combined.*GM_rep;
WM_csi=b0_csi_combined.*WM_rep;

u=1;
GM_first_point=zeros(1,1176);
GM_last_point=zeros(1,1176);
WM_first_point=zeros(1,1176);
WM_last_point=zeros(1,1176);

for x=1:32
for y=1:32
for z=1:21
if (sum(GM_csi(x,y,z,:,1))~=0)
GM_first_point(u,:)=GM_csi(x,y,z,:,1);
u=u+1;
end
end
end
end
u=1;

for x=1:32
for y=1:32
for z=1:21
if (sum(GM_csi(x,y,z,:,14))~=0)
GM_last_point(u,:)=GM_csi(x,y,z,:,14);
u=u+1;
end
end
end
end
u=1;

for x=1:32
for y=1:32
for z=1:21
if (sum(WM_csi(x,y,z,:,1))~=0)
WM_first_point(u,:)=WM_csi(x,y,z,:,1);
u=u+1;
end
end
end
end
u=1;

for x=1:32
for y=1:32
for z=1:21
if (sum(WM_csi(x,y,z,:,14))~=0)
WM_last_point(u,:)=WM_csi(x,y,z,:,14);
u=u+1;
end
end
end
end
% Save .mat files for later use in python
save(['GM_last_point.mat'],'GM_last_point');
save(['WM_last_point.mat'],'WM_last_point');
save(['WM_first_point.mat'],'WM_first_point');
save(['GM_first_point.mat'],'GM_first_point');

% GO TO PYTHON and use fn_fmrsi_mat2rda.py %Dont forget template .dat file
% for header
