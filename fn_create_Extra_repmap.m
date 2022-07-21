function fn_create_FWHM_repmap(nreps)

addpath(genpath('/Users/fabian/matlab/'));

for x=1:nreps
    
    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Orig/gl_u_o+gl_n_o_sd_map.nii'));
    Glx_o_sd_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    Glx_o_sd_map_flipped(isnan(Glx_o_sd_map_flipped))=0;
   % Glx_o_sd_map_flipped(Glx_o_sd_map_flipped>19.99)=0;

    Glx_o_sd_map(:,:,:,x)=temp.img;
    Glx_o_sd_map(isnan(Glx_o_sd_map))=0;
 %  Glx_o_sd_map(Glx_o_sd_map>19.99)=0;

    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Orig/gl_u_n+gl_n_n_sd_map.nii'));
    Glx_n_sd_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    Glx_n_sd_map_flipped(isnan(Glx_n_sd_map_flipped))=0;
 %  Glx_n_sd_map_flipped(Glx_n_sd_map_flipped>19.99)=0;
 
    Glx_n_sd_map(:,:,:,x)=temp.img;
    Glx_n_sd_map(isnan(Glx_n_sd_map))=0;
 %  Glx_o_sd_map(Glx_o_sd_map>19.99)=0;


    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Orig/Cr+PCr_sd_map.nii'));
    tCr_sd_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    tCr_sd_map_flipped(isnan(tCr_sd_map_flipped))=0;
 %  tCr_sd_map_flipped(tCr_sd_map_flipped>19.99)=0;
 
    tCr_sd_map(:,:,:,x)=temp.img;
    tCr_sd_map(isnan(tCr_sd_map))=0;
 %  Glx_o_sd_map(Glx_o_sd_map>19.99)=0;

    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Orig/NAA+NAAG_sd_map.nii'));
    tNAA_sd_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    tNAA_sd_map_flipped(isnan(tNAA_sd_map_flipped))=0;
 %  tNAA_sd_map_flipped(tNAA_sd_map_flipped>19.99)=0;
 
    tNAA_sd_map(:,:,:,x)=temp.img;
    tNAA_sd_map(isnan(tNAA_sd_map))=0;
 %  tNAA_sd_map(tNAA_sd_map>19.99)=0;

    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Orig/gl_cao+gl_cbo_sd_map.nii'));
    glc6_o_sd_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    glc6_o_sd_map_flipped(isnan(glc6_o_sd_map_flipped))=0;
 
 
    glc6_o_sd_map(:,:,:,x)=temp.img;
    glc6_o_sd_map(isnan(glc6_o_sd_map))=0;

    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Orig/gl_can+gl_cbn_sd_map.nii'));
    glc6_n_sd_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    glc6_n_sd_map_flipped(isnan(glc6_n_sd_map_flipped))=0;
 
 
    glc6_n_sd_map(:,:,:,x)=temp.img;
    glc6_n_sd_map(isnan(glc6_n_sd_map))=0;


    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Orig/gl_cbo_sd_map.nii'));
    glcb_o_sd_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    glcb_o_sd_map_flipped(isnan(glcb_o_sd_map_flipped))=0;
 
 
    glcb_o_sd_map(:,:,:,x)=temp.img;
    glcb_o_sd_map(isnan(glcb_o_sd_map))=0;

    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Orig/gl_cbn_sd_map.nii'));
    glcb_n_sd_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    glcb_n_sd_map_flipped(isnan(glcb_n_sd_map_flipped))=0;
 
 
    glcb_n_sd_map(:,:,:,x)=temp.img;
    glcb_n_sd_map(isnan(glcb_n_sd_map))=0;


    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Extra/FWHM_map.nii'));
    FWHM_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    FWHM_map_flipped(isnan(FWHM_map_flipped))=0;
   
    FWHM_map(:,:,:,x)=temp.img;
    FWHM_map(isnan(FWHM_map))=0;
    
    
    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Extra/SNR_PseudoReplica_spectral_map.nii'));
    SNR_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    SNR_map_flipped(isnan(SNR_map_flipped))=0;
    
    SNR_map(:,:,:,x)=temp.img;
    SNR_map(isnan(SNR_map))=0;
    
    temp=load_untouch_nii(strcat('maps/',int2str(x),'/Extra/0_pha_map.nii'));
    zero_pha_map_flipped(:,:,:,x)=flipud(fliplr(temp.img));
    zero_pha_map_flipped(isnan(zero_pha_map_flipped))=0;
    
    zero_pha_map(:,:,:,x)=temp.img;
    zero_pha_map(isnan(SNR_map))=0;





    
end

save(['FWHM_map_flipped.mat'],'FWHM_map_flipped')
save(['SNR_map_flipped.mat'],'SNR_map_flipped')
save(['zero_pha_map_flipped.mat'],'zero_pha_map_flipped')
save(['Glx_o_sd_map_flipped.mat'],'Glx_o_sd_map_flipped')
save(['Glx_n_sd_map_flipped.mat'],'Glx_n_sd_map_flipped')

save(['FWHM_map.mat'],'FWHM_map')
save(['SNR_map.mat'],'SNR_map')
save(['zero_pha_map.mat'],'zero_pha_map')
save(['Glx_o_sd_map.mat'],'Glx_o_sd_map')
save(['Glx_n_sd_map.mat'],'Glx_n_sd_map')
save(['Glc6_o_sd_map.mat'],'glc6_o_sd_map')
save(['Glc6_n_sd_map.mat'],'glc6_n_sd_map')


writeNPY(FWHM_map,'FWHM.npy');
writeNPY(SNR_map,'SNR.npy');
writeNPY(zero_pha_map,'zero_pha.npy');
writeNPY(Glx_o_sd_map,'Glx_o_sd_map.npy');
writeNPY(Glx_n_sd_map,'Glx_n_sd_map.npy');
writeNPY(glc6_o_sd_map,'Glc6_o_sd_map.npy');
writeNPY(glc6_n_sd_map,'Glc6_n_sd_map.npy');
writeNPY(glcb_o_sd_map,'Glcb_o_sd_map.npy');
writeNPY(glcb_n_sd_map,'Glcb_n_sd_map.npy');
writeNPY(tCr_sd_map,'tCr_sd_map.npy');
writeNPY(tNAA_sd_map,'tNAA_sd_map.npy');
end
