function fn_3Dnifti2text(niftifile,outputfile)

data_in=load_untouch_nii(niftifile);
data_in=data_in.img;

fid=fopen([outputfile],'wt');
for x=1:size(data_in,1)
for y=1:size(data_in,2)
for z=1:size(data_in,3)
fprintf(fid,'%d  %d  %d  %d  \n',x-1,y-1,z-1,data_in(x,y,z));
end 
end 
end

fclose(fid);

end
