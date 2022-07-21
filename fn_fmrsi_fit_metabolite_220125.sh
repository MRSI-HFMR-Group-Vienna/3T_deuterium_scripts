#!/bin/bash

BD=$PWD 


echo "Which metabolite?"
read m
echo "How many time points?"
read t
echo "Matrix Size?"
read r
read p
read z
echo "Scan Time?"
read time

mkdir -p fmrsi_fit/$m


for ((i=1;i<$t+1;i++)); do
	m_path="'maps/"$i"/Orig/"$m"_amp_map.nii'";
	#m_path="'maps/"$i"/QualityAndOutlier_Clip/"$m"_amp_map.nii'";
	#m_path="'maps/"$i"/Ratio/"$m"_RatToNAA+NAAG_map.nii'";
	#m_path="'maps/"$i"/Ratio/"$m"_RatToCr+PCr_map.nii'";
	o_path="'fmrsi_fit/"$m"/"$m"_"$i".txt'";
	 /Applications/MATLAB_R2021b.app/bin/matlab -nodisplay -r 'addpath(genpath("/Users/fabian/matlab/")); fn_3Dnifti2text('$m_path','$o_path'); exit; '
done

cd fmrsi_fit/$m
fn_2H_mrsi_220124.py -t $t -m $m -r $r -p $p -z $z -x $time


 /Applications/MATLAB_R2021b.app/bin/matlab -nodisplay < /Users/fabian/matlab/fn_numpy2nifti.m

cd ../../
