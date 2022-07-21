#!/bin/bash


bet2 magnitude.nii m_magnitude -f 0.33 -g 0
fast -t 1 -o masks -n 3 m_magnitude.nii.gz
gzip -d -f masks_pve_1.nii.gz
gzip -d -f masks_pve_2.nii.gz
nii2mnc masks_pve_1.nii masks_pve_1.mnc
nii2mnc masks_pve_2.nii masks_pve_2.mnc
mincresample -like ../mask.mnc masks_pve_1.mnc GM_mask.mnc -nearest_neighbour
mincresample -like ../mask.mnc masks_pve_2.mnc WM_mask.mnc -nearest_neighbour
mnc2nii GM_mask.mnc GM_mask.nii
mnc2nii WM_mask.mnc WM_mask.nii
