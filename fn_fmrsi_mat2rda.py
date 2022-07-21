#!/Users/fabian/opt/miniconda3/envs/py27/bin/python

import numpy as np
import fn_twix_parser_multi_raid_VE12 as fn_tp
import fn_spect as spect
import scipy.io as spio

scans=fn_tp.read('meas_MID00457_FID71746_fn_31P_staticDRESS_16meas_120V.dat');

header=spect.read_rda_header('meas_MID00457_FID71746_fn_31P_staticDRESS_16meas_120V.dat',False,200);
spect.rda_header_modify('VectorSize:',1176,header);
spect.rda_header_modify('MRFrequency:',123.252263,header);
spect.rda_header_modify('Nucleus:','1H',header);
spect.rda_header_modify('DwellTime:',754.8,header);
spect.rda_header_modify('MagneticFieldStrength:',2.894,header);


mat=spio.loadmat('GM_first_point.mat');
GM_first_point=mat['GM_first_point'];
mat=spio.loadmat('GM_last_point.mat');
GM_last_point=mat['GM_last_point'];
mat=spio.loadmat('WM_first_point.mat');
WM_first_point=mat['WM_first_point'];
mat=spio.loadmat('WM_last_point.mat');
WM_last_point=mat['WM_last_point'];



spect.rda_header_modify('CSIMatrixSize[0]:',np.shape(GM_first_point)[0],header);
spect.write_rda(GM_first_point, header, 'GM_first.rda')
spect.rda_header_modify('CSIMatrixSize[0]:',np.shape(GM_last_point)[0],header);
spect.write_rda(GM_last_point, header, 'GM_last.rda')

spect.rda_header_modify('CSIMatrixSize[0]:',np.shape(WM_first_point)[0],header);
spect.write_rda(WM_first_point, header, 'WM_first.rda')
spect.rda_header_modify('CSIMatrixSize[0]:',np.shape(WM_last_point)[0],header);
spect.write_rda(WM_last_point, header, 'WM_last.rda')