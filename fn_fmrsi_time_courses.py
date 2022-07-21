#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3




import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
from scipy.optimize import curve_fit
from scipy import stats
import scikit_posthocs as sp
import os

#os.chdir('/data/fniess/nmr/7T/mv_svs/')


def exp_decrease(t,m0,tau,const):
    m = m0 * np.exp(-t/tau) + const
    return m
    
    
def exp_increase(t,m0,tau,d,const):
    m = m0 * (1-d*np.exp(-t/tau)) + const
    return m    

def exp_increaseinv(t,m0,tau,d,const):
    m = m0 * np.exp(t/tau) + const
    return m  

def fit(function,time,signals,param):
    	
    	try:
    		popt, pcov = curve_fit(function, time, signals, param)
    		perr = np.sqrt(np.diag(pcov))
    		return [popt, perr]
    	except:
    		return [[0,0,0],[0,0,0]]
    		
		
######## Read in Subjects Folder
subjects=[]
subjects.append('Participant1');
subjects.append('Participant2');
subjects.append('Participant3');
subjects.append('Participant4');
subjects.append('Participant5');
subjects.append('Participant6');
subjects.append('Participant7');



###### ALL DATA POINTS glx_o, glx_n
###### Glx
GM_glx_o_all_points=[]
WM_glx_o_all_points=[]

GM_glx_n_all_points=[]
WM_glx_n_all_points=[]

###### Glc6
GM_glc6_o_all_points=[]
WM_glc6_o_all_points=[]

GM_glc6_n_all_points=[]
WM_glc6_n_all_points=[]

##### Glu
GM_glu_o_all_points=[]
WM_glu_o_all_points=[]

GM_glu_n_all_points=[]
WM_glu_n_all_points=[]

##### Gln
GM_gln_o_all_points=[]
WM_gln_o_all_points=[]

GM_gln_n_all_points=[]
WM_gln_n_all_points=[]

##### tCr
tCr_all_points=[]
GM_tCr_all_points=[]
WM_tCr_all_points=[]

##### NAA
tNAA_all_points=[]
GM_tNAA_all_points=[]
WM_tNAA_all_points=[]

GM_tNAA_tCr_all_points=[]
WM_tNAA_tCr_all_points=[]

###### ALL FITTING DATA glx_o
GM_glx_o_all_slopes=[]
WM_glx_o_all_slopes=[]

GM_glx_o_all_intercept=[]
WM_glx_o_all_intercept=[]

GM_glx_o_all_rvalue=[]
WM_glx_o_all_rvalue=[]

GM_glx_o_all_pvalue=[]
WM_glx_o_all_pvalue=[]


###### ALL FITTING DATA glc6_o
GM_glc6_o_all_slopes=[]
WM_glc6_o_all_slopes=[]

GM_glc6_o_all_intercept=[]
WM_glc6_o_all_intercept=[]

GM_glc6_o_all_rvalue=[]
WM_glc6_o_all_rvalue=[]

GM_glc6_o_all_pvalue=[]
WM_glc6_o_all_pvalue=[]


##### ALL PLOT DATA 
GM_glx_o_all_linfit=[]
WM_glx_o_all_linfit=[]

GM_glc6_o_all_linfit=[]
WM_glc6_o_all_linfit=[]



#### GLX/tCr Ratio

GM_glx_tCr_o_all_points=[]
WM_glx_tCr_o_all_points=[]

GM_glx_tCr_n_all_points=[]
WM_glx_tCr_n_all_points=[]

#### Glc6/tCr Ratio

GM_glc6_tCr_o_all_points=[]
WM_glc6_tCr_o_all_points=[]

GM_glc6_tCr_n_all_points=[]
WM_glc6_tCr_n_all_points=[]

#### count voxels
GM_glx_o_count=[]
WM_glx_o_count=[]

GM_glx_n_count=[]
WM_glx_n_count=[]

GM_glc6_o_count=[]
WM_glc6_o_count=[]

GM_glc6_n_count=[]
WM_glc6_n_count=[]

GM_tCr_count=[]
WM_tCr_count=[]

GM_tNAA_count=[]
WM_tNAA_count=[]

for i in np.arange(0,np.size(subjects)):

	os.chdir(subjects[i]+'/output_MRSI/')	 	
	glx_o=np.load('lessomit_combined/fmrsi_fit/gl_u_o+gl_n_o/mrsi_4D.npy');
	glx_n=np.load('lessomit_combined/fmrsi_fit/gl_u_n+gl_n_n/mrsi_4D.npy');

	glc6_o=np.load('glc6_combined/fmrsi_fit/glc/mrsi_4D.npy')
	glc6_n=np.load('glc6_combined/fmrsi_fit/glc/mrsi_4D.npy')


	tCr=np.load('lessomit_combined/fmrsi_fit/Cr+PCr/mrsi_4D.npy');
	tNAA=np.load('lessomit_combined/fmrsi_fit/NAA+NAAG/mrsi_4D.npy');


	zero_phase=np.load('lessomit_combined/zero_pha.npy');
	Glx_o_sd_map=np.load('lessomit_combined/Glx_o_sd_map.npy');
	Glx_n_sd_map=np.load('lessomit_combined/Glx_n_sd_map.npy');
	Glc6_o_sd_map=np.load('glc6_combined/Glc6_o_sd_map.npy');
	Glc6_n_sd_map=np.load('glc6_combined/Glc6_n_sd_map.npy');
	tCr_sd_map=np.load('lessomit_combined/tCr_sd_map.npy');
	tNAA_sd_map=np.load('lessomit_combined/tNAA_sd_map.npy');


	FWHM=np.load('lessomit_combined/FWHM.npy');
	SNR=np.load('lessomit_combined/SNR.npy');
	GM_mask=np.load('lessomit_combined/GM_mask.npy');
	WM_mask=np.load('lessomit_combined/WM_mask.npy');
	QC=(FWHM<0.1)*(SNR>15);
	time_real=np.arange(0,60.19999999999999,4.3);

	glx_o[np.isnan(glx_o)]=0;
	glx_n[np.isnan(glx_n)]=0;
	glc6_o[np.isnan(glc6_o)]=0;
	glc6_n[np.isnan(glc6_n)]=0;
	tCr[np.isnan(tCr)]=0;
	tNAA[np.isnan(tNAA)]=0;

	m_glx_o=glx_o*QC*(Glx_o_sd_map<25);
	m_glx_n=glx_n*QC*(Glx_n_sd_map<25);
	m_glc6_o=glc6_o*QC;
	m_glc6_n=glc6_n*QC;
	
	m_tCr=tCr*QC*(tCr_sd_map<25);
	m_NAA=tNAA*QC*(tNAA_sd_map<25);


###### count QC surviving voxels

	GM_glx_o_count.append(np.mean([np.count_nonzero(m_glx_o[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_glx_o,3))]))
	WM_glx_o_count.append(np.mean([np.count_nonzero(m_glx_o[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_glx_o,3))]))

	GM_glx_n_count.append(np.mean([np.count_nonzero(m_glx_n[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_glx_n,3))]))
	WM_glx_n_count.append(np.mean([np.count_nonzero(m_glx_n[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_glx_n,3))]))

	GM_glc6_o_count.append(np.mean([np.count_nonzero(m_glc6_o[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_glc6_o,3))]))
	WM_glc6_o_count.append(np.mean([np.count_nonzero(m_glc6_o[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_glc6_o,3))]))

	GM_glc6_n_count.append(np.mean([np.count_nonzero(m_glc6_n[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_glc6_n,3))]))
	WM_glc6_n_count.append(np.mean([np.count_nonzero(m_glc6_n[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_glc6_n,3))]))

	GM_tCr_count.append(np.mean([np.count_nonzero(m_tCr[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_tCr,3))]))
	WM_tCr_count.append(np.mean([np.count_nonzero(m_tCr[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_tCr,3))]))

	GM_tNAA_count.append(np.mean([np.count_nonzero(m_NAA[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_NAA,3))]))
	WM_tNAA_count.append(np.mean([np.count_nonzero(m_NAA[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_NAA,3))]))


#### glx
	GM_glx_o=[np.sum(m_glx_o[:,:,:,x] * GM_mask)/np.count_nonzero(m_glx_o[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_glx_o,3))]
	WM_glx_o=[np.sum(m_glx_o[:,:,:,x] * WM_mask)/np.count_nonzero(m_glx_o[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_glx_o,3))]

	GM_glx_n=[np.sum(m_glx_n[:,:,:,x] * GM_mask)/np.count_nonzero(m_glx_n[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_glx_n,3))]
	WM_glx_n=[np.sum(m_glx_n[:,:,:,x] * WM_mask)/np.count_nonzero(m_glx_n[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_glx_n,3))]

#### glc6
	GM_glc6_o=[np.sum(m_glc6_o[:,:,:,x] * GM_mask)/np.count_nonzero(m_glc6_o[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_glc6_o,3))]
	WM_glc6_o=[np.sum(m_glc6_o[:,:,:,x] * WM_mask)/np.count_nonzero(m_glc6_o[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_glc6_o,3))]

	GM_glc6_n=[np.sum(m_glc6_n[:,:,:,x] * GM_mask)/np.count_nonzero(m_glc6_n[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_glc6_n,3))]
	WM_glc6_n=[np.sum(m_glc6_n[:,:,:,x] * WM_mask)/np.count_nonzero(m_glc6_n[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_glc6_n,3))]

#### tCr 
	tCr=[np.sum(m_tCr[:,:,:,x] * Mask)/np.count_nonzero(m_tCr[:,:,:,x] * Mask) for x in np.arange(0,np.size(m_tCr,3))]
	GM_tCr=[np.sum(m_tCr[:,:,:,x] * GM_mask)/np.count_nonzero(m_tCr[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_tCr,3))]
	WM_tCr=[np.sum(m_tCr[:,:,:,x] * WM_mask)/np.count_nonzero(m_tCr[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_tCr,3))]

	tNAA=[np.sum(m_NAA[:,:,:,x] * Mask)/np.count_nonzero(m_NAA[:,:,:,x] * Mask) for x in np.arange(0,np.size(m_NAA,3))]
	GM_tNAA=[np.sum(m_NAA[:,:,:,x] * GM_mask)/np.count_nonzero(m_NAA[:,:,:,x] * GM_mask) for x in np.arange(0,np.size(m_NAA,3))]
	WM_tNAA=[np.sum(m_NAA[:,:,:,x] * WM_mask)/np.count_nonzero(m_NAA[:,:,:,x] * WM_mask) for x in np.arange(0,np.size(m_NAA,3))]

######## ratio for tables

	GM_glx_o_tCr=[GM_glx_o[x] / GM_tCr[x] for x in np.arange(0,np.size(GM_glx_o))]
	WM_glx_o_tCr=[WM_glx_o[x] / WM_tCr[x] for x in np.arange(0,np.size(WM_glx_o))]
	
	GM_glx_n_tCr=[GM_glx_n[x] / GM_tCr[x] for x in np.arange(0,np.size(GM_glx_n))]
	WM_glx_n_tCr=[WM_glx_n[x] / WM_tCr[x] for x in np.arange(0,np.size(WM_glx_n))]

	GM_glc6_o_tCr=[GM_glc6_o[x] / GM_tCr[x] for x in np.arange(0,np.size(GM_glc6_o))]
	WM_glc6_o_tCr=[WM_glc6_o[x] / WM_tCr[x] for x in np.arange(0,np.size(WM_glc6_o))]
	
	GM_glc6_n_tCr=[GM_glc6_n[x] / GM_tCr[x] for x in np.arange(0,np.size(GM_glc6_n))]
	WM_glc6_n_tCr=[WM_glc6_n[x] / WM_tCr[x] for x in np.arange(0,np.size(WM_glc6_n))]

	GM_tNAA_tCr=[GM_tNAA[x] / GM_tCr[x] for x in np.arange(0,np.size(GM_tNAA))]
	WM_tNAA_tCr=[WM_tNAA[x] / WM_tCr[x] for x in np.arange(0,np.size(WM_tNAA))]
	
#### Linear glx_o
	
	GM_glx_o_linreg=stats.linregress(time_real,GM_glx_o);
	WM_glx_o_linreg=stats.linregress(time_real,WM_glx_o);

	plot_GM_glx_o_linreg=time_fit_plot*GM_glx_o_linreg[0]+GM_glx_o_linreg[1];
	plot_WM_glx_o_linreg=time_fit_plot*WM_glx_o_linreg[0]+WM_glx_o_linreg[1];
	

#### Linear glc6_o
	GM_glc6_o_linreg=stats.linregress(time_real,GM_glc6_o);
	WM_glc6_o_linreg=stats.linregress(time_real,WM_glc6_o);
	plot_GM_glc6_o_linreg=time_fit_plot*GM_glc6_o_linreg[0]+GM_glc6_o_linreg[1];
	plot_WM_glc6_o_linreg=time_fit_plot*WM_glc6_o_linreg[0]+WM_glc6_o_linreg[1];
	

###### GLX

	GM_glx_o_all_points.append(GM_glx_o)
	WM_glx_o_all_points.append(WM_glx_o)
	
	GM_glx_n_all_points.append(GM_glx_n)
	WM_glx_n_all_points.append(WM_glx_n)
	
	GM_glx_o_all_slopes.append(GM_glx_o_linreg[0]);
	GM_glx_o_all_intercept.append(GM_glx_o_linreg[1]);
	GM_glx_o_all_rvalue.append(GM_glx_o_linreg[2]);
	GM_glx_o_all_pvalue.append(GM_glx_o_linreg[3]);
	WM_glx_o_all_slopes.append(WM_glx_o_linreg[0]);
	WM_glx_o_all_intercept.append(WM_glx_o_linreg[1]);
	WM_glx_o_all_rvalue.append(WM_glx_o_linreg[2]);
	WM_glx_o_all_pvalue.append(WM_glx_o_linreg[3]);

	GM_glx_o_all_linfit.append(plot_GM_glx_o_linreg);
	WM_glx_o_all_linfit.append(plot_WM_glx_o_linreg);


###### Glc6

	GM_glc6_o_all_points.append(GM_glc6_o)
	WM_glc6_o_all_points.append(WM_glc6_o)
	
	GM_glc6_n_all_points.append(GM_glc6_n)
	WM_glc6_n_all_points.append(WM_glc6_n)
	

	GM_glc6_o_all_slopes.append(GM_glc6_o_linreg[0]);
	GM_glc6_o_all_intercept.append(GM_glc6_o_linreg[1]);
	GM_glc6_o_all_rvalue.append(GM_glc6_o_linreg[2]);
	GM_glc6_o_all_pvalue.append(GM_glc6_o_linreg[3]);
	WM_glc6_o_all_slopes.append(WM_glc6_o_linreg[0]);
	WM_glc6_o_all_intercept.append(WM_glc6_o_linreg[1]);
	WM_glc6_o_all_rvalue.append(WM_glc6_o_linreg[2]);
	WM_glc6_o_all_pvalue.append(WM_glc6_o_linreg[3]);
	
	GM_glc6_o_all_linfit.append(plot_GM_glc6_o_linreg);
	WM_glc6_o_all_linfit.append(plot_WM_glc6_o_linreg);

##### tCr
	tCr_all_points.append(tCr)
	GM_tCr_all_points.append(GM_tCr)
	WM_tCr_all_points.append(WM_tCr)
##### tNAA
	tNAA_all_points.append(tNAA)
	GM_tNAA_all_points.append(GM_tNAA)
	WM_tNAA_all_points.append(WM_tNAA)



#### tNAA / tCr
	GM_tNAA_tCr_all_points.append(GM_tNAA_tCr)
	WM_tNAA_tCr_all_points.append(WM_tNAA_tCr)

###### GLX/tCr

	GM_glx_tCr_o_all_points.append(GM_glx_o_tCr)
	WM_glx_tCr_o_all_points.append(WM_glx_o_tCr)

	GM_glx_tCr_n_all_points.append(GM_glx_n_tCr)
	WM_glx_tCr_n_all_points.append(WM_glx_n_tCr)
	

###### Glc6/tCr

	GM_glc6_tCr_o_all_points.append(GM_glc6_o_tCr)
	WM_glc6_tCr_o_all_points.append(WM_glc6_o_tCr)

	GM_glc6_tCr_n_all_points.append(GM_glc6_n_tCr)
	WM_glc6_tCr_n_all_points.append(WM_glc6_n_tCr)
	

	os.chdir('../../')











########################################################### END DATA READ

####### averaging Glx & Glx/tCr for pooled exp fitting 

avg_GM_glx_o=np.mean(GM_glx_o_all_points,0)
avg_WM_glx_o=np.mean(WM_glx_o_all_points,0)

avg_GM_glx_n=np.mean(GM_glx_n_all_points,0)
avg_WM_glx_n=np.mean(WM_glx_n_all_points,0)

std_GM_glx_o=np.std(GM_glx_o_all_points,0)
std_WM_glx_o=np.std(WM_glx_o_all_points,0)

std_GM_glx_n=np.std(GM_glx_n_all_points,0)
std_WM_glx_n=np.std(WM_glx_n_all_points,0)

avg_GM_glc6_o=np.mean(GM_glc6_o_all_points,0)
avg_WM_glc6_o=np.mean(WM_glc6_o_all_points,0)

avg_GM_glc6_n=np.mean(GM_glc6_n_all_points,0)
avg_WM_glc6_n=np.mean(WM_glc6_n_all_points,0)

std_GM_glc6_o=np.std(GM_glc6_o_all_points,0)
std_WM_glc6_o=np.std(WM_glc6_o_all_points,0)

std_GM_glc6_n=np.std(GM_glc6_n_all_points,0)
std_WM_glc6_n=np.std(WM_glc6_n_all_points,0)

avg_GM_tCr=np.mean(GM_tCr_all_points,0)
avg_WM_tCr=np.mean(WM_tCr_all_points,0)

std_GM_tCr=np.mean(GM_tCr_all_points,0)
std_WM_tCr=np.mean(WM_tCr_all_points,0)

avg_GM_tNAA=np.mean(GM_tNAA_all_points,0)
avg_WM_tNAA=np.mean(WM_tNAA_all_points,0)

std_GM_tNAA=np.mean(GM_tNAA_all_points,0)
std_WM_tNAA=np.mean(WM_tNAA_all_points,0)


avg_GM_glx_tCr=[avg_GM_glx_o[x]/avg_GM_tCr[x] for x in np.arange(0,14)]
avg_WM_glx_tCr=[avg_WM_glx_o[x]/avg_WM_tCr[x] for x in np.arange(0,14)]

std_GM_glx_tCr=np.std(GM_glx_tCr_o_all_points,0)
std_WM_glx_tCr=np.std(WM_glx_tCr_o_all_points,0)

avg_GM_glc6_tCr=[avg_GM_glc6_o[x]/avg_GM_tCr[x] for x in np.arange(0,14)]
avg_WM_glc6_tCr=[avg_WM_glc6_o[x]/avg_WM_tCr[x] for x in np.arange(0,14)]

std_GM_glc6_tCr=np.std(GM_glc6_tCr_o_all_points,0)
std_WM_glc6_tCr=np.std(WM_glc6_tCr_o_all_points,0)

######## averaging linear 

GM_glx_avg_linfit=stats.linregress(time_real,avg_GM_glx_o)
WM_glx_avg_linfit=stats.linregress(time_real,avg_WM_glx_o)


GM_glx_avg_linfitplot=time_fit_plot*GM_glx_avg_linfit[0]+GM_glx_avg_linfit[1]
WM_glx_avg_linfitplot=time_fit_plot*WM_glx_avg_linfit[0]+WM_glx_avg_linfit[1]


GM_glc6_avg_linfit=stats.linregress(time_real,avg_GM_glc6_o)
WM_glc6_avg_linfit=stats.linregress(time_real,avg_WM_glc6_o)

GM_glc6_avg_linfitplot=time_fit_plot*GM_glc6_avg_linfit[0]+GM_glc6_avg_linfit[1]
WM_glc6_avg_linfitplot=time_fit_plot*WM_glc6_avg_linfit[0]+WM_glc6_avg_linfit[1]



########### SAVE ALL PLOTS

	

plt.figure('glx_linfit_time_courses_new');
#[plt.plot(time_real,GM_glx_o_all_points[x],marker='o',ls='none',label='Subject ' + str(x+1) + ': GM: ' + 'slope = ' + str(np.round(GM_glx_o_all_slopes[x],2))+' / min') for x in np.arange(0,np.size(GM_glx_o_all_slopes))]
#[plt.plot(time_real,WM_glx_o_all_points[x],marker='o',ls='none',label='Subject ' + str(x+1) + ': WM: ' + 'slope = ' + str(np.round(WM_glx_o_all_slopes[x],2))+' /  min') for x in np.arange(0,np.size(GM_glx_o_all_slopes))]
[plt.plot(time_fit_plot,GM_glx_o_all_linfit[x],lw=0.05,color='blue',label='#' + str(x+1) + ': GM ' + str(np.round(GM_glx_o_all_slopes[x],0))+'r='+str(np.round(GM_glx_o_all_rvalue[x],2))+' P<'+str(np.round(GM_glx_o_all_pvalue[x],4))) for x in np.arange(0,np.shape(GM_glx_o_all_linfit)[0])];
[plt.plot(time_fit_plot,GM_glx_o_all_linfit[x],lw=0.05,color='blue',label='#' + str(x+1) + ': WM' + str(np.round(WM_glx_o_all_slopes[x],0))+'r='+str(np.round(WM_glx_o_all_rvalue[x],2))+' P<'+str(np.round(WM_glx_o_all_pvalue[x],4))) for x in np.arange(0,np.shape(WM_glx_o_all_linfit)[0])];
[plt.plot(time_fit_plot,GM_glx_o_all_linfit[x],lw=0.05,color='blue') for x in np.arange(0,np.shape(GM_glx_o_all_linfit)[0])];
[plt.plot(time_fit_plot,WM_glx_o_all_linfit[x],lw=0.05,color='red') for x in np.arange(0,np.shape(WM_glx_o_all_linfit)[0])];
plt.errorbar(time_real,avg_GM_glx_o,yerr=std_GM_glx_o,capsize=3,elinewidth=0.5,color='blue',marker='o',ls='none',label='GM')
plt.errorbar(time_real,avg_WM_glx_o,yerr=std_WM_glx_o,capsize=3,elinewidth=0.5,color='red',marker='o',ls='none',label='WM')
plt.plot(time_fit_plot,GM_glx_avg_linfitplot,color='blue')
plt.plot(time_fit_plot,WM_glx_avg_linfitplot,color='red')
plt.legend()
plt.savefig('plots/glx_lin_fit_new.pdf');
	


plt.figure('glx_linfit_avg');
plt.errorbar(time_real,avg_GM_glx_o/100,yerr=std_GM_glx_o/100,capsize=3,elinewidth=0.5,color='blue',marker='o',ls='none',label='GM: ' + 'slopes = ' + str(np.round(np.mean(GM_glx_o_all_slopes),0))+'+-'+str(np.round(np.std(GM_glx_o_all_slopes)))+' a.u./min')
plt.errorbar(time_real,avg_WM_glx_o/100,yerr=std_WM_glx_o/100,capsize=3,elinewidth=0.5,color='red',marker='o',ls='none',label='WM: ' + 'slopes = ' + str(np.round(np.mean(WM_glx_o_all_slopes),0))+'+-'+str(np.round(np.std(WM_glx_o_all_slopes)))+' a.u./min')
plt.plot(time_fit_plot,GM_glx_avg_linfitplot/100,color='blue')
plt.plot(time_fit_plot,WM_glx_avg_linfitplot/100,color='red') 
plt.legend()
plt.savefig('plots/glx_linfit_avg.pdf');




plt.figure('glc6_linfit_avg');
plt.errorbar(time_real,avg_GM_glc6_o/100,yerr=std_GM_glc6_o/100,capsize=3,elinewidth=0.5,color='blue',marker='o',ls='none',label='GM: ' + 'slopes = ' + str(np.round(np.mean(GM_glc6_o_all_slopes),0))+'+-'+str(np.round(np.std(GM_glc6_o_all_slopes)))+' a.u./min')
plt.errorbar(time_real,avg_WM_glc6_o/100,yerr=std_WM_glc6_o/100,capsize=3,elinewidth=0.5,color='red',marker='o',ls='none',label='WM: ' + 'slopes = ' + str(np.round(np.mean(WM_glc6_o_all_slopes),0))+'+-'+str(np.round(np.std(WM_glc6_o_all_slopes)))+' a.u./min')
plt.plot(time_fit_plot,GM_glc6_avg_linfitplot/100,color='blue')
plt.plot(time_fit_plot,WM_glc6_avg_linfitplot/100,color='red') 
plt.ylim(0.15e4,1.05e4)
plt.legend()
plt.savefig('plots/glc6_linfit_avg.pdf');

	
####### statistics


GM_glx_o_sorted_stat=[[GM_glx_o_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]
WM_glx_o_sorted_stat=[[WM_glx_o_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]

GM_glx_n_sorted_stat=[[GM_glx_n_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]
WM_glx_n_sorted_stat=[[WM_glx_n_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]

GM_glc6_o_sorted_stat=[[GM_glc6_o_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]
WM_glc6_o_sorted_stat=[[WM_glc6_o_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]

GM_glc6_n_sorted_stat=[[GM_glc6_n_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]
WM_glc6_n_sorted_stat=[[WM_glc6_n_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]

sorted_tCr=[[tCr_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]

GM_sorted_tCr=[[GM_tCr_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]
WM_sorted_tCr=[[WM_tCr_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]

sorted_tNAA=[[tNAA_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]

GM_sorted_tNAA=[[GM_tNAA_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]
WM_sorted_tNAA=[[WM_tNAA_all_points[x][y] for x in np.arange(0,np.size(subjects))] for y in np.arange(0,14)]




data_bar_mean=[[np.mean(GM_glx_o_sorted_stat[0]),np.mean(GM_glx_n_sorted_stat[0]),np.mean(GM_glc6_o_sorted_stat[0]),np.mean(GM_glc6_n_sorted_stat[0]),np.mean(GM_sorted_tCr[0]),np.mean(GM_sorted_tNAA[0])],[np.mean(GM_glx_o_sorted_stat[13]),np.mean(GM_glx_n_sorted_stat[13]),np.mean(GM_glc6_o_sorted_stat[13]),np.mean(GM_glc6_n_sorted_stat[13]),np.mean(GM_sorted_tCr[13]),np.mean(GM_sorted_tNAA[13])],[np.mean(WM_glx_o_sorted_stat[0]),np.mean(WM_glx_n_sorted_stat[0]),np.mean(WM_glc6_o_sorted_stat[0]),np.mean(WM_glc6_n_sorted_stat[0]),np.mean(WM_sorted_tCr[0]),np.mean(WM_sorted_tNAA[0])],[np.mean(WM_glx_o_sorted_stat[13]),np.mean(WM_glx_n_sorted_stat[13]),np.mean(WM_glc6_o_sorted_stat[13]),np.mean(WM_glc6_n_sorted_stat[13]),np.mean(WM_sorted_tCr[13]),np.mean(WM_sorted_tNAA[13])]]
data_bar_std=[[np.std(GM_glx_o_sorted_stat[0]),np.std(GM_glx_n_sorted_stat[0]),np.std(GM_glc6_o_sorted_stat[0]),np.std(GM_glc6_n_sorted_stat[0]),np.std(GM_sorted_tCr[0]),np.std(GM_sorted_tNAA[0])],[np.std(GM_glx_o_sorted_stat[13]),np.std(GM_glx_n_sorted_stat[13]),np.std(GM_glc6_o_sorted_stat[13]),np.std(GM_glc6_n_sorted_stat[13]),np.std(GM_sorted_tCr[13]),np.std(GM_sorted_tNAA[13])],[np.std(WM_glx_o_sorted_stat[0]),np.std(WM_glx_n_sorted_stat[0]),np.std(WM_glc6_o_sorted_stat[0]),np.std(WM_glc6_n_sorted_stat[0]),np.std(WM_sorted_tCr[0]),np.std(WM_sorted_tNAA[0])],[np.std(WM_glx_o_sorted_stat[13]),np.std(WM_glx_n_sorted_stat[13]),np.std(WM_glc6_o_sorted_stat[13]),np.std(WM_glc6_n_sorted_stat[13]),np.std(WM_sorted_tCr[13]),np.std(WM_sorted_tNAA[13])]]

plt.figure('barplot')

xpos=np.arange(6)

plt.bar(xpos+0.00,data_bar_mean[0],yerr=data_bar_std[0],color='blue',width=0.2,label='GM first point')
plt.bar(xpos+0.2,data_bar_mean[1],yerr=data_bar_std[1],color='none',edgecolor='blue',width=0.2,label='GM last point')
plt.bar(xpos+0.4,data_bar_mean[2],yerr=data_bar_std[2],color='red',width=0.2,label='WM first point')
plt.bar(xpos+0.6,data_bar_mean[3],yerr=data_bar_std[3],color='none',edgecolor='red',width=0.2,label='WM last point')
plt.xticks(xpos,['Glx$_4$','Glx$_{2+3}$','Glc$_6$','Glc$_{1-5}$','tCr','tNAA'])
plt.ylabel('metabolite concentrations [a.u.]')
plt.legend()

plt.savefig('plots/barplot_metabolites.pdf')


######### Blood glucose

blood_glc=[]
blood_glc_time=[]
subject_names=['P3','P5','P6','P7','P8']
DEU_03_glc=[80,101,91,88,84,90]
DEU_03_time=[0,17,36,51,63,100]

DEU_05_glc=[76,109,130,110,112,86]
DEU_05_time=[0,19,32,50,63,101]

DEU_06_glc=[97,92,108,120,141,132]
DEU_06_time=[0,15,34,47,60,98]

DEU_07_glc=[90,102,164,149,147,119]
DEU_07_time=[0,17,36,51,63,100]

DEU_08_glc=[85,115,135,108,109,105]
DEU_08_time=[0,18,32,50,62,99]

blood_glc.append(DEU_03_glc)
blood_glc.append(DEU_05_glc)
blood_glc.append(DEU_06_glc)
blood_glc.append(DEU_07_glc)
blood_glc.append(DEU_08_glc)

blood_glc_time.append(DEU_03_time)
blood_glc_time.append(DEU_05_time)
blood_glc_time.append(DEU_06_time)
blood_glc_time.append(DEU_07_time)
blood_glc_time.append(DEU_08_time)

blood_glc_mean=np.mean(blood_glc,0)
blood_glc_std=np.std(blood_glc,0)
blood_glc_time_mean=np.mean(blood_glc_time,0)
plt.figure('blood glucose')
plt.errorbar(blood_glc_time_mean,blood_glc_mean,yerr=blood_glc_std,label='average blood glucose over subjects')
plt.legend()
plt.xlabel('time / [min]')
plt.ylabel('blood glucose [mg/dl]')
plt.ylim(60,160)

plt.savefig('plots/blood_glucose.pdf')


##### Some calculations for Table 1


##concentration average over subjects
np.mean(GM_glx_tCr_o_all_points,0)[0]
np.std(GM_glx_tCr_o_all_points,0)[0]
np.mean(WM_glx_tCr_o_all_points,0)[0]
np.std(WM_glx_tCr_o_all_points,0)[0]
np.mean(GM_glx_tCr_n_all_points,0)[0]
np.std(GM_glx_tCr_n_all_points,0)[0]
np.mean(WM_glx_tCr_n_all_points,0)[0]
np.std(WM_glx_tCr_n_all_points,0)[0]


np.mean(GM_glc6_tCr_o_all_points,0)[0]
np.std(GM_glc6_tCr_o_all_points,0)[0]
np.mean(WM_glc6_tCr_o_all_points,0)[0]
np.std(WM_glc6_tCr_o_all_points,0)[0]
np.mean(GM_glc6_tCr_n_all_points,0)[0]
np.std(GM_glc6_tCr_n_all_points,0)[0]
np.mean(WM_glc6_tCr_n_all_points,0)[0]
np.std(WM_glc6_tCr_n_all_points,0)[0]

np.mean(GM_tNAA_tCr_all_points,0)[0]
np.std(GM_tNAA_tCr_all_points,0)[0]
np.mean(WM_tNAA_tCr_all_points,0)[0]
np.std(WM_tNAA_tCr_all_points,0)[0]



## CV over 60 min averaged over subjects
np.mean(np.std(GM_glx_tCr_o_all_points,1)/np.mean(GM_glx_tCr_o_all_points,1)*100)
np.std(np.std(GM_glx_tCr_o_all_points,1)/np.mean(GM_glx_tCr_o_all_points,1)*100)
np.mean(np.std(WM_glx_tCr_o_all_points,1)/np.mean(WM_glx_tCr_o_all_points,1)*100)
np.std(np.std(WM_glx_tCr_o_all_points,1)/np.mean(WM_glx_tCr_o_all_points,1)*100)
np.mean(np.std(GM_glx_tCr_n_all_points,1)/np.mean(GM_glx_tCr_n_all_points,1)*100)
np.std(np.std(GM_glx_tCr_n_all_points,1)/np.mean(GM_glx_tCr_n_all_points,1)*100)
np.mean(np.std(WM_glx_tCr_n_all_points,1)/np.mean(WM_glx_tCr_n_all_points,1)*100)
np.std(np.std(WM_glx_tCr_n_all_points,1)/np.mean(WM_glx_tCr_n_all_points,1)*100)

np.mean(np.std(GM_glc6_tCr_o_all_points,1)/np.mean(GM_glc6_tCr_o_all_points,1)*100)
np.std(np.std(GM_glc6_tCr_o_all_points,1)/np.mean(GM_glc6_tCr_o_all_points,1)*100)
np.mean(np.std(WM_glc6_tCr_o_all_points,1)/np.mean(WM_glc6_tCr_o_all_points,1)*100)
np.std(np.std(WM_glc6_tCr_o_all_points,1)/np.mean(WM_glc6_tCr_o_all_points,1)*100)
np.mean(np.std(GM_glc6_tCr_n_all_points,1)/np.mean(GM_glc6_tCr_n_all_points,1)*100)
np.std(np.std(GM_glc6_tCr_n_all_points,1)/np.mean(GM_glc6_tCr_n_all_points,1)*100)
np.mean(np.std(WM_glc6_tCr_n_all_points,1)/np.mean(WM_glc6_tCr_n_all_points,1)*100)
np.std(np.std(WM_glc6_tCr_n_all_points,1)/np.mean(WM_glc6_tCr_n_all_points,1)*100)


np.mean(np.std(GM_tNAA_tCr_all_points,1)/np.mean(GM_tNAA_tCr_all_points,1)*100)
np.std(np.std(GM_tNAA_tCr_all_points,1)/np.mean(GM_tNAA_tCr_all_points,1)*100)
np.mean(np.std(WM_tNAA_tCr_all_points,1)/np.mean(WM_tNAA_tCr_all_points,1)*100)
np.std(np.std(WM_tNAA_tCr_all_points,1)/np.mean(WM_tNAA_tCr_all_points,1)*100)


### Drop after 60 min (last/first time point)
np.mean([100-(GM_glx_tCr_o_all_points[x][13]/GM_glx_tCr_o_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(GM_glx_tCr_o_all_points[x][13]/GM_glx_tCr_o_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.mean([100-(WM_glx_tCr_o_all_points[x][13]/WM_glx_tCr_o_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(WM_glx_tCr_o_all_points[x][13]/WM_glx_tCr_o_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])

np.mean([100-(GM_glx_tCr_n_all_points[x][13]/GM_glx_tCr_n_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(GM_glx_tCr_n_all_points[x][13]/GM_glx_tCr_n_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.mean([100-(WM_glx_tCr_n_all_points[x][13]/WM_glx_tCr_n_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(WM_glx_tCr_n_all_points[x][13]/WM_glx_tCr_n_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])


np.mean([100-(GM_glc6_tCr_o_all_points[x][13]/GM_glc6_tCr_o_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(GM_glc6_tCr_o_all_points[x][13]/GM_glc6_tCr_o_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.mean([100-(WM_glc6_tCr_o_all_points[x][13]/WM_glc6_tCr_o_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(WM_glc6_tCr_o_all_points[x][13]/WM_glc6_tCr_o_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])

np.mean([100-(GM_glc6_tCr_n_all_points[x][13]/GM_glc6_tCr_n_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(GM_glc6_tCr_n_all_points[x][13]/GM_glc6_tCr_n_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.mean([100-(WM_glc6_tCr_n_all_points[x][13]/WM_glc6_tCr_n_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(WM_glc6_tCr_n_all_points[x][13]/WM_glc6_tCr_n_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])

np.mean([100-(GM_tNAA_tCr_all_points[x][13]/GM_tNAA_tCr_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(GM_tNAA_tCr_all_points[x][13]/GM_tNAA_tCr_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.mean([100-(WM_tNAA_tCr_all_points[x][13]/WM_tNAA_tCr_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])
np.std([100-(WM_tNAA_tCr_all_points[x][13]/WM_tNAA_tCr_all_points[x][0]*100) for x in np.arange(0,np.size(subjects))])


###### Correlation between Glx and Glc


GM_glx_corr=np.array(GM_glx_o_sorted_stat).flatten()
WM_glx_corr=np.array(WM_glx_o_sorted_stat).flatten()

GM_glc6_corr=np.array(GM_glc6_o_sorted_stat).flatten()
WM_glc6_corr=np.array(WM_glc6_o_sorted_stat).flatten()

GM_regress=stats.linregress(GM_glx_corr,GM_glc6_corr)
WM_regress=stats.linregress(WM_glx_corr,WM_glc6_corr)

GM_plt_regress=GM_regress[0]*GM_glx_corr+GM_regress[1]
WM_plt_regress=WM_regress[0]*WM_glx_corr+WM_regress[1]

plt.figure('Correlation Glx Glc')
plt.plot(GM_glx_corr,GM_glc6_corr,ls='none',marker='o',color='blue',label='GM')
plt.plot(WM_glx_corr,WM_glc6_corr,ls='none',marker='o',color='red',label='WM')
plt.plot(GM_glx_corr,GM_plt_regress,color='blue')
plt.plot(WM_glx_corr,WM_plt_regress,color='red')
plt.legend()
plt.xlabel('Glx4 concentration')
plt.ylabel('Glc6 concentration')

plt.savefig('plots/glx_glc_correlation.pdf')


