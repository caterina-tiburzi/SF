#!/usr/bin/env python

import numpy as np
import os, argparse, math, sys
import matplotlib.pyplot as plt
from sklearn import preprocessing
import numpy.lib.recfunctions as rfn



parser = argparse.ArgumentParser(description="Compute SF from ASCII file containing MJD, Value, ValueError")
parser.add_argument('-f', '--infile', type=str,
                    help="File containing MJD, Value, ValueError")
parser.add_argument('-sh', '--skip_header', type=int, default=0,
                    help="Number of header lines to skip [default=0]")
parser.add_argument('-e', '--eb_cut', type=float, default=0.0001,
                    help="Maximum accepted uncertainty [default=1e-4]")
parser.add_argument('-s', '--sampling', type=str, default='lin',
                    help='Sampling type (linear=lin, logarithmic=log) [default=lin]')
parser.add_argument('--binsize', type=float, default=7.,
                    help="In the case of LINEAR BINNING: size of the sampling bin; in the case of LOGARITHMIC BINNING: size of the initial bin [default=7]")
parser.add_argument('--logarithmic_binning_factor', type=float, default=1.5,
                    help='t_(n+1) = 10^(lbf) x t_(n), only required for LOGARITHMIC BINNING [default=0.5]')
parser.add_argument('-n', '--montecarlo_iterations', type=int, default=1000,
                    help="Iteration in the Monte Carlo to compute error bars on the SF values [default=1000]")
parser.add_argument('-w', '--weighted', type=str, default="y",
                    help="Apply weights to the SF [default=y]")
parser.add_argument('-wn', '--white_noise', type=str, default="y",
                    help="Apply white noise subtraction [default=y]")
parser.add_argument('-p', '--plot', type=str, default="n",
                    help="Show plots [default=n]")


args = parser.parse_args()


data = np.genfromtxt(args.infile,skip_header=args.skip_header,names='MJD, value, evalue')
data = data[data['evalue']<args.eb_cut]


plt.figure(1)
plt.errorbar(data['MJD'],data['value'],yerr=data['evalue'],fmt='ko',zorder=0)
plt.xlabel("MJD")
plt.ylabel("Value")
plt.savefig("timeseries.eps")

if args.plot != "n":
    plt.show()
    
plt.close(1)



#SF computation 

tobs = data['MJD'][-1]-data['MJD'][0]
nrep = args.montecarlo_iterations


#start Monte Carlo

for h in range(0,nrep):
    
    print "Iteration number",h

    newVALUE = []
    for z in range(len(data)):
        newVALUE.append(np.random.normal(data['value'][z], data['evalue'][z], 1))#create a new data array through randomly sampling the normal distributions of the data points
    
    tlv, sfv, wev = [], [], []#prepare the arrays for the individual values of lag, SF and weighting factor
        
    for i in range(len(data)):#fill the arrays
        for j in range(i+1,len(data)):
            sfv.append(np.power(newVALUE[i]-newVALUE[j],2))
            tlv.append(np.absolute(data['MJD'][i]-data['MJD'][j]))
            wev.append(np.power(data['evalue'][i],2)+np.power(data['evalue'][j],2))#weighting factor comes from a modification of the error propagation -- suggested by Julian
        
    tlv = np.core.records.fromrecords(tlv, names='tlv',formats='f8')#transform the arrays in structured numpy arrays
    sfv = np.core.records.fromrecords(sfv, names='sfv',formats='f8')
    wev = np.core.records.fromrecords(wev, names='wev',formats='f8')
    
    #create a joint structured array
    pairs = tlv.copy()
    pairs = rfn.append_fields(pairs, ['sfv','wev'], [sfv['sfv'],wev['wev']], usemask = False)
    pairs.sort(order='tlv')

    sfa, sfaw, lags, npairs = [], [], [], []#initialize arrays for the average SF 

    
    if args.sampling == 'lin':

        if (args.binsize):
            binsize = args.binsize

        #print "tobs",tobs, "binsize",binsize, "number:",int(np.floor(tobs/binsize))
        center = binsize/2.

    	for k in range(0,int(np.floor(tobs/binsize))):
            #print "Iteration:",k,"\n", "center:",center,"\n", "center+binsize",center + binsize, "\n", "tobs",tobs
    	    if k==0:
    	        lags.append(center)
    	        binsf = pairs[pairs['tlv']<binsize]
    	        binsf = binsf[binsf['tlv']>0]
                wnarray = pairs[pairs['tlv']<center]
                wn = np.sum(wnarray['sfv']*1./wnarray['wev'])/np.sum(1./wnarray['wev'])
    	    else:
    	        center=center + binsize
    	        lags.append(center)
    	        binsf = pairs[pairs['tlv']<center+binsize/2.]
    	        binsf = binsf[binsf['tlv']>center-binsize/2.]    	    
                
            sfa.append(np.mean(binsf['sfv']))#SFA contains the averaged, non-weighted SF
            sfaw.append(np.sum(binsf['sfv']/binsf['wev'])/np.sum(1./binsf['wev']))#SF contains the averaged, weighted SF
            npairs.append(len(binsf))

        

    elif args.sampling == 'log':
        
        if (args.binsize):
            binend = args.binsize
        if (args.logarithmic_binning_factor):
            lbf = args.logarithmic_binning_factor

        binstart = 0

        while (binend < tobs):
            
            if binstart == 0:
                center = 10**(np.log10(np.sqrt(binend)))
                lags.append(center)
                binsf = pairs[pairs['tlv']<binend]
    	        binsf = binsf[binsf['tlv']>0]
                wnarray = pairs[pairs['tlv']<center]
                wn = np.sum(wnarray['sfv']*1./wnarray['wev'])/np.sum(1./wnarray['wev'])
            
            else:
                center = 10.**(np.log10(np.sqrt(binend*binstart)))
                lags.append(center)
                binsf = pairs[pairs['tlv']<binend]
    	        binsf = binsf[binsf['tlv']>binstart]

            #print "start:",binstart, "center:",center, "end:",binend, "tobs:",tobs
            binstart = binend
            binend = binend*10**lbf 

            sfa.append(np.mean(binsf['sfv']))#SFA contains the averaged, non-weighted SF
            sfaw.append(np.sum(binsf['sfv']/binsf['wev'])/np.sum(1./binsf['wev']))#SF contains the averaged, weighted SF
            npairs.append(len(binsf))
            
    sfawn = sfa - wn#SFAWN contains the averaged, non-weighted, wn-subtracted SF
    sfawwn = sfaw - wn#SFAWN contains the averaged, weighted, wn-subtracted SF
    

    npairs = np.core.records.fromrecords(npairs, names='npairs',formats='f8')#transform everything in a structured numpy array
    lags = np.core.records.fromrecords(lags, names='lags',formats='f8')
    sfa = np.core.records.fromrecords(sfa, names='sfa',formats='f8')
    sfaw = np.core.records.fromrecords(sfaw, names='sfaw',formats='f8')
    
    sfawn = np.core.records.fromrecords(sfawn, names='sfawn',formats='f8')
    sfawwn = np.core.records.fromrecords(sfawwn, names='sfawwn',formats='f8')
    
        
    if h == 0:
        sfmat = lags.copy()

    if args.sampling == "lin":
	    if args.weighted == "y":
	        if args.white_noise == "y":
	            sfmat = rfn.append_fields(sfmat, 'sf%s'%(h),sfawwn['sfawwn'], usemask = False)
	            lab = '_linsampl_wt_wnsub'
	        else:
	            sfmat = rfn.append_fields(sfmat, 'sf%s'%(h),sfaw['sfaw'], usemask = False)
	            lab = '_linsampl_wt'
	            
	    else:
	        if args.white_noise == "y":
	            sfmat = rfn.append_fields(sfmat, 'sf%s'%(h),sfawn['sfawn'], usemask = False)
	            lab = '_linsampl_wnsub'
	        else:
	            sfmat = rfn.append_fields(sfmat, 'sf%s'%(h),sfa['sfa'], usemask = False)
	            lab = '_linsampl'

    elif args.sampling == "log":
	    if args.weighted == "y":
	        if args.white_noise == "y":
	            sfmat = rfn.append_fields(sfmat, 'sf%s'%(h),sfawwn['sfawwn'], usemask = False)
	            lab = '_logsampl_wt_wnsub'
	        else:
	            sfmat = rfn.append_fields(sfmat, 'sf%s'%(h),sfaw['sfaw'], usemask = False)
	            lab = '_logsampl_wt'
	            
	    else:
	        if args.white_noise == "y":
	            sfmat = rfn.append_fields(sfmat, 'sf%s'%(h),sfawn['sfawn'], usemask = False)
	            lab = '_logsampl_wnsub'
	        else:
	            sfmat = rfn.append_fields(sfmat, 'sf%s'%(h),sfa['sfa'], usemask = False)
	            lab = '_logsampl'

fig = plt.figure(2)

for h in range(0,nrep):
    plt.loglog(sfmat['lags'],sfmat['sf%s'%h])

plt.xlabel("Lag [days]")
plt.ylabel(r"SF [u$^2$]")
plt.grid(which='both')

#insert Kolmogorov #WORK ON THIS
#tau_diff = 5./(60.*24.)#days
#lags_kolm = np.arange(tau_diff, 4000.)
#f = 345.*10**6 #Hz
#K = 2.41*10**(-4)*10**(-12) #Hz^-2 cm^-3 pc s^-1
#SF_kolm= (f*K/(2*np.pi))**2 * (lags_kolm/tau_diff)**(5./3.)
#plt.loglog(lags_kolm,SF_kolm,color='g',linewidth=2,label='Kolmogorov')
plt.legend(loc='best')

plt.savefig("SF_all%s.eps"%(lab))
if args.plot != "n":
    plt.show()

plt.close(fig)


format = str("%.16f'")
for h in range(0,nrep):
    if h !=nrep-1:
        format = format + str(",'%.16f'")
    else:
        format = format + str(",'%.16f")



exec"np.savetxt(\"SF_all%s.dat\", sfmat, fmt =[b'%s'])"%(lab,format)

sfmat_all = np.genfromtxt("SF_all%s.dat"%(lab))
lags = sfmat_all[:,0]
sfmat_all = sfmat_all[:,1:]
sfmat_all_m = np.mean(sfmat_all,axis=1)
sfmat_all_s = np.sort(sfmat_all,axis=1)

minussigma =sfmat_all_s[:,int(round(nrep*0.16))]
plussigma =sfmat_all_s[:,int(round(nrep*0.84))]

fig = plt.figure(3)
plt.xlabel('Lag [days]')
plt.ylabel(r'SF [u$^2$]')

plt.loglog(lags,sfmat_all_m,'k-',label='Mean')
plt.loglog(lags,minussigma,'r-',label=r'+/-1 $\sigma$')
plt.loglog(lags,plussigma,'r-')
plt.grid(which='both')

#insert Kolmogorov
#tau_diff = 15./(60.*24.)#days
#lags_kolm = np.arange(tau_diff, 4000.)
#f = 345.*10**6 #Hz
#K = 2.41*10**(-4)*10**(-12) #Hz^-2 cm^-3 pc s^-1
#SF_kolm= (f*K/(2*np.pi))**2 * (lags_kolm/tau_diff)**(5./3.)  
#plt.loglog(lags_kolm,SF_kolm,color='g',label='Kolmogorov')

plt.legend(loc='best')
plt.tight_layout()

plt.savefig("SF%s.eps"%(lab))

if args.plot != "n":
        plt.show()
        
plt.close(fig)

f = open("SF%s.dat"%(lab),'w')
#print sfmat_all_m
for i in range(len(sfmat_all_m)):
    f.write("%.16f %.16f %.16f %.16f %d\n"%(lags[i],minussigma[i],sfmat_all_m[i],plussigma[i],npairs['npairs'][i]))

f.close()


#"""
