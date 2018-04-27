#!/usr/bin/env python

import numpy as np
import os, argparse, math, sys
import matplotlib.pyplot as plt
from sklearn import preprocessing
import numpy.lib.recfunctions as rfn

parser = argparse.ArgumentParser(description="Compute SF")
parser.add_argument('-f', '--dmfile', type=str, nargs=1,
                    help="File containing DM variations, as MJD, DM, eDM")
parser.add_argument('--eb_cut', type=float, nargs='?', default=0.0001,
                    help="error bar threshold in the DM plot")
parser.add_argument('--a_cut', type=float, nargs='?',
                    help="Window on the Solar angle")
parser.add_argument('-binsize', '--binsize', type=float, nargs=1,
                    help="Length of the bin in the SF")

args = parser.parse_args()

dmfile = args.dmfile[0]
if (args.eb_cut):
    cut = args.eb_cut
if (args.a_cut):
    acut = args.a_cut
if (args.binsize):
    binsize = args.binsize[0]

data = np.genfromtxt(dmfile,names='name, date, MJD, DM, edm, angle')
data = data[data['angle']>acut]
data = data[data['edm']<cut]

tobs = data['MJD'][-1]-data['MJD'][0]

#meandt = np.ediff1d(data['MJD']).mean(axis=0)

sfv, tcv, va = [], [], []

for i in range(len(data)):
    for j in range(i+1,len(data)):
        sfv.append(np.power(data['DM'][i]-data['DM'][j],2))
        tcv.append(np.absolute(data['MJD'][i]-data['MJD'][j]))
        va.append(np.power(data['edm'][i],2)+np.power(data['edm'][j],2))

tcv = np.core.records.fromrecords(tcv, names='tcv',formats='f8')
sfv = np.core.records.fromrecords(sfv, names='sfv',formats='f8')
va = np.core.records.fromrecords(va, names='va',formats='f8')


pairs = tcv
pairs = rfn.append_fields(pairs, ['sfv','va'], [sfv['sfv'],va['va']], usemask = False)
pairs.sort(order='tcv')

wn = pairs[pairs['tcv']<binsize/2.]
wn = np.sum(wn['sfv']*1./wn['va'])/np.sum(1./wn['va'])

sca, scaw, lags = [], [], []
edge=binsize/2

for k in range(0,int(np.floor(tobs/binsize))):
    if k==0:
        center = 0
        lags.append(center)
        binsf = pairs[pairs['tcv']<edge]
        binsf = binsf[binsf['tcv']>0]
    else:
        center=k*binsize
        lags.append(center)
        binsf = pairs[pairs['tcv']<center+edge]
        binsf = binsf[binsf['tcv']>center-edge]                             
        sca.append(np.mean(binsf['sfv']))
    scaw.append(np.sum(binsf['sfv']/binsf['va'])/np.sum(1./binsf['va']))

scawn = sca - wn
scawwn = scaw - wn

lags = np.core.records.fromrecords(lags, names='lags',formats='f8')
sca = np.core.records.fromrecords(sca, names='sca',formats='f8')
scaw = np.core.records.fromrecords(scaw, names='scaw',formats='f8')

scawn = np.core.records.fromrecords(scawn, names='scawn',formats='f8')
scawwn = np.core.records.fromrecords(scawwn, names='scawwn',formats='f8')



sf = lags
sf = rfn.append_fields(sf, ['scaw','scawwn','sca','scawn'], [scaw['scaw'],scawwn['scawwn'],sca['sca'],scawn['scawn']], usemask = False)

sfcorr = sf[np.invert(np.isnan(sf['scaw']))]
sfcorr = sfcorr[np.invert(np.isnan(sfcorr['sca']))]

sfcorr = sf[np.invert(np.isnan(sf['scawwn']))]
sfcorr = sfcorr[np.invert(np.isnan(sfcorr['scawn']))]

plt.figure()
plt.loglog(sfcorr['lags'],sfcorr['scawwn'],color='k',linestyle='-')
plt.loglog(sfcorr['lags'],sfcorr['scaw'],color='r',linestyle='-')
plt.grid(which='both')
plt.show()

np.savetxt("SF_acut%s.dat"%(acut),sf, fmt=[b'%.16f','%.16f', '%.16f','%.16f','%.16f'])
