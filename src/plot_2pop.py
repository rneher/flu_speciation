from __future__ import division, print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, glob
plt.ion()


def read_file(fname):
	bname = os.path.basename(fname)
	N0, s, u, tcross = map(float, bname[:-4].split('_')[3::2])
	var = []
	corr = []
	t_ext = []
	count = 0
	with open(fname) as ifile:
		for line in ifile:
			if line[0]!='#':
				entries = map(lambda x:x.strip(), line.strip().split(','))
				dn, df, var1, var2, n1, n2 = map(float, entries[:6])
				count+=1
				var.extend([var1, var2])
				if entries[6]!='--':
					tmpt = int(entries[6])
					t_ext.append(abs(tmpt))
					corr.append([dn,df, var1-var2, n1-n2, np.sign(tmpt)*1/(np.abs(tmpt)+1)])
	corr = np.array(corr)
	if len(t_ext):
		return (N0, s, u, tcross,np.mean(var),np.mean(np.abs(corr[:,3])),  len(t_ext)/count, np.mean(t_ext), count), corr
	else:
		return (N0, s, u, tcross, 0, 0, 0, np.nan, count), corr

if __name__=="__main__":
	files = sorted(glob.glob('data/*dat'))
	res = []
	corrs = []
	fs=16
	for fname in files:
		a,corr = read_file(fname)
		res.append(a)
		corrs.append(corr)
	res = np.array(res)

	N0_vals = np.unique(res[:,0])
	s_vals = np.unique(res[:,1])
	mu_vals = np.unique(res[:,2])
	cols = ['g', 'r', 'c', 'b', 'm']
	ls = ['-', '--', ':', '-.', '-']
	marker = ['o', 's', 'v', '<', '>']
	plt.figure(figsize=(12,8))
	for ni,N0 in enumerate(N0_vals):
		for si, s in enumerate(s_vals):
			for ui, u in enumerate(mu_vals):
				ind = (res[:,0]==N0)&(res[:,1]==s)&(res[:,2]==u)
				ind_good = ind&(res[:,-1]>400)
				avg_fit_std = np.mean(res[ind_good,4])**0.5
				avg_nose_diff = np.mean(res[ind_good,5])

				# plt.scatter(res[ind_good,3]*avg_fit_std**1*(u*s**2)**-0.33/np.log(N0*s**0)**4,1-res[ind_good,5],
				# 			 c=cols[si], marker=marker[ui]) #, ls=ls[ni])
#				Tc = np.log(N0*s**2)**0.5/avg_fit_std
				# plt.scatter(res[ind_good,3]/Tc**1.0,1-res[ind_good,6],
				# 			 c=cols[ni], marker=marker[si]) #, ls=ls[ni])
				Tc = np.log(N0*s**2*u)/(s**2*u)**0.333
				plt.scatter(res[ind_good,3]/Tc**1.0,1-res[ind_good,6],
							 c=cols[si], marker=marker[ni]) #, ls=ls[ni])

				#plt.scatter(res[ind_good,3]*(u*s**2)**0.33/np.log(N0*s),1-res[ind_good,5],
				#			c=cols[si], marker=marker[ui])
	plt.yscale('log')
	plt.ylabel('Both populations survive', fontsize=fs)
	plt.xlabel('Cross immunity x fitness standard deviation $t_c \sigma/\log Ns^2$', fontsize=fs)
	plt.savefig("figures/p_survival.pdf")

	corr_coeffs = np.ma.masked_invalid([np.concatenate((v, np.corrcoef(a.T)[:,-1])) for v, a in zip(res,corrs) if len(a)>20])
	print("nose difference survival correlation:", np.ma.median(corr_coeffs[:,-2]))