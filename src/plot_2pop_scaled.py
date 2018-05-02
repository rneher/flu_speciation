from __future__ import division, print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, glob
plt.ion()


def read_file(fname):
	bname = os.path.basename(fname)
	N0, s, logsou, tcross = map(float, bname[:-4].split('_')[3::2])
	var = []
	avg_N = []
	corr = []
	t_ext = []
	count = 0
	with open(fname) as ifile:
		for line in ifile:
			if line[0]!='#':
				entries = map(lambda x:x.strip(), line.strip().split(','))
				N, dn, df, var1, var2, n1, n2 = map(float, entries[:7])
				count+=1
				var.extend([var1, var2])
				avg_N.append(np.log(N))
				if entries[7]!='--':
					tmpt = int(entries[7])
					t_ext.append(abs(tmpt))
					corr.append([dn,df, var1-var2, n1-n2, np.sign(tmpt)*1/(np.abs(tmpt)+1)])
	corr = np.array(corr)
	if len(t_ext):
		return (N0, s, logsou, tcross,np.mean(var),np.exp(np.mean(avg_N)),  len(t_ext)/count, np.mean(t_ext), count), corr
	else:
		return (N0, s, logsou, tcross, 0, np.exp(np.mean(avg_N)), 0, np.nan, count), corr

if __name__=="__main__":
	files = sorted(glob.glob('data/*logsou*dat'))
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
	logsou = np.unique(res[:,2])
	cols = ['g', 'r', 'c', 'b', 'm']
	ls = ['-', '--', ':', '-.', '-']
	marker = ['o', 's', 'v', '<', '>']
	plt.figure(figsize=(12,8))
	for ni,N0 in enumerate(N0_vals):
		for si, s in enumerate(s_vals):
			for ui, lsou in enumerate(logsou):
				ind = (res[:,0]==N0)&(res[:,1]==s)&(res[:,2]==lsou)
				ind_good = ind&(res[:,-1]>400)
				if ind_good.sum():
					avg_fit_std = np.mean(res[ind_good,4])**0.5
					avg_N = np.mean(res[ind_good,5])
					print(np.log(avg_N*s)/lsou)
					Tc = lsou**0.5/s*np.log(avg_N)
					plt.scatter(res[ind_good,3]/Tc**1.0,1-res[ind_good,6],
								 c=cols[ni], marker=marker[si]) #, ls=ls[ni])

	plt.plot([0,10],[1,np.exp(-10/4)])
	plt.yscale('log')
	plt.ylabel('Both populations survive', fontsize=fs)
	plt.xlabel('Cross immunity x fitness standard deviation $t_c/T_c$', fontsize=fs)
	plt.savefig("figures/p_survival.pdf")

	corr_coeffs = np.ma.masked_invalid([np.concatenate((v, np.corrcoef(a.T)[:,-1])) for v, a in zip(res,corrs) if len(a)>20])
	print("nose difference survival correlation:", np.ma.median(corr_coeffs[:,-2]))