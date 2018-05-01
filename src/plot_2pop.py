from __future__ import division, print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, glob

def read_file(fname):
	bname = os.path.basename(fname)
	N0, s, u, tcross = map(float, bname[:-4].split('_')[3::2])
	t_ext = []
	var = []
	count = 0
	with open(fname) as ifile:
		for line in ifile:
			if line[0]!='#':
				entries = map(lambda x:x.strip(), line.strip().split(','))
				dn, df, var1, var2, n1, n2 = map(float, entries[:6])
				count+=1
				var.extend([var1, var2])
				if entries[6]!='--':
					t_ext.append(abs(int(entries[6])))

	if len(t_ext):
		return N0, s, u, tcross,np.mean(var), len(t_ext)/count, np.mean(t_ext), count
	else:
		return N0, s, u, tcross, 0, 0, np.nan, count

if __name__=="__main__":
	files = sorted(glob.glob('data/*dat'))
	res = []
	for fname in files:
		res.append(read_file(fname))
	res = np.array(res)

	N0_vals = np.unique(res[:,0])
	s_vals = np.unique(res[:,1])
	mu_vals = np.unique(res[:,2])
	cols = ['g', 'r', 'c', 'b', 'm']
	ls = ['-', '--', ':', '-.']
	marker = ['o', 's', 'v', '<', '>']
	plt.figure()
	for ni,N0 in enumerate(N0_vals):
		for si, s in enumerate(s_vals):
			for ui, u in enumerate(mu_vals):
				ind = (res[:,0]==N0)&(res[:,1]==s)&(res[:,2]==u)
				print(N0, s, u, res[ind,-1].mean())
				ind_good = ind&(res[:,-1]>100)
				avg_fit_std = np.mean(res[ind_good,4])**0.5
				#plt.plot(res[ind,3]*avg_fit_std,1-res[ind,5], c=cols[si], ls=ls[ui])
				plt.scatter(res[ind_good,3]*avg_fit_std/np.log(N0),1-res[ind_good,5], c=cols[ni], marker=marker[ui])


