from __future__ import division
import os
niter = 30
nrep = 20
for logsou in [1,2,3,4,5]:
	for tc in [20, 50, 100, 200, 500, 1000, 2000]:
		for N0 in [1e4, 1e5, 1e6]:
			for s in [0.3/tc, 1.0/tc, 2.0/tc]:
				cmd = ["src/pop2pop_cluster.py", "--logsou", logsou, "-s", s, "--cross", tc,
					   "--N0", N0, "--niter", niter, "--nrep", nrep]
				cmd_str = "sbatch submit_script.sh " + " ".join(map(str, cmd))
				print(cmd_str)
				os.system(cmd_str)
