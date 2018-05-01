import os
niter = 30
nrep = 20
for logsou in [1,2,3,4,5]:
	for logNs in [3,4,5,6,7]:
		for N0 in [1e8, 1e10]:
			for tc in [20, 50, 100, 200, 500, 1000, 2000]:
				cmd = ["src/pop2pop_cluster.py", "--logsou", logsou, "--logNs", logNs, "--cross", tc,
					   "--N0", N0, "--niter", niter, "--nrep", nrep]
				cmd_str = "sbatch submit_script.sh " + " ".join(map(str, cmd))
				print(cmd_str)
				os.system(cmd_str)
