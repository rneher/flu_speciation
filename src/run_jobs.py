import os
niter = 10
nrep = 20
for s in [0.005, 0.01, 0.02, 0.05]:
	for u in [0.001, 0.002, 0.005]:
		for N0 in [1e8, 1e9, 1e10]:
			for tc in [20, 50, 100, 200, 500]:
				cmd = ["src/pop2pop_cluster.py", "-u", u, "-s", s, "--cross", tc,
					   "--N0", N0, "--niter", niter, "--nrep", nrep]
				cmd_str = "sbatch submit_script.sh " + " ".join(map(str, cmd))
				print(cmd_str)
				os.system(cmd_str)
