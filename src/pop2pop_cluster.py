from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.special import airy
import sys
plt.ion()



def generation(pop, mfit, mut):
    pop *= 1+(fitgrid - mfit - mut)
    pop[1:] += mut*pop[:-1]

    big_ii = pop>1000
    pop[big_ii] += np.sqrt(pop[big_ii])*np.random.normal(size=big_ii.sum())

    small_ii = (pop<=1000)&(pop>0)
    pop[small_ii] = np.random.poisson(pop[small_ii])
    pop[pop<0]=0
    return pop

def variance(pop):
    N = pop.sum()
    if N:
        avg = (pop*fitgrid).sum()/N
        return (pop*(fitgrid-avg)**2).sum()/N, avg
    else:
        return 0, 0

def corrfunc(a,b, tmax, both=True):
    n = len(a)
    ma = np.mean(a)
    mb = np.mean(b)
    acorr_fwd = np.array([np.mean(a[:n-ii]*b[ii:]) for ii in range(min(tmax,n/2))]) - ma*mb
    if both:
        acorr_bwd = np.array([np.mean(b[:n-ii]*a[ii:]) for ii in range(min(tmax,n/2))]) - ma*mb
        res = np.concatenate((acorr_bwd[::-1][:-1], acorr_fwd))
        return res
    else:
        return acorr_fwd

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("model divergence between cross-immune pathogen populations")
    parser.add_argument("-s", type=float, default=None, help="selection coefficient")
    parser.add_argument("-u", type=float, default=3e-3, help="mutation rate")
    parser.add_argument("--logNs", type=float, default=6, help="log Ns")
    parser.add_argument("--logsou", type=float, default=2, help="log s over mu")
    parser.add_argument("--N0", type=float, default=1e9, help="census population size")
    parser.add_argument("--nrep", type=int, default=3, help="repetitions with same burnin")
    parser.add_argument("--niter", type=int, default=3, help="iteration with separate burnin")
    parser.add_argument("--cross", type=float, help="time after which immunity is lost")
    parser.add_argument("--trun", type=int, default=10000, help="run time of simulation")
    parser.add_argument("--tburn", type=int, default=10000, help="burnin time of simulation")
    args = parser.parse_args()

    # s = 1.0/cross
    # mut = s*np.exp(-args.lnsu)
    # N0 = np.exp(lnNs)/s**2
    N0 = args.N0
    t_cross = args.cross
    if args.s:
        s=args.s
        mut=args.u
    else:
        s = np.exp(args.logNs)/N0*t_cross**2
        mut = np.exp(-args.logsou)*s

    if s<1e-3 or s>0.2 or u>0.05:
        print("out of range")
        sys.exit(0)

    n=int(1.0/s)*2
    fmax = 1.0
    fitgrid = np.linspace(-fmax,fmax,n+1)
    dx = fitgrid[1]-fitgrid[0]

    pop0=np.zeros((2,n+1))
    print('mu=%f, s=%f, N_0=%d'%(mut, dx, N0))

    t_burn = args.tburn
    t_ext = []
    for nii in range(args.niter):
        # set up initial condition for niter separate runs
        for i in [0,1]:
            pop0[i] = 0.00001*N0*np.exp(-fitgrid**2/0.05**2)
        mfit = np.array([fitgrid[n//2], fitgrid[n//2]])
        for t in xrange(t_burn):
            for i in [0,1]:
                pop0[i] = generation(pop0[i], mfit[i], mut)
            N = pop0.sum(axis=1)
            mfit += N/N0*0.5

            for i in [0,1]:
                if mfit[i]>dx:
                    ii = int(mfit[i]/dx)
                    pop0[i,:-ii] = pop0[i,ii:,]
                    pop0[i,-ii:]=0
                    mfit[i]-=dx*ii

        if any(N==0):
            print("initial population extinct")
            continue
        init_ratio = np.log(N[1]/N[0])
        init_N = N.mean()
        dmfit = mfit[1] - mfit[0]
        ssq_init = [variance(pop0[i])[0] for i in [0,1]]
        nose_pos_init = [np.max((pop0[i]>0)*fitgrid) for i in [0,1]]
        # repeat nrep times with the same initial condition
        for ri in range(args.nrep):
            mfit_list = [[],[]]
            pop = pop0.copy()
            mfit[:] = 0
            one_extinct = "--"
            both_extinct = "--"
            # the actual simulation
            for t in xrange(args.trun):
                for i in [0,1]:
                    pop[i] = generation(pop[i], mfit[i], mut)
                N = pop.sum(axis=1)

                mfit[0] += (N[0]+N[1]*np.exp(-1.0*t/t_cross))/N0
                mfit[1] += (N[1]+N[0]*np.exp(-1.0*t/t_cross))/N0
                for i in [0,1]:
                    ssq, avg = variance(pop[i])
                    nose_pos = np.max((pop[i]>0)*fitgrid)
                    mfit_list[i].append([mfit[i], avg, N[i], ssq, nose_pos-avg])

                if np.max(mfit)>dx:
                    ii = int(np.max(mfit)/dx)
                    for i in [0,1]:
                        pop[i,:-ii] = pop[i,ii:]
                        pop[i,-ii:]=0
                        mfit[i]-=dx*ii

                if any(N==0) and one_extinct=="--":
                    one_extinct = t if N[0] else -t
                if all(N==0) and both_extinct=="--":
                    both_extinct=t
                    break

            t_ext.append((init_N, init_ratio, dmfit, ssq_init[0], ssq_init[1],
                          nose_pos_init[0], nose_pos_init[1], str(one_extinct), str(both_extinct)))


    if args.s:
        fname = 'data/two_pop_N0_%1.1e_s_%f_u_%f_tc_%f.dat'%(N0,s,mut,t_cross)
    else:
        fname = 'data/two_pop_N0_%1.1e_logNs_%f_logsou_%f_tc_%f.dat'%(N0,args.logNs,args.logsou,t_cross)

    with open(fname, 'w') as ofile:
        ofile.write('#N, log(N2/N1), dmfit, var1, var2, nose1, nose2, first extinction, second extinction\n')
        for a in t_ext:
            ofile.write("%f, %f, %f, %f, %f, %f, %f, %s, %s\n"%a)
