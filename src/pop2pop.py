from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.special import airy
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
    n=50
    fmax = 1.0
    fitgrid = np.linspace(-fmax,fmax,n)
    dx = fitgrid[1]-fitgrid[0]

    pop0=np.zeros((2,n))
    N0 = 1000000000
    mut=0.003
    print('mu=%f, s=%f, N_0=%d'%(mut, dx, N0))

    t_burn = 10000
    res = []
    t_cross_times = np.array([10,20,50,100,200,300,500,1000,2000,5000,10000])
    for t_cross in t_cross_times:
        t_ext = []
        #print('t_cross=%f, t_sw=%f'%(t_cross, 4/dx))
        for nii in range(10):
            for i in [0,1]:
                pop0[i] = 0.00001*N0*np.exp(-fitgrid**2/0.05**2)
            mfit = np.array([fitgrid[n/2], fitgrid[n/2]])
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

            for ni in range(10):
                mfit_list = [[],[]]
                pop = pop0.copy()
                mfit[:] = 0
                extinct=False
                for t in xrange(10000):
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

                    if any(N==0):
                        t_ext.append(t)
                        extinct=True
                        break
                if not extinct:
                    t_ext.append(t)
                mfit_list = np.array(mfit_list)
        res.append(np.array(t_ext))
        r=res[-1]
        print("t_cross %f, split: %f, avg_ext: %f"%(t_cross, np.mean(r>9998),np.mean(r[r<9999])))

    for r in res:
        print("split: %f, avg_ext: %f"%(np.mean(r>9998),np.mean(r[r<9999])))

    plt.figure()
    split_prob = [np.mean(r>9998) for r in res]
    plt.plot(t_cross_times*dx, split_prob)
    plt.xlabel('$\sigma$')
    plt.ylabel('$P_{split}$')

    plt.figure()
    ext_time = [np.mean(r[r<9998]) for r in res]
    plt.plot(t_cross_times*dx, ext_time)
    plt.xlabel('$\sigma$')
    plt.ylabel('$T_{ext}$')
