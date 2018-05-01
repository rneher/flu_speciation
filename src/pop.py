import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.special import airy
plt.ion()



def generation(pop, mfit, mut):
    pop *= 1+(fitgrid - mfit)
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
    parser.add_argument("--lnNs", "log of Ns")
    parser.add_argument("--lnsu", "log of s/u")
    parser.add_argument("--cross", "time after which immunity is lost")
    args = parser.parse_args()

    s = 1.0/cross
    mut = s*np.exp(-args.lnsu)
    N0 = np.exp(lnNs)/s**2

    n=int(1.0/s)
    fmax = 1.0
    fitgrid = np.linspace(-fmax,fmax,n)
    dx = fitgrid[1]-fitgrid[0]

    pop=np.zeros(n)
    # N0 = 1000000000
    # mut=0.003
    print('mu=%f, s=%f, N_0=%d'%(mut, dx, N0))

    pop[:] = 0.00001*N0*np.exp(-fitgrid**2/0.05**2)
    mfit = fitgrid[n/2]
    mfit_list = []
    for t in xrange(100000):
        pop = generation(pop, mfit, mut)
        N = pop.sum()
        mfit += N/N0

        if t>10000:
            ssq, avg = variance(pop)
            nose_pos = np.max((pop>0)*fitgrid)
            mfit_list.append([mfit, avg,N, ssq, nose_pos-avg])
            if t%10000==0:
                print(t)
        if mfit>dx:
            ii = int(mfit/dx)
            #print("backshift",t,mfit, ii,variance(pop), N)
            pop[:-ii] = pop[ii:]
            pop[-ii:]=0
            mfit-=dx*ii

    mfit_list=np.array(mfit_list)
    print("fitness stddev: %1.3e, %1.3e, %1.3e"%(np.sqrt(np.min(mfit_list[:,3])), np.sqrt(np.mean(mfit_list[:,3])), np.sqrt(np.max(mfit_list[:,3]))))
    print("Infected fraction: %1.3e, %1.3e, %1.3e"%(np.min(mfit_list[:,2])/N0, np.mean(mfit_list[:,2])/N0, np.max(mfit_list[:,2])/N0))
    print("Nose: %1.3e, %1.3e, %1.3e"%(np.min(mfit_list[:,-1]), np.mean(mfit_list[:,-1]), np.max(mfit_list[:,-1])))

    n_grid=20
    vec_field = np.zeros((n_grid,n_grid,2), dtype=float)
    vec_field_count = np.zeros((n_grid,n_grid), dtype=int)
    N_bins = np.linspace(4,9,n_grid)
    sq_bins = np.linspace(np.min(mfit_list[:,3]),np.max(mfit_list[:,3]),n_grid)
    dN=np.diff(np.log10(mfit_list[:,2]), axis=0)
    dsq=np.diff(mfit_list[:,3], axis=0)

    sq_ii = sq_bins.searchsorted(mfit_list[:,3])
    n_ii = N_bins.searchsorted(np.log10(mfit_list[:,2]+1))
    for ii,val in enumerate(dN):
        vec_field[sq_ii[ii],n_ii[ii]]+=[val, dsq[ii]]
        vec_field_count[sq_ii[ii],n_ii[ii]]+=1
    vec_field[:,:,0]/=(vec_field_count+.01)
    vec_field[:,:,1]/=(vec_field_count+.01)
    vec_field[:,:,0]/=vec_field[:,:,0].max()
    vec_field[:,:,1]/=vec_field[:,:,1].max()
    plt.figure()
    plt.quiver(N_bins, sq_bins, vec_field[:,:,0], vec_field[:,:,1])

    # plt.figure()
    # plt.plot(mfit_list[:,1] - mfit_list[:,0])
    # plt.plot(np.sqrt(mfit_list[:,3]))
    # plt.plot(np.sqrt(mfit_list[:,4]))

    # plt.figure()
    # plt.plot(np.log(mfit_list[:,2]))



    plt.figure()
    tmax=int(5/dx)
    plt.plot([0,0],[-1,1], lw=2, c='k', alpha=0.5)
    plt.plot([-tmax,tmax],[0,0], lw=2, c='k', alpha=0.5)
    dt = np.arange(-tmax+1, tmax)
    acorr = corrfunc(mfit_list[:,2], mfit_list[:,2],tmax)
    plt.plot(dt, acorr/np.abs(acorr).max(), label="N*N(t+tau)")

    acorr = corrfunc(mfit_list[:,3], mfit_list[:,2],tmax)
    plt.plot(dt, acorr/np.abs(acorr).max(), label="sigma^2(t)*N(t+tau)")

    acorr = corrfunc(mfit_list[:,1]-mfit_list[:,0], mfit_list[:,2],tmax)
    plt.plot(dt, acorr/np.abs(acorr).max(), label="dfit(t)*N(t+tau)")

    acorr = corrfunc(mfit_list[:,4], mfit_list[:,2],tmax)
    plt.plot(dt, acorr/np.abs(acorr).max(), label="nose(t)*N(t+tau)")
    plt.legend()

    plt.figure()
    bins = np.logspace(1,np.log10(N0), 100)
    y,x = np.histogram(mfit_list[:,2], bins=bins)
    bc = 0.5*(bins[:-1]+bins[1:])
    bw = np.diff(bins)
    plt.plot(bc,y)
    plt.xscale('log')
