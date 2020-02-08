import numpy as np
from cndreg2d import estimator
import matplotlib.pyplot as plt


def single_map(map_path, mask_path):
    obj = estimator('CND_REG2D', map_path, mask_path)
    tmp = obj.run()
    '''
    returning std_shredshold, v0, v1, v2
    '''
    return tmp[:,0], tmp[:,1], tmp[:,2], tmp[:,3]


def routine_single_example():
    mask_name = './data/mask_stripe_60deg_128.fits'
    map_prefix = './data/rnd_g1_001x_k10_'
    nside = 128
    N = 10  # ensemble size
    
    lv = 100  # threshold levels, follows CND_REG2D
    stot = np.zeros(N)
    v0tot = np.zeros((N,lv))
    v1tot = np.zeros((N,lv))
    v2tot = np.zeros((N,lv))
    for i in range(N):
        map_name = map_prefix+str(i)+'.fits'
        stot, v0tot[i], v1tot[i], v2tot[i] = single_map(map_name, mask_name, nside)
    
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(9,9), sharex=True)
    ax[0].errorbar(stot, np.mean(v0tot, axis=0), np.std(v0tot, axis=0), marker='.',linestyle='dotted',color='firebrick')
    ax[1].errorbar(stot, np.mean(v1tot, axis=0), np.std(v1tot, axis=0), marker='.',linestyle='dotted',color='firebrick')
    ax[2].errorbar(stot, np.mean(v2tot, axis=0), np.std(v2tot, axis=0), marker='.',linestyle='dotted',color='firebrick')
    
    ax[2].xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
    ax[2].tick_params(axis='both', which='major', labelsize='20')
    ax[2].set_xlabel(r'threshold', fontsize=20)
    ax[0].set_ylabel(r'MFs $v_0$', fontsize=20)
    ax[1].set_ylabel(r'MFs $v_1$', fontsize=20)
    ax[2].set_ylabel(r'MFs $v_2$', fontsize=20)
    
    plt.savefig('single_example.png')

    
def routine_multiple_example():
    mask_name = './data/mask_stripe_60deg_128.fits'
    map_prefix = './data/rnd_g1_001x_k10_'
    nside = 128
    N = 10  # ensemble size
    
    lv = 100  # threshold levels
    stot = np.zeros(N)
    v0tot = np.zeros((N,lv))
    v1tot = np.zeros((N,lv))
    v2tot = np.zeros((N,lv))
    for i in range(N):
        map_name = map_prefix+str(i)+'.fits'
        stot, v0tot[i], v1tot[i], v2tot[i] = single_map(map_name, mask_name, nside)

    map_prefix = './data/rnd_g10_001x_k10_'
    stot2 = np.zeros(N)
    v0tot2 = np.zeros((N,lv))
    v1tot2 = np.zeros((N,lv))
    v2tot2 = np.zeros((N,lv))
    for i in range(N):
        map_name = map_prefix+str(i)+'.fits'
        stot2, v0tot2[i], v1tot2[i], v2tot2[i] = single_map(map_name, mask_name, nside)

    import matplotlib.ticker as ticker
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(9,9), sharex=True)
    ax[0].errorbar(stot, np.mean(v0tot, axis=0), np.std(v0tot, axis=0), marker='.',linestyle='dotted',color='firebrick', label=r'$\rho=1$')
    ax[1].errorbar(stot, np.mean(v1tot, axis=0), np.std(v1tot, axis=0), marker='.',linestyle='dotted',color='firebrick')
    ax[2].errorbar(stot, np.mean(v2tot, axis=0), np.std(v2tot, axis=0), marker='.',linestyle='dotted',color='firebrick')
    ax[0].errorbar(stot2, np.mean(v0tot2, axis=0), np.std(v0tot2, axis=0), marker='.',linestyle='dotted',color='steelblue', label=r'$\rho=10$')
    ax[1].errorbar(stot2, np.mean(v1tot2, axis=0), np.std(v1tot2, axis=0), marker='.',linestyle='dotted',color='steelblue')
    ax[2].errorbar(stot2, np.mean(v2tot2, axis=0), np.std(v2tot2, axis=0), marker='.',linestyle='dotted',color='steelblue')
    
    ax[0].legend(loc=0,fontsize='20')
    ax[2].xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
    ax[2].tick_params(axis='both', which='major', labelsize='20')
    ax[2].set_xlabel(r'threshold', fontsize=20)
    ax[0].set_ylabel(r'MFs $v_0$', fontsize=20)
    ax[1].set_ylabel(r'MFs $v_1$', fontsize=20)
    ax[2].set_ylabel(r'MFs $v_2$', fontsize=20)
    
    plt.savefig('multiple_example.png')

if __name__=='__main__':
    routine_single_example()
