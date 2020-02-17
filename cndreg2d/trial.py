import numpy as np
from cndreg2d import estimator
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def single_map(map_path, mask_path, map_nside, mfs_threshold):
    obj = estimator('CND_REG2D', map_path, mask_path, map_nside, mfs_threshold)
    tmp = obj.run()
    '''
    returning std_shredshold, v0, v1, v2
    '''
    return tmp[:,0], tmp[:,1], tmp[:,2], tmp[:,3]


def test_example():
	map_name = 'map_cmb_ns64.fits'
	mask_name = 'mask_gal_fsky0.80_ns64.fits'
	nside = 64
	lv = 30
	stot, v0tot, v1tot, v2tot = single_map(map_name, mask_name, nside, lv)

	fit, ax = plt.subplots(nrows=3, ncols=1, figsize=(9,9), sharex=True)
	ax[0].plot(stot, v0tot)
	ax[1].plot(stot, v1tot)
	ax[2].plot(stot, v2tot)
	ax[2].xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
	ax[0].tick_params(axis='both', which='major', labelsize='15')
	ax[1].tick_params(axis='both', which='major', labelsize='15')
	ax[2].tick_params(axis='both', which='major', labelsize='15')
	ax[2].set_xlabel(r'threshold', fontsize=15)
	ax[0].set_ylabel(r'MFs $v_0$', fontsize=15)
	ax[1].set_ylabel(r'MFs $v_1$', fontsize=15)
	ax[2].set_ylabel(r'MFs $v_2$', fontsize=15)
	plt.savefig('test_example.png')


if __name__=='__main__':
	test_example()
