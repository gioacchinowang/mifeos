"""
Author: Jiaxin Wang (SJTU)
Email: jiaxin.wang@sjtu.edu.cn
"""
import numpy as np
import healpy as hp
from mifeos.tools.icy_decorator import icy


@icy
class hmfs(object):
    
    def __init__(self, maps, mask=None, norm=True, nlv=None, width=None):
        """HEALPix-based partial-sky edge-trimed Minkowski Functional Estimator.
        
        Parameters:
        -----------
        
        maps : numpy.ndarray
            1D, single map;
            2D, with row # of maps (nsamp,npix).

        norm : bool
            switch for pixel value normalization.
            
        nlv : int, (list/tuple/array)
            int: number of MFs thresholds, by default 20, ranging from [-3,+3].
            (list/tuple/array): a list of MFs thresholds.
            
        width : float
            width for delta function,
            dimensionless, by default 1.e-4
            
        mask : (list/tuple/numpy.ndarray)
            1D array of mask map.
        """
        self.maps = maps
        self.norm = norm
        self.mask = mask
        self.nlv = nlv
        self.width = width

    @property
    def maps(self):
        return self._maps
    
    @property
    def mask(self):
        return self._mask

    @property
    def norm(self):
        return self._norm
    
    @property
    def nlv(self):
        return self._nlv
    
    @property
    def levels(self):
        return self._levels
    
    @property
    def theta(self):
        return self._theta

    @property
    def phi(self):
        return self._phi
    
    @property
    def width(self):
        return self._width

    @property
    def mean(self):
        return self._mean

    @property
    def std(self):
        return self._std

    @maps.setter
    def maps(self, maps):
        assert isinstance(maps, np.ndarray)
        # reshape maps and register info
        if (len(maps.shape) == 1):
            self._nsamp = 1
            self._npix = len(maps)
        elif (len(maps.shape) == 2):
            self._nsamp, self._npix = maps.shape
        self._neff = self._npix
        self._nside = int(np.sqrt(self._npix//12))
        self._maps = (maps.reshape(self._nsamp,self._npix).copy()).astype(np.float32)
        # prepare pixel positions
        self._theta = np.zeros(self._npix,dtype=np.float32)
        self._phi = np.zeros(self._npix,dtype=np.float32)
        for i in range(self._npix):
            self._theta[i],self._phi[i] = hp.pix2ang(self._nside,i)

    @norm.setter
    def norm(self, norm):
        assert isinstance(norm, bool)
        self._norm = norm

    @mask.setter
    def mask(self, mask):
        if mask is None:
            self._mask = [*range(self._npix)]
            self._neff = self._npix
        elif isinstance(mask, (list,tuple,np.ndarray)):
            assert (len(mask) == self._npix)
            self._mask = list(np.arange(self._npix)[mask>0])
            self._neff = len(self._mask)
        else:
            raise TypeError('unsupported type {}'.format(type(mask)))
        # normalize maps
        self._mean = np.mean(self._maps[:,self._mask])
        self._std = np.std(self._maps[:,self._mask])
        if self._norm:
            self._maps = (self._maps - self._mean)/self._std
            self._mean = 0.
            self._std = 1.

    @nlv.setter
    def nlv(self, nlv):
        if nlv is None:
            self._nlv = 20
            self._levels = np.linspace(-3.,3.,self._nlv)*self._std + self._mean
        elif isinstance(nlv, int):
            self._nlv = nlv
            self._levels = np.linspace(-3.,3.,self._nlv)*self._std + self._mean
        elif isinstance(nlv, (list,tuple,np.ndarray)):
            self._levels = np.array(nlv,dtype=np.float32)
            self._nlv = len(self._levels)
        else:
            raise TypeError('unsupported type {}'.format(type(nlv)))

    @width.setter
    def width(self, width):
        if width is None:
            self._width = 1.e-4*self._std
        else:
            if (width < 0.):
                raise ValueError('width = {}, cannot be negative'.format(width))
            elif (width > 6.*self._std/(self._nlv-1)):
                raise ValueError('width = {}, better smaller than {}'.format(width,6.*self._std/(self._nlv-1)))
            else:
                self._width = width
 
    #def trim_mask(self):
    #    """find 1st and 2nd mask boundary layers"""
    #    _f = np.arange(self._npix)[self._mask>0]
    #    _b = list()
    #    for i in _f:
    #        _n = hp.get_all_neighbours(self._nside,self._theta[i],self._phi[i])
    #        if not self._mask[_n[_n>0]].all():
    #            _b.append(i)
    #    self._mask[_b] = 0
    #    _f = np.arange(self._npix)[self._mask>0]
    #    _b.clear()
    #    for i in _f:
    #        _n = hp.get_all_neighbours(self._nside,self._theta[i],self._phi[i])
    #        if not self._mask[_n[_n>0]].all():
    #            _b.append(i)
    #    self._mask[_b] = 0

    def sigmoid(self, x, a):
        return 1./(1.+np.exp(-x/a))

    #def dirac(self, x, a):
    #    return np.nan_to_num( np.exp(-(x/a)**2)/(a*np.sqrt(np.pi)) )

    #def d11(self, input_dtheta):
    #    """2st order map derivative wrt theta then theta"""
    #    _alm = hp.map2alm(input_dtheta)
    #    return hp.alm2map_der1(_alm,self._nside)[1]

    #def d22(self, input_dtheta,input_dphi):
    #    """2nd order map derivative wrt. theta then phi"""
    #    _alm = hp.map2alm(input_dphi)
    #    return hp.alm2map_der1(_alm,self._nside)[2]+input_dtheta*np.cos(self._theta)/np.sin(self._theta)

    #def d12(self, input_dphi):
    #    """2nd order map derivative wrt phi then phi"""
    #    _alm = hp.map2alm(input_dphi*np.sin(self._theta))
    #    return hp.alm2map_der1(_alm,self._nside)[1]/np.sin(self._theta)-input_dphi*np.cos(self._theta)/np.sin(self._theta)
    
    def gmfs(self, input_map):
        """calculates mask-edge-trimed MFs for a single field"""
        _alm = hp.map2alm(input_map)
        _map, _Dthe, _Dphi = hp.alm2map_der1(_alm,self._nside)
        #_Dthethe = self.d11(_Dthe)
        #_Dthephi = self.d12(_Dphi)
        #_Dphiphi = self.d22(_Dthe,_Dphi)
        #v2_num = (2.*_Dthe*_Dphi*_Dthephi-(_Dthe**2)*_Dphiphi-(_Dphi**2)*_Dthethe)
        v2_denom = (_Dthe**2+_Dphi**2)
        v0 = np.zeros(self._nlv)
        #v1 = np.zeros_like(v0)
        #v2 = np.zeros_like(v0)
        for i in range(self._nlv):
            tmp_sigmoid = self.sigmoid(input_map-self._levels[i],self._width)
            v0[i] = np.sum( tmp_sigmoid[self._mask] )/self._neff
            #tmp_delta = self.dirac(input_map-self._levels[i],self._width)
            #v1[i] = 0.25*np.sum( (tmp_delta*np.sqrt(v2_denom))[self._mask>0] )/self._neff
            #v2[i] = (0.5/np.pi)*np.sum( (tmp_delta*v2_num/v2_denom)[self._mask>0] )/self._neff
        return v0 #(v0, v1, v2)

    def run(self):
        """
        Returns:
        --------
            ( MFs thresholds, MFs V0, MFs V1, MFs V2)
        """
        v0tot = np.zeros((self._nsamp,self._nlv),dtype=np.float32)
        #v1tot = np.zeros((self._nsamp,self._nlv),dtype=np.float32)
        #v2tot = np.zeros((self._nsamp,self._nlv),dtype=np.float32)
        for i in range(self._nsamp):
            v0tot[i] = self.gmfs(self._maps[i])
            #v0tot[i], v1tot[i], v2tot[i] = self.gmfs(self._maps[i])
        return (self._levels, v0tot) #, v1tot, v2tot)
