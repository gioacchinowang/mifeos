"""
ez CND_REG2D python wrapper
"""

import os
import subprocess
import healpy as hp
import numpy as np

class estimator(object):

    def __init__(self,
                 exe_path=None,
                 map_path=None,
                 msk_path=None,
                 nside=None):
        # current working directory
        self.wk_dir = os.getcwd()
        self.exe_path = exe_path
        self.map_path = map_path
        self.msk_path = msk_path
        self.nside = nside
        
    @property
    def wk_dir(self):
        return self._wk_dir

    @wk_dir.setter
    def wk_dir(self, wk_dir):
        self._wk_dir = wk_dir
        
    @property
    def exe_path(self):
        return self._exe_path
        
    @property
    def map_path(self):
        return self._map_path
        
    @property
    def msk_path(self):
        return sef._msk_path
        
    @property
    def nside(self):
        return self._nside
        
    @exe_path.setter
    def exe_path(self, exe_path):
        """
        by default hammurabiX executable "hamx" should be available in 'PATH'
        """
        assert isinstance(exe_path, str)
        self._exe_path = os.path.abspath(exe_path)
        self._executable = self._exe_path
    
    @map_path.setter
    def map_path(self, map_path):
        assert isinstance(map_path, str)
        self._map_path = os.path.abspath(map_path)
        assert os.path.isfile(self._map_path)
        
    @msk_path.setter
    def msk_path(self, msk_path):
        assert isinstance(msk_path, str)
        self._msk_path = os.path.abspath(msk_path)
        assert os.path.isfile(self._msk_path)
        
    @nside.setter
    def nside(self, nside):
        assert isinstance(nside, int)
        self._nside = nside
        
    def run(self):
        dumpfile = os.path.join(self._wk_dir, 'mf.dat')
        if os.path.isfile(dumpfile):
            os.remove(dumpfile)
        
        temp_process = subprocess.Popen([self._exe_path,
                                         self._map_path,
                                         self._msk_path,
                                         str(self._nside)],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT)
        temp_process.wait()
        if temp_process.returncode != 0:
            last_call_log, last_call_err = temp_process.communicate()
            print(last_call_log)
            print(last_call_err)
        
        result = np.loadtxt(dumpfile)
        os.remove(dumpfile)
        return result
