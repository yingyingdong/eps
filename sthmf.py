import numpy as np
import ctypes
#from ctypes import *
import os
PS=ctypes.cdll.LoadLibrary('pshmf.so')

class DATA(ctypes.Structure):
     _fields_ = [("n", ctypes.c_int),
                ("z1", ctypes.POINTER(ctypes.c_double)),
                ("DZ", ctypes.POINTER(ctypes.c_double)),
                ("M", ctypes.POINTER(ctypes.c_double)),
                ("sigm", ctypes.POINTER(ctypes.c_double))]

z1   = 0.45
M    = np.arange(11,15,0.01)
sigm = np.zeros(M.shape[0],np.float)
redshift1 = np.ones(M.shape[0])*z1
data = DATA(M.shape[0],
            np.ctypeslib.as_ctypes(redshift1),
            np.ctypeslib.as_ctypes(np.zeros(M.shape[0])),
            np.ctypeslib.as_ctypes(M),
            np.ctypeslib.as_ctypes(sigm))



PS.prepare()
PS.main_power.restype = None
PS.main_power(ctypes.byref(data))
print(data.sigm[0:20])
