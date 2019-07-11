
fname = '../waves.dft/test.dat'
with open(fname,mode='rb') as bfile:
    data = bfile.read()

import numpy as np
dd = np.fromfile(fname,dtype=np.float32)
N = 64
for ii in range(0,len(dd),N**2+1):
    print '%d: %f'%(ii,dd[ii])
