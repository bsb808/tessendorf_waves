from math import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

def phil(k,V,A):
    g = 9.81
    L = V**2/g
    
    return A*exp(-1.0/((k*L)**2))/k**4
    

figure(1)
clf()

V = range(1,10)
for v in V:
    P = []
    K = linspace(0.01,1,100)
    for k in K:
        P.append(phil(k,v,1.0))

    plot(K,P,label='V=%f'%v)
legend()
show()

N = 64
M = 64
Lx = 64.0
Lz = 64.0

ho = zeros((N,M),dtype=complex)
wx = 1.0
wz = 0.0
for n,ii in zip(arange(-N/2,N/2),range(N)):
    for m,jj in zip(arange(-M/2,M/2),range(M)):
        kx = 2.0*pi*n/Lx
        kz = 2.0*pi*m/Lz
        
        k = sqrt(kx**2+kz**2)
        if abs(k) < 1.0e-6:
            k = 1.0e-6
        A = 1.0
        V = 10.0
        g = 9.81
        L = V**2/g
        P = A*exp(-1.0/((k*L)**2))/k**4 * abs(dot([kx,kz],[wx,wz]))**2
        
        zr =normal(0,1)
        zi =normal(0,1)
        ho[ii,jj] =( 1/sqrt(2.0) * complex(zr,zi) * sqrt(P))
        

# over spatial
h = zeros((N,M),dtype=complex)
xx = zeros((N,M),dtype=complex)
# Loop through space
for n,ii in zip(arange(-N/2,N/2),range(N)):
    for m,jj in zip(arange(-M/2,M/2),range(M)):
        xx[ii,jj] = complex(n*Lx/N,m*Lz/M)
        # Loop throuh wave numbers
        sx = 0.0
        for nn,iii in zip(arange(-N/2,N/2),range(N)):
            for mm,jjj in zip(arange(-M/2,M/2),range(M)):
                kx = 2.0*pi*nn/Lx
                kz = 2.0*pi*mm/Lz
                k = complex(kx,kz)
                ee = exp(dot(complex(0,1)*k,xx[ii,jj]))
                sx+=ho[iii,jjj]*ee
        h[ii,jj] = abs(sx)


        
fig=figure(2)
clf()
ax = fig.gca(projection='3d')
ax.plot_surface(real(xx),imag(xx),h, color='b')       
show()

        

                   
        
        


