from math import *

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

N = 4
M = 4
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
        P = A*exp(-1.0/((k*L)**2))/k**4 * abs(dot([kx,kz],[wx,wz]))
        
        zr =normal(0,1)
        zi =normal(0,1)
        ho[ii,jj] =( 1/sqrt(2.0) * complex(zr,zi) * sqrt(P))


        

                   
        
        


