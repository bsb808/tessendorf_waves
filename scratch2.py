from math import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np


class vertex_ocean:
    def __init__(self):
        self.x = 0.
        self.y = 0.
        self.z = 0.    # location
        self.nx = 0.
        self.ny = 0.
        self.nz = 0. # normal
        self.a = 0.
        self.b = 0.
        self.c = 0.    # htilde0
        self._a = 0.
        self._b = 0.
        self._c = 0.    # htilde0mk conjugate
        self.ox = 0.
        self.oy = 0.
        self.oz = 0.    # orig position

class complex_vector_normal:
    def __init__(self):
        self.h = complex(0.,0.)   # wave height
        self.D = array([0.,0.])   # displacement
        self.n = array([0.,0.,0.]) # normal

class cOcean:
    def __init__(self,N,A,w,length):
        self.N = N
        self.Nplus1 = self.N+1
        self.A = A
        self.w = w
        self.length = length
        self.geometry = False

        self.g = 9.81

        self.h_tilde = zeros((N,N),dtype=complex)
        self.h_tilde_slopex = zeros((N,N),dtype=complex)
        self.h_tilde_slopez = zeros((N,N),dtype=complex)
        self.h_tilde_slopex = zeros((N,N),dtype=complex)
        self.h_tilde_dx = zeros((N,N),dtype=complex)
        self.h_tilde_dz = zeros((N,N),dtype=complex)

        self.vertices = []
        for ii in range(self.Nplus1*self.Nplus1):
            self.vertices.append(vertex_ocean())
        '''
        for ii in range(self.Nplus1):
            v  = []
            for jj in range(self.Nplus1):
                v.append(vertex_ocean())
            self.vertices.append(v)
        '''
        self.indices = zeros(self.Nplus1*self.Nplus1*10)

        htilde0 = complex(0.,0.)
        htilde0mk_conj = complex(0.,0.)
        for m in range(self.Nplus1):
            for n in range(self.Nplus1):
                index = m*self.Nplus1+n
                htilde0 = self.hTilde_0(n,m)
                htilde0mk_conj = conj(self.hTilde_0(-n,-m))
                #print('%d of %d'%(index,len(self.vertices)))
                self.vertices[index].a = real(htilde0)
                self.vertices[index].b = imag(htilde0)
                self.vertices[index]._a = real(htilde0mk_conj)
                self.vertices[index]._b = imag(htilde0mk_conj)
                x = (n-self.N/2.0)*self.length/self.N
                self.vertices[index].ox = x
                self.vertices[index].x = x 
                y = 0.0
                self.vertices[index].oy = y
                self.vertices[index].y =  y
                z = (m-self.N/2.0)*self.length/self.N
                self.vertices[index].oz = z
                self.vertices[index].z =  z
                self.vertices[index].nx = 0.0;
                self.vertices[index].ny = 1.0;
                self.vertices[index].nz = 0.0;    

    def phillips(self,n,m):
        k = array([pi*(2*n-self.N)/self.length,
                   pi*(2*m-self.N)/self.length])
        k_length = np.linalg.norm(k)
        w_length = np.linalg.norm(self.w)
        if (k_length < 1e-6):
            return 0.0
        
        L = w_length**2 / self.g
        damping = 0.001
        l = (self.length*damping)**2
        k_dot_w = dot(k/np.linalg.norm(k),self.w/w_length)
        return (self.A*exp(-1.0/(k_length*L)**2)/(k_length**4) 
                * (k_dot_w)**2 *exp(-1.*(k_length)**2 * l**2))

    def hTilde_0(self,n,m):
        r = complex(normal(0,1),normal(0,1))
        p = self.phillips(n,m)
        return r * sqrt(p)/2.0

    def hTilde(self,t,n,m):
        index = m * self.Nplus1 + n
        htilde0 = complex(self.vertices[index].a,
                          self.vertices[index].b)
        htilde0mkconj = complex(self.vertices[index]._a,
                                self.vertices[index]._b)
        omegat = 0
        cos_ = cos(omegat)
        sin_ = sin(omegat)
        c0 = complex(cos_,sin_)
        c1 = complex(cos_,-1.0*sin_)
        res = htilde0*c0 + htilde0mkconj*c1
        return res

    def h_D_and_n(self,x,t):
        hh = complex(0.,0.)
        DD = array([0.,0.])
        nn = array([0.,0.,0.])

        for m in range(self.N):
            kz = 2.0*pi*(m-N/2.0)/self.length
            for n in range(self.N):
                kx = 2.0*pi*(n-N/2.0)/self.length
                k = array([kx,kz])
                
                k_length = np.linalg.norm(k)
                k_dot_x = dot(k,x)
                
                c = complex(cos(k_dot_x),sin(k_dot_x))
                htilde_c = self.hTilde(t,n,m)*c
                hh = hh + htilde_c
                nn = nn + array([-kx*imag(htilde_c),0.0,-kz*imag(htilde_c)])
            cvn = complex_vector_normal()
            cvn.h = hh
            cvn.D = DD
            cvn.n = nn
            return cvn
    
    def evaluateWaves(self,t):
        lamb= -1.0
        for m in range(self.N):
            for n in range(self.N):
                index = m*self.Nplus1 + n
                x = array([self.vertices[index].x,
                           self.vertices[index].z])
                hdn = self.h_D_and_n(x,t)
                self.vertices[index].y = float(real(hdn.h))

def getXZY(ocean):
    x = []
    z = []
    y = []
    for m in range(ocean.N):
        for n in range(ocean.N):
            index = m*ocean.Nplus1 + n
            x.append(ocean.vertices[index].x)
            z.append(ocean.vertices[index].z)
            y.append(ocean.vertices[index].y)

    return x,z,y

def getho(ocean):
    ho_real = []
    ho_imag = []
    for m in range(ocean.N):
        for n in range(ocean.N):
            index = m*ocean.Nplus1 + n
            ho_real.append(float(ocean.vertices[index].a))
            ho_imag.append(float(ocean.vertices[index].b))
    return ho_real, ho_imag

N = 32
Lx = 64.0
w = array([0.,5.])
A = 0.0005
ocean = cOcean(N,A,w,Lx)
x,z,y = getXZY(ocean)
#print((x,z,y))

ocean.evaluateWaves(0.0)
x,z,y = getXZY(ocean)


fig=figure(2)
clf()
ax = fig.gca(projection='3d')
ax.plot_surface(x,z,y, color='b')       


fig=figure(3)
clf()
ax = fig.gca(projection='3d')
ax.plot_wireframe(x,z,y)       

fig=figure(4)
clf()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x,z,y)       
ax.set_xlabel('X [m]')
ax.set_ylabel('Z [m]')

hr, hi = getho(ocean)

fig=figure(5)
clf()
ax = fig.gca(projection='3d')
ax.plot_trisurf(x,z,hr)   



show()

