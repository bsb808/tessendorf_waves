

def pm(omega,omega_p):
    ''' 
    Return spectrum value as a function of freq and peak freq
    '''
    alpha = 0.0081
    beta = 0.74
    g = 9.81
    return alpha*g**2/(omega**5)*exp(-(5.0/4.0)*(omega_p/omega)**4)




figure(10)
clf()
Tps = [1]#2,4,6]
for Tp in Tps:
    omega_p = 2*pi/Tp
    omegas = linspace(omega_p*.25,omega_p*4,100)
    s = []
    for o in omegas:
        s.append(pm(o,omega_p))
    plot(omegas,s)

    scale = 1.25
    omegas = [omega_p/scale, omega_p, omega_p*scale]
    delos = [omega_p*(1-1/scale), omega_p*(scale-1/scale)/2.0, omega_p*(scale-1)]
    s = []
    a = []
    for o,d in zip(omegas,delos):
        s.append(pm(o,omega_p))
        a.append(sqrt(2*pm(o,omega_p)*d))
    print 2*pi/omega_p
    print 2*pi/array(omegas)
    print a     
    plot(omegas,s,'o')

grid(True)
xlabel('Omega [rad/s]')
ylabel('S(omega)')

             

print "----------------------"
figure(1)
clf()
figure(2)
clf()
Tos = [2,3,4,5]

for To in Tos:
    omegao = 2*pi/To
    omegas = linspace(0.1,5,100)
    s = []
    for omega in omegas:
        s.append(pm(omega,omegao))

    figure(1)
    plot(omegas,s,label='To=%f, Oo=%f'%(To,omegao))
    figure(2)
    plot(2*pi/omegas,s,label='To=%f, Oo=%f'%(To,omegao))

figure(1)
xlabel('omega [rad/s]')
ylabel('S(omega)')
legend()
grid(True)
figure(2)
xlabel('T [s]')
ylabel('S(omega)')
legend()
grid(True)
show()


for To in Tos:
    omegao = 2*pi/To
    omegas = [omegao/2.0, omegao, omegao*2]
    delo = [omegao/2.0, 3.0*omegao/4.0, omegao]
    a = []
    T = []
    for o,d in zip(omegas,delo):
        a.append(sqrt(2*pm(o,omegao)*d))
        T.append(2*pi/o)
    print To
    print T
    print a
    
    
