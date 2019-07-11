'''
Plot a few spectra
'''
# Wind velocity
U = 10.0


A = 0.00348
g = 9.81
L = U**2/g
# Phillips
k = linspace(0.01,1,100)
ph = (A/k**4)*exp(-1/(k*L)**2)


figure(1)
clf()
plot(k,ph,label='Ph')
xlabel('wave number [rad/m]')

w = linspace(0.01,5,100)
phw = A*pi*g**2/(w**5)*exp((-g**4)/((w*U)**4))

# Pierson-Moskowitz as a function of omega
alpha = 0.0081
beta = 0.74
pmw = alpha*g**2/(w**5)*exp(-beta*(g/(U*w))**4)


figure(2)
clf()
plot(w,phw,label='Ph')
plot(w,pmw,label='PM')
legend()
xlabel('omega [rad/s]')

H3 = 1 
gamma = 0.032
pmh = alpha*g**2/(w**5)*exp(-gamma*g**2/(H3**2*w**4))
U = sqrt((H3/sqrt(gamma))*g)

A = alpha/pi
phh = A*pi*g**2/(w**5)*exp((-g**4)/((w*U)**4))

figure(3)
clf()
plot(w,pmh,label='PM')
plot(w,phh,label='Ph')
legend()

show()
