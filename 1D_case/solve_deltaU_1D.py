"""
Semi-analytic solving of a hyperelastic bar impact (1D)

Analytical solution for vibrating parts (force without its constant value)
Assumption of small strains is made for vibrating part. Then, a linear PDE
is obtain. This equation is analytically solved in Fourier space
"""

#%% Librairies import
from scipy.fft import fft, ifft, fftfreq
from scipy.signal import find_peaks
from matplotlib.pyplot import plot, figure, xlabel, ylabel, grid, legend, show
from numpy import linspace, cos, sin, sqrt, pi, size, zeros, abs, ones

#%% External force definition
F0 = 0.01 #N
T = 0.1 #s
tf = 0.5 #s
n = 1000 #Number of points
t = linspace(0,tf,n)

#%% External force computation
f = zeros(size(t))

for i in range(len(t)):
    if t[i] <= T:
        f[i] = (F0/2)*(1 - cos(2*pi*t[i]/T))
##        f[i] = F0
    else:
        continue

#Plot force
figure(1)
plot(t,f,'-',label='External force')
xlabel('Time [s]')
ylabel('Force [N]')
grid(which='both')
legend(loc = 'best')

#%% Geometrical values
L = 0.1 #Beam length [m]
d0 = 0.01 #Beam diameter [m]
S0 = (pi/4)*d0**2 #Cross-sectional area [m²]

#%% Material parameters (from Kleiven, 2007)
mu1 = 107.6; mu2 = -240.8 #[Pa]
alpha1 = 10.1; alpha2 = -12.9 #[1]
mu = 0.5*(mu1*alpha1 + mu2*alpha2)
rho = 1.04e-3 #Brain density [kg/m^3]
# Deduce wave celerity
c = sqrt(3*mu/rho); print(c)
#%% Execute code to extract vibrating part of external force from FFT
with open("FFT_external_force_1D.py") as r:
    exec(r.read())

#%% Compute beam function H(X,omega) for different values of X
def H1d(X,omega):
    """
    This function computes "beam function" value for position X and pulsation omega

    INPUTS :
    - X [array or float]: position at initial configuration in m
    - omega [array or float]: pulsation in rad/s
    OUTPUT :
    - H [array or float]: "beam function" value in m/N
    """
    H = zeros(size(omega))
##    for i in range(len(omega)):
##        if abs(omega[i]) < 1e-6:
##            H[i] = (X - L)/(3*mu*S0)
##        else:
##            H[i] = (c/(3*mu*S0*omega[i]))*(sin(omega[i]*(X - L)/c)/cos(omega[i]*L/c))
    H[1:-1] = (c/(3*mu*S0*omega[1:-1]))*(sin(omega[1:-1]*(X - L)/c)/cos(omega[1:-1]*L/c))
    H[0] = (X - L)/(3*mu*S0)

    return H

# Compute for specific values of X (0, L/4, L/2, 3L/4, L)
omega = 2*pi*fftfreq(n, T)

H0 = H1d(0,omega) #X = 0
H1 = H1d(L/4,omega) #X = L/4
H2 = H1d(L/2,omega) #X = L/2
H3 = H1d(3*L/4,omega) #X = 3L/4
H4 = H1d(L,omega) #X = L
##for i in range(size(omega,1)):
##    if abs(omega[i]) < 1e-6:
##        H[i] = (0 - L)/(3*mu*S0)
##    else:
##        H[i] = (c/(3*mu*S0*omega[i]))*(sin(omega[i]*(0 - L)/c)/cos(omega[i]*L/c))


figure(3)
##plot(omega, 2.0/n * abs(dynamicFFreq[0:n//2]),'k-', label = 'Vibrating component')
plot(omega[:n//2],2.0/n*abs(H0[:n//2]),label = r"$H(0,\omega)$")
plot(omega[:n//2],2.0/n*abs(H1[:n//2]),label = r"$H(L/4,\omega)$")
plot(omega[:n//2],2.0/n*abs(H2[:n//2]),label = r"$H(L/2,\omega)$")
plot(omega[:n//2],2.0/n*abs(H3[:n//2]),label = r"$H(3L/4,\omega)$")
plot(omega[:n//2],2.0/n*abs(H4[:n//2]),label = r"$H(L,\omega)$")
xlabel(r'\omega')
legend(loc = 'best')

#%% Compute resultant deltaU spectra
U0 = H0*dynamicFFreq
U1 = H1*dynamicFFreq
U2 = H2*dynamicFFreq
U3 = H3*dynamicFFreq
U4 = H4*dynamicFFreq

figure(4)
plot(omega[:n//2],2.0/n*abs(U0[:n//2]),label = r"$\delta \tilde{u}(0,\omega)$")
plot(omega[:n//2],2.0/n*abs(U1[:n//2]),label = r"$\delta \tilde{u}(L/4,\omega)$")
plot(omega[:n//2],2.0/n*abs(U2[:n//2]),label = r"$\delta \tilde{u}(L/2,\omega)$")
plot(omega[:n//2],2.0/n*abs(U3[:n//2]),label = r"$\delta \tilde{u}(3L/4,\omega)$")
plot(omega[:n//2],2.0/n*abs(U4[:n//2]),label = r"$\delta \tilde{u}(L,\omega)$")
xlabel(r'\omega')
legend(loc = 'best')

#%% Then time-trajectories
u0 = ifft(U0); u0 = u0.real
u1 = ifft(U1); u1 = u1.real
u2 = ifft(U2); u2 = u2.real
u3 = ifft(U3); u3 = u3.real
u4 = ifft(U4); u4 = u4.real

figure(5)
plot(t,u0,label = r"$\delta u(0,t)$")
plot(t,u1,label = r"$\delta u(L/4,t)$")
plot(t,u2,label = r"$\delta u(L/2,t)$")
plot(t,u3,label = r"$\delta u(3L/4,t)$")
plot(t,u4,label = r"$\delta u(L,t)$")
xlabel('time')


show()
