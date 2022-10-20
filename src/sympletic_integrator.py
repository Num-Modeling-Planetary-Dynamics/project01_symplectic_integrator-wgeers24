%matplotlib inline
from orbitgeometry import *
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib im
import astropy as astro


#Initial Conditions

#Global constants in SI units
G= 6.6743E-11
m_Sun=1.9885E+30
mu= -G*(m_Sun)
Gmass=(m_Neptune + m_Pluto)

JD2S=86400
YR2S=np.longdouble(365.25*JD2S)
#Neptune Initial Conditions
m_Neptune= 102.409E+24
x_Neptune=2.266223889656152E+08
y_Neptune=4.461841849742382E+09
z_Neptune=-9.710614902790475E+07
vx_Neptune=-5.447617894349518E+00
vy_Neptune= 3.044797889906726E-01
vz_Neptune= 1.187834180290304E-01

#Pluto Initial Conditions
m_Pluto= 1.307E+22
x_Pluto= 1.540431633806554E+09
y_Pluto= 6.753704837972438E+09
z_Pluto=-1.168941164341141E+09
vx_Pluto=-3.741653796726164E+00
vy_Pluto= 3.621999156857085E-01
vz_Pluto= 1.046185526637238E+00


def xv2el (mu,x,y,z,vx,vy,vz):
#input position and velocity vectors in cartesian coordinates
    r_vec=np.array([x,y,z])
    r_mag=np.linalg.norm(r_vec)
    v_vec=np.array([vx,vy,vz])
    v_mag=np.linalg.norm(v_vec)
    h_vec=np.cross(r,v)
    h_mag=np.linalg.norm(h_vec)
    r_r_dot=np.vdot(r_vec,v_vec)
    r_dot=np.sign(r_r_dot)*np.sqrt((v_mag**2)-(h_mag**2/r_mag**2))


    a=((2/r_mag)-(v_mag**2/mu))**-1
    ecc=np.sqrt(1-(h_mag**2/(mu*a)))
    I=np.arccos(hz/h_mag)

    sO=(np.sign(h[2])*h[0]/(h_mag*np.sin(I)))
    cO=(np.sign(h[2])*h[1])/(h_mag*np.sin(I)))
    omega=np.arctan2(sO,cO)

    sof=z/(r_mag*np.sin(I))
    cof=((x/r_mag)+np.sin(omega)*sof*np.cos(I))/np.cos(omega)
    pf=np.arctan2(sof,cof)

    sf=((a*(1-ecc**2))/(h_mag*ecc))*(r_dot)
    cf=(1/ecc)*((a*(1-e**2)/r_mag)-1)
    f=np.arctan2(sf,cf)
    peri=pf-f



#output orbital elements
#Use equations from sxn 2.8 pg 52 and 53
    return a,ecc,I,omega,peri,f

def vh2vb (vhvec, mu, Gmass):
    #Converts heliocentric velocity vector to barycentric
    Gmtot= mu + np.sum(Gmass)
    vbcb = -np.sum(Gmass* vhvec,axis=0)/Gmtot
    vbvec= vhvec + vbcb
    #vhvec: heliocentric velocity vectors of all n bodies
    #Gmass: masses of n bodies in system
    #vb: barycentric velocity vectors of n bodies
    #mu: central body gravitational parameter
    #vbcv: barycentric velocity vector of central body
    return vbvec, vbcb

def vb2vh (vbvec, mu, Gmass):
    Gmtot= mu + np.sum(Gmass)
    vbcb = -np.sum(Gmass * vbvec, axis=0)/ mu
    vhvec= vbvec - vbcb
    # vhvec: heliocentric velocity vectors of all n bodies
    # Gmass: masses of n bodies in system
    # vb: barycentric velocity vectors of n bodies
    # mu: central body gravitational parameter
    # vbcv: barycentric velocity vector of central body
    return vhvec



def Kep_drift (,dt):
    #Uses Kepler's equation for eccentric anomaly E

    E = [M]
    e=ecc
    # Equation 2.49 solve for E
    def Danby(M,ecc):

        def f(E):
            return E - e * np.sin(E) - M

        def fP(E):
            return 1 - e * np.cos(E)

        def fP2(E):
            return e * np.sin(E)

        def fP3(E):
            return e * np.cos(E)

        def D1(E):
            return -(f(E)) / (fP(E))

        def D2(E):
            return -(f(E)) / (fP(E) + 0.5 * (D1(E) * fP2(E)))

        def D3(E):
            return -(f(E)) / (fP(E) + 0.5 * (D2(E) * fP2(E)) + (1 / 6) * ((D2(E) ** 2) * fP3(E)))

        while True:
            Enew = E[-1] + D3(E[-1])
            E.append(Enew)
            if np.abs((E[-1] - E[-2])/E[-2])<accuracy:
                break
    #utilizes f and g functions
    #Danby method
    if ecc < np.finfo(np.float64).tiny:
        E0=0.0
    else:
        E0=np.arccos(-(rmag0-a)/(a*ecc))

    if np.sign(np.vdot(r0,v0))<0.0
        E0= 2*np.pi-E0

    M0=E0-ecc*np.sin(E0)

    M-M0+n*dt
    E=danby(M,ecc)


    f=a/rmag0*(np.cos(dE)-1.0)+1.0
    g=dt+1.0/n*(np.sin(dE)-dE)

    r_mag=f*r_mag0+g*v_mag0
    return r_mag
    #Define fdot and gdot for vmag using eq 2.71 on pg 37
def kick ():
    #E_int step
    #uses the Hamiltonian for the interaction step
    #half-time step (solar drift/linear drift)


def Sun_drift(self,dt):
    #updates barycenter of Sun
    pt = np.sum(self.Gmass * self.vbvec, axis=0) / self.param['GMcb']
    self.rhvec += pt*dt

def step(dt):
    #Advances simulation by dt

    dth=0.5*dt

def drift_one(mu,x,y,z,vx,vy,vz,dt):
    #inputs cartesian coordinates

    #Outputs new position and velocity cartesian coordinates

    def calc_energy_one(vhvec, rhvec, Gmass, mu):
        # Computes system energy at particular state of simulation
        # rhvec: heliocentric position vectors of n bodies
        # vhvec: heliocentric velocity vectors of n bodies
        # Gmass: masses of n bodies in system
        # mu: central body gravitational parameter
        # vbvec,vbcb: from vh2vb
        vbmag2 = np.einsum("ij,ij->i", vbvec, vbvec)
        irh = 1.0 / np.linalg.norm(rhvec, axis=1)
        ke = 0.5 * (mu * np.vdot(vbcb, vbcb) + np.sum(Gmass * vbmag2)) / GC

        def pe_one(Gm, rvec, i):
            drvec = r_vec[i + 1:, :] - rvec[i, :]
            irij = np.linalg.norm(drvec, axis=1)
            irij = Gm[i + 1:] * Gm[i] / irij
            return np.sum(drvec.T * irij, axis=1)

        n = rhvec.shape[0]
        pe = (-mu * np.sum(Gmass * irh) - np.sum([pe_one(Gmass.flatten(), rhvec, i) for i in range(n - 1)])) / GC

        return ke + pe

